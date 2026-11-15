//! Chunked, compressed, checksummed columnar trajectory format (version 1).
//!
//! A trajectory is a sequence of [`Frame`]s (step, time, orthorhombic box and
//! positions for a fixed number of atoms). Frames are grouped into chunks of
//! `chunk_frames` frames. Inside a chunk every quantity is stored as its own
//! column over the frames of that chunk, so similar values sit next to each
//! other and compress well. The index is written after the last chunk, which
//! lets a writer stream frames without seeking and lets a reader jump to any
//! frame.
//!
//! All integers are little-endian.
//!
//! ```text
//! file header (28 bytes)
//!   magic       8 bytes   b"TPTCTRJ1"
//!   version     u16       1
//!   precision   u8        1 = f32 positions/box/time, 2 = f64
//!   reserved    u8        0 (non-zero is rejected)
//!   n_atoms     u64
//!   chunk_frames u32      frames per full chunk (>= 1)
//!   header_crc  u32       CRC-32 of the 24 bytes before it
//!
//! chunk record, repeated (one per chunk, in frame order)
//!   tag         4 bytes   b"CHNK"
//!   n_frames    u32       frames in this chunk (>= 1, <= chunk_frames)
//!   payload_len u64       bytes of payload that follow
//!   payload     payload_len bytes (layout below)
//!   crc         u32       CRC-32 over n_frames, payload_len and payload
//!
//! index block (after the last chunk record)
//!   tag         4 bytes   b"INDX"
//!   n_chunks    u64
//!   entries     n_chunks × { record_offset u64, record_len u64,
//!                            first_frame u64, n_frames u32 }
//!   crc         u32       CRC-32 over tag, n_chunks and entries
//!
//! footer (24 bytes, always the last bytes of the file)
//!   index_offset u64      byte offset of the index block
//!   total_frames u64
//!   end_magic    8 bytes  b"TPTCEND1"
//! ```
//!
//! # Chunk payload
//!
//! The payload holds `5 + 3·n_atoms` columns in this fixed order: `step`,
//! `time`, `box_x`, `box_y`, `box_z`, then `x, y, z` of atom 0, of atom 1, and so
//! on. The payload starts with a table of `ncols` `u32` byte lengths, followed
//! by the column bytes back to back. The table lets a reader decode only the
//! columns it needs.
//!
//! # Compression (lossless step)
//!
//! Each column is encoded as follows, with no external dependency:
//!
//! 1. Every value is mapped to an unsigned *key*. For `step` the key is the
//!    integer itself. For floating-point values the key is the IEEE bit pattern
//!    with the sign handled so that numeric order equals key order (negative
//!    values have all bits flipped, non-negative values have only the sign bit
//!    set). Adjacent floats therefore have adjacent keys, and this mapping is a
//!    bijection, so NaN payloads and `-0.0` survive.
//! 2. Consecutive keys are differenced (wrapping, starting from zero) and the
//!    differences are zig-zag encoded so small negative steps are small.
//! 3. The results are written as LEB128 varints.
//!
//! For a smooth trajectory, consecutive frames share high-order key bits, so the
//! deltas are short. Constant columns (for example a fixed box) cost one byte
//! per frame.
//!
//! # Precision
//!
//! With [`Precision::F32`] the writer rounds positions, box lengths and time to
//! single precision once, when a frame is pushed. Everything after that is
//! lossless: the reader returns exactly the stored `f32` values, widened to
//! `f64`. With [`Precision::F64`] values are stored exactly.
//!
//! # Random access and bounded memory
//!
//! [`TrajectoryReader::open`] reads the header, footer and index only. Reading
//! frame `i` binary-searches the index, then reads and verifies one chunk
//! record. [`TrajectoryReader::frames`] streams frames in order and keeps one
//! decoded chunk in memory at a time, so memory use depends on `chunk_frames`
//! and `n_atoms`, not on the trajectory length.
//!
//! [`LazyTrajectory`] exposes the same data as a lazy-engine [`Source`]. Each
//! column is named and projection-pushed into the reader, so unused columns are
//! never decoded.
//!
//! # Robustness
//!
//! Every header, record, index and footer field is covered by a checksum or a
//! consistency check. Allocations are bounded by bytes actually present or by
//! values validated against the file length. Corrupt input returns a
//! [`TrajectoryError`]; nothing panics.
//!
//! ```
//! use std::io::Cursor;
//! use tpt_dynamics_data::trajectory::{Frame, Precision, TrajectoryReader, TrajectoryWriter, FrameSink};
//!
//! let mut w = TrajectoryWriter::new(Vec::new(), 2, Precision::F64, 4).unwrap();
//! for step in 0..10u64 {
//!     let f = Frame {
//!         step,
//!         time: 0.002 * step as f64,
//!         box_lengths: [3.0, 3.0, 3.0],
//!         positions: vec![[0.1 * step as f64, 0.0, 0.0], [1.0, 1.0, 1.0]],
//!     };
//!     w.push_frame(&f).unwrap();
//! }
//! let bytes = w.finish().unwrap();
//! let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
//! assert_eq!(r.n_frames(), 10);
//! assert_eq!(r.frame(7).unwrap().step, 7);
//! ```

use std::fmt;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::sync::Mutex;

use crate::crc32::{checksum, Crc32};
use crate::lazy::{Column, DType as LazyDType, DataFrame, Field, LazyError, Schema, Source};

/// File magic.
pub const MAGIC: &[u8; 8] = b"TPTCTRJ1";
/// Footer end magic.
pub const END_MAGIC: &[u8; 8] = b"TPTCEND1";
/// Format version written by this crate.
pub const VERSION: u16 = 1;

const CHUNK_TAG: &[u8; 4] = b"CHNK";
const INDEX_TAG: &[u8; 4] = b"INDX";
const HEADER_LEN: u64 = 28;
const RECORD_HEADER_LEN: u64 = 16;
const RECORD_TRAILER_LEN: u64 = 4;
const INDEX_HEAD_LEN: u64 = 12;
const INDEX_ENTRY_LEN: u64 = 28;
const FOOTER_LEN: u64 = 24;

/// Columns before the coordinates: `step`, `time`, `box_x`, `box_y`, `box_z`.
const N_LEADING_COLUMNS: usize = 5;

/// Numeric precision used for floating-point quantities in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Precision {
    /// Single precision (values are rounded once, at push time).
    F32,
    /// Double precision (exact).
    #[default]
    F64,
}

impl Precision {
    fn code(self) -> u8 {
        match self {
            Self::F32 => 1,
            Self::F64 => 2,
        }
    }

    fn from_code(c: u8) -> Option<Self> {
        match c {
            1 => Some(Self::F32),
            2 => Some(Self::F64),
            _ => None,
        }
    }

    /// Maps a value to its ordered key (see the module docs).
    fn key(self, v: f64) -> u64 {
        match self {
            Self::F32 => {
                let b = (v as f32).to_bits();
                u64::from(if b & SIGN32 != 0 { !b } else { b | SIGN32 })
            }
            Self::F64 => {
                let b = v.to_bits();
                if b & SIGN64 != 0 {
                    !b
                } else {
                    b | SIGN64
                }
            }
        }
    }

    /// Inverse of [`Precision::key`]. Keys above the type width are not valid
    /// for this precision; the reader rejects them before calling this.
    fn value(self, k: u64) -> f64 {
        match self {
            Self::F32 => {
                let e = k as u32;
                let b = if e & SIGN32 != 0 { e & !SIGN32 } else { !e };
                f64::from(f32::from_bits(b))
            }
            Self::F64 => {
                let b = if k & SIGN64 != 0 { k & !SIGN64 } else { !k };
                f64::from_bits(b)
            }
        }
    }

    fn max_key(self) -> u64 {
        match self {
            Self::F32 => u64::from(u32::MAX),
            Self::F64 => u64::MAX,
        }
    }
}

const SIGN32: u32 = 0x8000_0000;
const SIGN64: u64 = 0x8000_0000_0000_0000;

/// One trajectory frame.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Frame {
    /// Integer step index.
    pub step: u64,
    /// Simulation time in ps.
    pub time: f64,
    /// Orthorhombic box edge lengths in nm.
    pub box_lengths: [f64; 3],
    /// Positions in nm, one `[x, y, z]` per atom.
    pub positions: Vec<[f64; 3]>,
}

/// Errors from writing or reading the trajectory format.
#[derive(Debug)]
pub enum TrajectoryError {
    /// Underlying I/O failure.
    Io(io::Error),
    /// The input ended before the expected data.
    Truncated,
    /// The file does not start with [`MAGIC`].
    BadMagic,
    /// The footer does not end with [`END_MAGIC`].
    BadEndMagic,
    /// The format version is not supported by this build.
    UnsupportedVersion(u16),
    /// A header field has an invalid value.
    BadHeader(&'static str),
    /// A checksum did not match. `what` names the failing part.
    ChecksumMismatch {
        /// Which part failed verification.
        what: String,
    },
    /// Structure is inconsistent (bad index, bad counts, bad varint, ...).
    Corrupt(String),
    /// A frame has a different number of atoms than the file.
    AtomCountMismatch {
        /// Atoms per frame in the file.
        expected: usize,
        /// Atoms in the offered frame.
        found: usize,
    },
    /// A frame index is past the end of the trajectory.
    FrameOutOfRange {
        /// Requested frame.
        index: u64,
        /// Frames in the trajectory.
        n_frames: u64,
    },
    /// Frames were pushed after the writer was closed.
    Closed,
    /// A writer argument was invalid.
    InvalidInput(&'static str),
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Truncated => write!(f, "trajectory data is truncated"),
            Self::BadMagic => write!(f, "not a tpt-dynamics compressed trajectory (bad magic)"),
            Self::BadEndMagic => write!(f, "trajectory footer is missing (file incomplete?)"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported trajectory format version {v}"),
            Self::BadHeader(m) => write!(f, "invalid trajectory header: {m}"),
            Self::ChecksumMismatch { what } => write!(f, "checksum mismatch in {what}"),
            Self::Corrupt(m) => write!(f, "corrupt trajectory: {m}"),
            Self::AtomCountMismatch { expected, found } => {
                write!(f, "frame has {found} atoms, trajectory expects {expected}")
            }
            Self::FrameOutOfRange { index, n_frames } => {
                write!(f, "frame {index} out of range (trajectory has {n_frames})")
            }
            Self::Closed => write!(f, "trajectory writer is already closed"),
            Self::InvalidInput(m) => write!(f, "invalid trajectory argument: {m}"),
        }
    }
}

impl std::error::Error for TrajectoryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for TrajectoryError {
    fn from(e: io::Error) -> Self {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            Self::Truncated
        } else {
            Self::Io(e)
        }
    }
}

/// Destination for frames produced by a running simulation or a converter.
///
/// The simulation's run loop can call [`FrameSink::push_frame`] at each output
/// step and [`FrameSink::close`] at the end of the run.
pub trait FrameSink {
    /// Appends one frame.
    ///
    /// # Errors
    /// Returns an error on I/O failure, on an atom-count mismatch, or if the
    /// sink is closed.
    fn push_frame(&mut self, frame: &Frame) -> Result<(), TrajectoryError>;

    /// Flushes buffered frames and writes the index. Calling it again is a no-op.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    fn close(&mut self) -> Result<(), TrajectoryError>;
}

// ---------------------------------------------------------------------------
// Varint, zig-zag and column codecs
// ---------------------------------------------------------------------------

fn zigzag(d: i64) -> u64 {
    ((d << 1) ^ (d >> 63)) as u64
}

fn unzigzag(z: u64) -> i64 {
    (z >> 1).cast_signed() ^ -(z & 1).cast_signed()
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let b = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(b);
            return;
        }
        out.push(b | 0x80);
    }
}

fn get_varint(buf: &[u8], pos: &mut usize) -> Result<u64, TrajectoryError> {
    let mut v = 0u64;
    let mut shift = 0u32;
    loop {
        let b = *buf
            .get(*pos)
            .ok_or_else(|| TrajectoryError::Corrupt("varint runs past end of column".into()))?;
        *pos += 1;
        // The tenth byte may only carry bit 63; anything else overflows u64.
        if shift == 63 && b > 1 {
            return Err(TrajectoryError::Corrupt("varint overflows u64".into()));
        }
        v |= u64::from(b & 0x7f) << shift;
        if b & 0x80 == 0 {
            return Ok(v);
        }
        shift += 7;
    }
}

fn encode_column(keys: &[u64], out: &mut Vec<u8>) {
    let mut prev = 0u64;
    for &k in keys {
        put_varint(out, zigzag(k.wrapping_sub(prev).cast_signed()));
        prev = k;
    }
}

fn decode_column(buf: &[u8], n: usize, max_key: u64) -> Result<Vec<u64>, TrajectoryError> {
    let mut pos = 0usize;
    let mut prev = 0u64;
    let mut out = Vec::with_capacity(n.min(buf.len()));
    for _ in 0..n {
        let d = unzigzag(get_varint(buf, &mut pos)?);
        prev = prev.wrapping_add(d as u64);
        if prev > max_key {
            return Err(TrajectoryError::Corrupt("key exceeds precision width".into()));
        }
        out.push(prev);
    }
    if pos != buf.len() {
        return Err(TrajectoryError::Corrupt("column has trailing bytes".into()));
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Chunk records
// ---------------------------------------------------------------------------

fn ncols(n_atoms: usize) -> usize {
    N_LEADING_COLUMNS + 3 * n_atoms
}

/// Encodes frames into a chunk payload (length table followed by column bytes).
fn encode_payload(frames: &[Frame], n_atoms: usize, precision: Precision) -> Vec<u8> {
    let mut columns: Vec<Vec<u8>> = Vec::with_capacity(ncols(n_atoms));
    let mut push_col = |keys: &[u64]| {
        let mut bytes = Vec::with_capacity(keys.len());
        encode_column(keys, &mut bytes);
        columns.push(bytes);
    };
    let steps: Vec<u64> = frames.iter().map(|f| f.step).collect();
    push_col(&steps);
    let times: Vec<u64> = frames.iter().map(|f| precision.key(f.time)).collect();
    push_col(&times);
    for d in 0..3 {
        let keys: Vec<u64> = frames.iter().map(|f| precision.key(f.box_lengths[d])).collect();
        push_col(&keys);
    }
    for a in 0..n_atoms {
        for c in 0..3 {
            let keys: Vec<u64> = frames.iter().map(|f| precision.key(f.positions[a][c])).collect();
            push_col(&keys);
        }
    }

    let table_len = 4 * columns.len();
    let body: usize = columns.iter().map(Vec::len).sum();
    let mut payload = Vec::with_capacity(table_len + body);
    for col in &columns {
        payload.extend_from_slice(&(col.len() as u32).to_le_bytes());
    }
    for col in &columns {
        payload.extend_from_slice(col);
    }
    payload
}

/// A parsed chunk payload: per-column byte ranges into the payload.
struct PayloadView<'a> {
    payload: &'a [u8],
    ranges: Vec<(usize, usize)>,
}

impl<'a> PayloadView<'a> {
    fn parse(payload: &'a [u8], n_cols: usize) -> Result<Self, TrajectoryError> {
        let table_len = 4 * n_cols;
        if payload.len() < table_len {
            return Err(TrajectoryError::Corrupt("chunk payload shorter than its table".into()));
        }
        let mut ranges = Vec::with_capacity(n_cols);
        let mut offset = table_len;
        for i in 0..n_cols {
            let mut b = [0u8; 4];
            b.copy_from_slice(&payload[4 * i..4 * i + 4]);
            let len = u32::from_le_bytes(b) as usize;
            let end = offset
                .checked_add(len)
                .filter(|&e| e <= payload.len())
                .ok_or_else(|| TrajectoryError::Corrupt("column length past payload".into()))?;
            ranges.push((offset, end));
            offset = end;
        }
        if offset != payload.len() {
            return Err(TrajectoryError::Corrupt("payload has unaccounted bytes".into()));
        }
        Ok(Self { payload, ranges })
    }

    fn column(&self, i: usize) -> &'a [u8] {
        let (a, b) = self.ranges[i];
        &self.payload[a..b]
    }
}

fn write_record<W: Write>(
    w: &mut W,
    n_frames: u32,
    payload: &[u8],
) -> Result<u64, TrajectoryError> {
    let payload_len = payload.len() as u64;
    let n_bytes = n_frames.to_le_bytes();
    let len_bytes = payload_len.to_le_bytes();
    let mut crc = Crc32::new();
    crc.update(&n_bytes);
    crc.update(&len_bytes);
    crc.update(payload);
    w.write_all(CHUNK_TAG)?;
    w.write_all(&n_bytes)?;
    w.write_all(&len_bytes)?;
    w.write_all(payload)?;
    w.write_all(&crc.finish().to_le_bytes())?;
    Ok(RECORD_HEADER_LEN + payload_len + RECORD_TRAILER_LEN)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct ChunkEntry {
    record_offset: u64,
    record_len: u64,
    first_frame: u64,
    n_frames: u32,
}

/// Streaming writer. Frames are buffered until a chunk is full, then encoded
/// and written. Nothing needs to be seeked, so any [`Write`] works.
#[derive(Debug)]
pub struct TrajectoryWriter<W: Write> {
    inner: W,
    n_atoms: usize,
    precision: Precision,
    chunk_frames: usize,
    buffer: Vec<Frame>,
    index: Vec<ChunkEntry>,
    offset: u64,
    frames_pushed: u64,
    closed: bool,
}

impl<W: Write> TrajectoryWriter<W> {
    /// Writes the header and returns a writer for `n_atoms` atoms per frame,
    /// with `chunk_frames` frames per chunk.
    ///
    /// # Errors
    /// Returns [`TrajectoryError::InvalidInput`] if `chunk_frames` is zero, or
    /// an I/O error if the header cannot be written.
    pub fn new(
        mut inner: W,
        n_atoms: usize,
        precision: Precision,
        chunk_frames: u32,
    ) -> Result<Self, TrajectoryError> {
        if chunk_frames == 0 {
            return Err(TrajectoryError::InvalidInput("chunk_frames must be at least 1"));
        }
        let mut header = Vec::with_capacity(HEADER_LEN as usize);
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&VERSION.to_le_bytes());
        header.push(precision.code());
        header.push(0);
        header.extend_from_slice(&(n_atoms as u64).to_le_bytes());
        header.extend_from_slice(&chunk_frames.to_le_bytes());
        let crc = checksum(&header);
        header.extend_from_slice(&crc.to_le_bytes());
        inner.write_all(&header)?;
        Ok(Self {
            inner,
            n_atoms,
            precision,
            chunk_frames: chunk_frames as usize,
            buffer: Vec::with_capacity(chunk_frames.min(4096) as usize),
            index: Vec::new(),
            offset: HEADER_LEN,
            frames_pushed: 0,
            closed: false,
        })
    }

    /// Frames accepted so far (including ones still buffered).
    pub fn frames_pushed(&self) -> u64 {
        self.frames_pushed
    }

    /// Closes the trajectory and returns the underlying writer.
    ///
    /// # Errors
    /// Returns an error if the final chunk or the index cannot be written.
    pub fn finish(mut self) -> Result<W, TrajectoryError> {
        self.close()?;
        Ok(self.inner)
    }

    fn flush_chunk(&mut self) -> Result<(), TrajectoryError> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let n = self.buffer.len() as u32;
        let payload = encode_payload(&self.buffer, self.n_atoms, self.precision);
        let record_len = write_record(&mut self.inner, n, &payload)?;
        self.index.push(ChunkEntry {
            record_offset: self.offset,
            record_len,
            first_frame: self.frames_pushed - u64::from(n),
            n_frames: n,
        });
        self.offset += record_len;
        self.buffer.clear();
        Ok(())
    }
}

impl<W: Write> FrameSink for TrajectoryWriter<W> {
    fn push_frame(&mut self, frame: &Frame) -> Result<(), TrajectoryError> {
        if self.closed {
            return Err(TrajectoryError::Closed);
        }
        if frame.positions.len() != self.n_atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.n_atoms,
                found: frame.positions.len(),
            });
        }
        let frame = match self.precision {
            Precision::F64 => frame.clone(),
            Precision::F32 => Frame {
                step: frame.step,
                time: f64::from(frame.time as f32),
                box_lengths: frame.box_lengths.map(|v| f64::from(v as f32)),
                positions: frame.positions.iter().map(|p| p.map(|v| f64::from(v as f32))).collect(),
            },
        };
        self.buffer.push(frame);
        self.frames_pushed += 1;
        if self.buffer.len() == self.chunk_frames {
            self.flush_chunk()?;
        }
        Ok(())
    }

    fn close(&mut self) -> Result<(), TrajectoryError> {
        if self.closed {
            return Ok(());
        }
        self.flush_chunk()?;
        let index_offset = self.offset;
        let mut block = Vec::with_capacity(
            (INDEX_HEAD_LEN + INDEX_ENTRY_LEN * self.index.len() as u64) as usize,
        );
        block.extend_from_slice(INDEX_TAG);
        block.extend_from_slice(&(self.index.len() as u64).to_le_bytes());
        for e in &self.index {
            block.extend_from_slice(&e.record_offset.to_le_bytes());
            block.extend_from_slice(&e.record_len.to_le_bytes());
            block.extend_from_slice(&e.first_frame.to_le_bytes());
            block.extend_from_slice(&e.n_frames.to_le_bytes());
        }
        let crc = checksum(&block);
        block.extend_from_slice(&crc.to_le_bytes());
        self.inner.write_all(&block)?;
        let mut footer = Vec::with_capacity(FOOTER_LEN as usize);
        footer.extend_from_slice(&index_offset.to_le_bytes());
        footer.extend_from_slice(&self.frames_pushed.to_le_bytes());
        footer.extend_from_slice(END_MAGIC);
        self.inner.write_all(&footer)?;
        self.inner.flush()?;
        self.closed = true;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// Header fields of an open trajectory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Header {
    /// Atoms per frame.
    pub n_atoms: usize,
    /// Precision of stored floating-point values.
    pub precision: Precision,
    /// Frames per full chunk.
    pub chunk_frames: u32,
}

/// Index entry for one chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChunkInfo {
    /// Index of the first frame in this chunk.
    pub first_frame: u64,
    /// Frames in this chunk.
    pub n_frames: u32,
    /// Byte offset of the chunk record in the file.
    pub record_offset: u64,
    /// Total bytes of the chunk record (header, payload and checksum).
    pub record_len: u64,
}

/// Counters for chunk reads, used to verify memory and I/O bounds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReadStats {
    /// Chunk records read and decoded.
    pub chunks_loaded: u64,
    /// Bytes of chunk records read.
    pub bytes_read: u64,
}

/// Random-access and streaming reader for the compressed trajectory format.
#[derive(Debug)]
pub struct TrajectoryReader<R: Read + Seek> {
    inner: R,
    header: Header,
    chunks: Vec<ChunkInfo>,
    n_frames: u64,
    stats: ReadStats,
}

impl<R: Read + Seek> TrajectoryReader<R> {
    /// Opens a trajectory by reading its header, footer and index. Chunk data
    /// is not read.
    ///
    /// # Errors
    /// Returns an error for a bad magic, unsupported version, checksum failure,
    /// or inconsistent index.
    pub fn open(mut inner: R) -> Result<Self, TrajectoryError> {
        let file_len = inner.seek(SeekFrom::End(0))?;
        if file_len < HEADER_LEN + INDEX_HEAD_LEN + RECORD_TRAILER_LEN + FOOTER_LEN {
            return Err(TrajectoryError::Truncated);
        }

        inner.seek(SeekFrom::Start(0))?;
        let mut hdr = [0u8; HEADER_LEN as usize];
        inner.read_exact(&mut hdr)?;
        if &hdr[..8] != MAGIC {
            return Err(TrajectoryError::BadMagic);
        }
        let version = u16::from_le_bytes([hdr[8], hdr[9]]);
        if version != VERSION {
            return Err(TrajectoryError::UnsupportedVersion(version));
        }
        if checksum(&hdr[..24]) != read_u32(&hdr[24..28]) {
            return Err(TrajectoryError::ChecksumMismatch { what: "header".into() });
        }
        let precision = Precision::from_code(hdr[10])
            .ok_or(TrajectoryError::BadHeader("unknown precision code"))?;
        if hdr[11] != 0 {
            return Err(TrajectoryError::BadHeader("reserved byte is non-zero"));
        }
        let n_atoms = usize::try_from(read_u64(&hdr[12..20]))
            .map_err(|_| TrajectoryError::BadHeader("atom count does not fit in usize"))?;
        let chunk_frames = read_u32(&hdr[20..24]);
        if chunk_frames == 0 {
            return Err(TrajectoryError::BadHeader("chunk_frames is zero"));
        }

        inner.seek(SeekFrom::End(-FOOTER_LEN.cast_signed()))?;
        let mut footer = [0u8; FOOTER_LEN as usize];
        inner.read_exact(&mut footer)?;
        if &footer[16..24] != END_MAGIC {
            return Err(TrajectoryError::BadEndMagic);
        }
        let index_offset = read_u64(&footer[0..8]);
        let total_frames = read_u64(&footer[8..16]);

        let index_end = file_len - FOOTER_LEN;
        if index_offset < HEADER_LEN || index_offset + INDEX_HEAD_LEN + 4 > index_end {
            return Err(TrajectoryError::Corrupt("index offset out of range".into()));
        }
        inner.seek(SeekFrom::Start(index_offset))?;
        let mut head = [0u8; INDEX_HEAD_LEN as usize];
        inner.read_exact(&mut head)?;
        if &head[..4] != INDEX_TAG {
            return Err(TrajectoryError::Corrupt("index tag missing".into()));
        }
        let n_chunks = read_u64(&head[4..12]);
        let body_len = n_chunks
            .checked_mul(INDEX_ENTRY_LEN)
            .ok_or_else(|| TrajectoryError::Corrupt("chunk count overflows".into()))?;
        if index_offset + INDEX_HEAD_LEN + body_len + 4 != index_end {
            return Err(TrajectoryError::Corrupt("index length does not match footer".into()));
        }
        let mut block = Vec::with_capacity((INDEX_HEAD_LEN + body_len) as usize);
        block.extend_from_slice(&head);
        let body_start = block.len();
        block.resize(body_start + body_len as usize, 0);
        inner.read_exact(&mut block[body_start..])?;
        let mut crc_bytes = [0u8; 4];
        inner.read_exact(&mut crc_bytes)?;
        if checksum(&block) != u32::from_le_bytes(crc_bytes) {
            return Err(TrajectoryError::ChecksumMismatch { what: "index".into() });
        }

        let mut chunks = Vec::with_capacity(n_chunks as usize);
        let mut expected_first = 0u64;
        for e in block[INDEX_HEAD_LEN as usize..].chunks(INDEX_ENTRY_LEN as usize) {
            let info = ChunkInfo {
                record_offset: read_u64(&e[0..8]),
                record_len: read_u64(&e[8..16]),
                first_frame: read_u64(&e[16..24]),
                n_frames: read_u32(&e[24..28]),
            };
            if info.n_frames == 0
                || info.n_frames > chunk_frames
                || info.first_frame != expected_first
                || info.record_len < RECORD_HEADER_LEN + RECORD_TRAILER_LEN
                || info.record_offset < HEADER_LEN
                || info.record_offset.checked_add(info.record_len).is_none_or(|e| e > index_offset)
            {
                return Err(TrajectoryError::Corrupt("inconsistent chunk index entry".into()));
            }
            expected_first += u64::from(info.n_frames);
            chunks.push(info);
        }
        if expected_first != total_frames {
            return Err(TrajectoryError::Corrupt("frame count disagrees with index".into()));
        }

        Ok(Self {
            inner,
            header: Header { n_atoms, precision, chunk_frames },
            chunks,
            n_frames: total_frames,
            stats: ReadStats::default(),
        })
    }

    /// Header fields.
    pub fn header(&self) -> Header {
        self.header
    }

    /// Atoms per frame.
    pub fn n_atoms(&self) -> usize {
        self.header.n_atoms
    }

    /// Number of frames.
    pub fn n_frames(&self) -> u64 {
        self.n_frames
    }

    /// The chunk index, in frame order.
    pub fn chunks(&self) -> &[ChunkInfo] {
        &self.chunks
    }

    /// Chunk-read counters since the reader was opened or last reset.
    pub fn stats(&self) -> ReadStats {
        self.stats
    }

    /// Resets the chunk-read counters to zero.
    pub fn reset_stats(&mut self) {
        self.stats = ReadStats::default();
    }

    /// Reads and verifies the record of chunk `chunk` and returns its frame
    /// count and payload.
    fn read_payload(&mut self, chunk: usize) -> Result<(u32, Vec<u8>), TrajectoryError> {
        let info = self.chunks[chunk];
        self.inner.seek(SeekFrom::Start(info.record_offset))?;
        let mut head = [0u8; RECORD_HEADER_LEN as usize];
        self.inner.read_exact(&mut head)?;
        if &head[..4] != CHUNK_TAG {
            return Err(TrajectoryError::Corrupt(format!("chunk {chunk}: bad record tag")));
        }
        let n_bytes = &head[4..8];
        let len_bytes = &head[8..16];
        let n_frames = read_u32(n_bytes);
        let payload_len = read_u64(len_bytes);
        if n_frames != info.n_frames
            || payload_len.checked_add(RECORD_HEADER_LEN + RECORD_TRAILER_LEN)
                != Some(info.record_len)
        {
            return Err(TrajectoryError::Corrupt(format!(
                "chunk {chunk}: header disagrees with index"
            )));
        }
        // Read through `take` so a bogus length cannot allocate more than the file holds.
        let mut payload = Vec::new();
        (&mut self.inner).take(payload_len).read_to_end(&mut payload)?;
        if payload.len() as u64 != payload_len {
            return Err(TrajectoryError::Truncated);
        }
        let mut crc_bytes = [0u8; 4];
        self.inner.read_exact(&mut crc_bytes)?;
        let mut crc = Crc32::new();
        crc.update(n_bytes);
        crc.update(len_bytes);
        crc.update(&payload);
        if crc.finish() != u32::from_le_bytes(crc_bytes) {
            return Err(TrajectoryError::ChecksumMismatch { what: format!("chunk {chunk}") });
        }
        self.stats.chunks_loaded += 1;
        self.stats.bytes_read += info.record_len;
        Ok((n_frames, payload))
    }

    /// Decodes the listed columns of a chunk. Returns the keys per requested
    /// column, in the order given.
    fn read_columns(
        &mut self,
        chunk: usize,
        wanted: &[usize],
    ) -> Result<Vec<Vec<u64>>, TrajectoryError> {
        let (n_frames, payload) = self.read_payload(chunk)?;
        let view = PayloadView::parse(&payload, ncols(self.header.n_atoms))?;
        let max_key = self.header.precision.max_key();
        wanted.iter().map(|&c| decode_column(view.column(c), n_frames as usize, max_key)).collect()
    }

    /// Index of the chunk holding frame `index`, which must be in range.
    fn chunk_of(&self, index: u64) -> usize {
        self.chunks.partition_point(|c| c.first_frame <= index) - 1
    }

    /// Reads frame `index` by seeking straight to its chunk.
    ///
    /// # Errors
    /// Returns [`TrajectoryError::FrameOutOfRange`] past the end, or an error if
    /// the chunk is corrupt.
    pub fn frame(&mut self, index: u64) -> Result<Frame, TrajectoryError> {
        if index >= self.n_frames {
            return Err(TrajectoryError::FrameOutOfRange { index, n_frames: self.n_frames });
        }
        let chunk = self.chunk_of(index);
        let row = (index - self.chunks[chunk].first_frame) as usize;
        let all: Vec<usize> = (0..ncols(self.header.n_atoms)).collect();
        let cols = self.read_columns(chunk, &all)?;
        Ok(self.frame_from_columns(&cols, row))
    }

    /// Streams all frames in order, holding one decoded chunk at a time.
    pub fn frames(&mut self) -> FrameIter<'_, R> {
        FrameIter { reader: self, next_chunk: 0, current: Vec::new().into_iter(), failed: false }
    }

    fn frame_from_columns(&self, cols: &[Vec<u64>], row: usize) -> Frame {
        let p = self.header.precision;
        let n = self.header.n_atoms;
        let positions = (0..n)
            .map(|a| {
                let base = N_LEADING_COLUMNS + 3 * a;
                [
                    p.value(cols[base][row]),
                    p.value(cols[base + 1][row]),
                    p.value(cols[base + 2][row]),
                ]
            })
            .collect();
        Frame {
            step: cols[0][row],
            time: p.value(cols[1][row]),
            box_lengths: [p.value(cols[2][row]), p.value(cols[3][row]), p.value(cols[4][row])],
            positions,
        }
    }

    /// Decodes a whole chunk into frames.
    fn load_chunk(&mut self, chunk: usize) -> Result<Vec<Frame>, TrajectoryError> {
        let n_frames = self.chunks[chunk].n_frames;
        let all: Vec<usize> = (0..ncols(self.header.n_atoms)).collect();
        let cols = self.read_columns(chunk, &all)?;
        Ok((0..n_frames as usize).map(|row| self.frame_from_columns(&cols, row)).collect())
    }
}

/// Streaming iterator over the frames of a [`TrajectoryReader`].
///
/// At most one decoded chunk is resident at a time. See
/// [`FrameIter::resident_frames`].
#[derive(Debug)]
pub struct FrameIter<'a, R: Read + Seek> {
    reader: &'a mut TrajectoryReader<R>,
    next_chunk: usize,
    current: std::vec::IntoIter<Frame>,
    failed: bool,
}

impl<R: Read + Seek> FrameIter<'_, R> {
    /// Frames of the current chunk that are decoded and not yet yielded.
    pub fn resident_frames(&self) -> usize {
        self.current.len()
    }
}

impl<R: Read + Seek> Iterator for FrameIter<'_, R> {
    type Item = Result<Frame, TrajectoryError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        loop {
            if let Some(f) = self.current.next() {
                return Some(Ok(f));
            }
            if self.next_chunk >= self.reader.chunks.len() {
                return None;
            }
            match self.reader.load_chunk(self.next_chunk) {
                Ok(frames) => {
                    self.next_chunk += 1;
                    self.current = frames.into_iter();
                }
                Err(e) => {
                    self.failed = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Lazy-engine source
// ---------------------------------------------------------------------------

/// Trajectory exposed to the lazy engine as one row per frame.
///
/// Columns: `step` (`I64`), `time`, `box_x`, `box_y`, `box_z`, then `x{a}`,
/// `y{a}`, `z{a}` for each atom `a` (all `F64`). Projection pushdown means only
/// the requested columns are decoded, and only for the chunks a row range
/// touches.
#[derive(Debug)]
pub struct LazyTrajectory<R: Read + Seek + Send> {
    reader: Mutex<TrajectoryReader<R>>,
    schema: Schema,
}

impl<R: Read + Seek + Send> LazyTrajectory<R> {
    /// Wraps an opened reader.
    pub fn new(reader: TrajectoryReader<R>) -> Self {
        let n_atoms = reader.n_atoms();
        let mut schema = vec![
            Field { name: "step".into(), dtype: LazyDType::I64 },
            Field { name: "time".into(), dtype: LazyDType::F64 },
            Field { name: "box_x".into(), dtype: LazyDType::F64 },
            Field { name: "box_y".into(), dtype: LazyDType::F64 },
            Field { name: "box_z".into(), dtype: LazyDType::F64 },
        ];
        for a in 0..n_atoms {
            for c in ["x", "y", "z"] {
                schema.push(Field { name: format!("{c}{a}"), dtype: LazyDType::F64 });
            }
        }
        Self { reader: Mutex::new(reader), schema }
    }

    /// Recovers the reader.
    ///
    /// # Errors
    /// Returns an error if a previous read panicked while holding the lock.
    pub fn into_reader(self) -> Result<TrajectoryReader<R>, TrajectoryError> {
        self.reader
            .into_inner()
            .map_err(|_| TrajectoryError::Corrupt("trajectory reader lock poisoned".into()))
    }

    fn column_index(&self, name: &str) -> Option<usize> {
        self.schema.iter().position(|f| f.name == name)
    }
}

impl<R: Read + Seek + Send + 'static> Source for LazyTrajectory<R> {
    fn schema(&self) -> Schema {
        self.schema.clone()
    }

    fn num_rows(&self) -> usize {
        match self.reader.lock() {
            Ok(r) => r.n_frames() as usize,
            Err(_) => 0,
        }
    }

    fn read(
        &self,
        rows: std::ops::Range<usize>,
        projection: &[String],
    ) -> Result<DataFrame, LazyError> {
        let mut reader = self
            .reader
            .lock()
            .map_err(|_| LazyError::Invalid("trajectory reader lock poisoned".into()))?;
        let total = reader.n_frames() as usize;
        let end = rows.end.min(total);
        let start = rows.start.min(end);

        let mut wanted = Vec::with_capacity(projection.len());
        for name in projection {
            let idx =
                self.column_index(name).ok_or_else(|| LazyError::UnknownColumn(name.clone()))?;
            wanted.push(idx);
        }
        let mut out: Vec<Vec<u64>> = vec![Vec::with_capacity(end - start); wanted.len()];
        if start < end {
            let chunks = reader.chunks().to_vec();
            let mut c = chunks.partition_point(|ch| ch.first_frame <= start as u64) - 1;
            while c < chunks.len() && (chunks[c].first_frame as usize) < end {
                let first = chunks[c].first_frame as usize;
                let n = chunks[c].n_frames as usize;
                let lo = start.max(first) - first;
                let hi = end.min(first + n) - first;
                let cols = reader
                    .read_columns(c, &wanted)
                    .map_err(|e| LazyError::Invalid(e.to_string()))?;
                for (dst, src) in out.iter_mut().zip(&cols) {
                    dst.extend_from_slice(&src[lo..hi]);
                }
                c += 1;
            }
        }
        let p = reader.header().precision;
        let columns = projection.iter().zip(out).zip(&wanted).map(|((name, keys), &idx)| {
            let col = if idx == 0 {
                Column::I64(keys.iter().map(|&k| k.cast_signed()).collect())
            } else {
                Column::F64(keys.iter().map(|&k| p.value(k)).collect())
            };
            (name.clone(), col)
        });
        DataFrame::with_rows(end - start, columns)
    }
}

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

fn read_u32(b: &[u8]) -> u32 {
    let mut a = [0u8; 4];
    a.copy_from_slice(&b[..4]);
    u32::from_le_bytes(a)
}

fn read_u64(b: &[u8]) -> u64 {
    let mut a = [0u8; 8];
    a.copy_from_slice(&b[..8]);
    u64::from_le_bytes(a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SplitMix64;
    use std::io::Cursor;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::Arc;

    /// Counts bytes returned by `read` and bytes skipped by `seek` calls.
    struct Counting<R> {
        inner: R,
        bytes: Arc<AtomicU64>,
    }

    impl<R: Read> Read for Counting<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.bytes.fetch_add(n as u64, Ordering::SeqCst);
            Ok(n)
        }
    }

    impl<R: Seek> Seek for Counting<R> {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            self.inner.seek(pos)
        }
    }

    fn smooth_frames(n_frames: usize, n_atoms: usize) -> Vec<Frame> {
        (0..n_frames)
            .map(|s| {
                let t = s as f64 * 0.002;
                Frame {
                    step: s as u64 * 10,
                    time: t,
                    box_lengths: [3.0, 3.1, 3.2],
                    positions: (0..n_atoms)
                        .map(|a| {
                            let a = a as f64;
                            [
                                0.3 * a + 0.01 * (t + a).sin(),
                                0.2 * a + 0.01 * (t * 1.3 + a).cos(),
                                0.1 * a + 0.005 * t,
                            ]
                        })
                        .collect(),
                }
            })
            .collect()
    }

    fn write_all(frames: &[Frame], n_atoms: usize, p: Precision, k: u32) -> Vec<u8> {
        let mut w = TrajectoryWriter::new(Vec::new(), n_atoms, p, k).unwrap();
        for f in frames {
            w.push_frame(f).unwrap();
        }
        w.finish().unwrap()
    }

    fn f64s(c: &Column) -> &[f64] {
        match c {
            Column::F64(v) => v,
            _ => panic!("expected an f64 column"),
        }
    }

    fn i64s(c: &Column) -> &[i64] {
        match c {
            Column::I64(v) => v,
            _ => panic!("expected an i64 column"),
        }
    }

    fn bits(f: &Frame) -> Vec<u64> {
        let mut v = vec![f.step, f.time.to_bits()];
        v.extend(f.box_lengths.iter().map(|x| x.to_bits()));
        for p in &f.positions {
            v.extend(p.iter().map(|x| x.to_bits()));
        }
        v
    }

    #[test]
    fn zigzag_and_varint_round_trip_edges() {
        for d in [0i64, 1, -1, 63, -64, 64, i64::MAX, i64::MIN] {
            assert_eq!(unzigzag(zigzag(d)), d);
        }
        for v in [0u64, 1, 127, 128, 300, u64::MAX / 2, u64::MAX] {
            let mut buf = Vec::new();
            put_varint(&mut buf, v);
            let mut pos = 0;
            assert_eq!(get_varint(&buf, &mut pos).unwrap(), v);
            assert_eq!(pos, buf.len());
        }
        // Eleven continuation bytes overflow u64.
        let bad = [0xffu8; 11];
        assert!(get_varint(&bad, &mut 0).is_err());
    }

    #[test]
    fn float_keys_are_order_preserving_bijections() {
        let samples = [f64::NEG_INFINITY, -1.5, -0.0, 0.0, 1e-300, 2.0, f64::INFINITY];
        let keys: Vec<u64> = samples.iter().map(|&v| Precision::F64.key(v)).collect();
        assert!(keys.windows(2).all(|w| w[0] < w[1]));
        for &v in &samples {
            assert_eq!(Precision::F64.value(Precision::F64.key(v)).to_bits(), v.to_bits());
            let r = Precision::F32.value(Precision::F32.key(v));
            assert_eq!(r.to_bits(), f64::from(v as f32).to_bits());
        }
        let nan = f64::from_bits(0x7ff8_0000_0000_1234);
        assert_eq!(Precision::F64.value(Precision::F64.key(nan)).to_bits(), nan.to_bits());
    }

    #[test]
    fn round_trip_is_bit_exact_for_f64() {
        let mut frames = smooth_frames(7, 5);
        frames[2].positions[1][0] = -0.0;
        frames[3].positions[0][2] = f64::NAN;
        frames[4].time = f64::MIN_POSITIVE / 4.0;
        frames[5].positions[4][1] = -1.0e300;
        let bytes = write_all(&frames, 5, Precision::F64, 3);
        let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        assert_eq!(r.n_frames(), 7);
        assert_eq!(r.chunks().len(), 3);
        let got: Vec<Frame> = r.frames().map(|f| f.unwrap()).collect();
        assert_eq!(got.len(), frames.len());
        for (a, b) in frames.iter().zip(&got) {
            assert_eq!(bits(a), bits(b));
        }
    }

    #[test]
    fn f32_precision_rounds_once_then_is_exact() {
        let frames = smooth_frames(4, 3);
        let bytes = write_all(&frames, 3, Precision::F32, 2);
        let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        assert_eq!(r.header().precision, Precision::F32);
        for (i, f) in frames.iter().enumerate() {
            let got = r.frame(i as u64).unwrap();
            assert_eq!(got.positions[2][0], f64::from(f.positions[2][0] as f32));
            assert_eq!(got.time, f64::from(f.time as f32));
        }
    }

    #[test]
    fn random_access_matches_streaming() {
        let frames = smooth_frames(40, 4);
        let bytes = write_all(&frames, 4, Precision::F64, 7);
        let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        let streamed: Vec<Frame> = r.frames().map(|f| f.unwrap()).collect();
        for i in [0u64, 6, 7, 13, 27, 39] {
            assert_eq!(r.frame(i).unwrap(), streamed[i as usize]);
        }
        assert!(matches!(
            r.frame(40),
            Err(TrajectoryError::FrameOutOfRange { index: 40, n_frames: 40 })
        ));
    }

    #[test]
    fn empty_trajectory_is_valid() {
        let bytes =
            TrajectoryWriter::new(Vec::new(), 3, Precision::F64, 4).unwrap().finish().unwrap();
        let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        assert_eq!(r.n_frames(), 0);
        assert!(r.frames().next().is_none());
    }

    #[test]
    fn writer_validates_input() {
        assert!(TrajectoryWriter::new(Vec::new(), 1, Precision::F64, 0).is_err());
        let mut w = TrajectoryWriter::new(Vec::new(), 2, Precision::F64, 4).unwrap();
        let bad = Frame { positions: vec![[0.0; 3]], ..Frame::default() };
        assert!(matches!(
            w.push_frame(&bad),
            Err(TrajectoryError::AtomCountMismatch { expected: 2, found: 1 })
        ));
        w.close().unwrap();
        w.close().unwrap();
        assert!(matches!(w.push_frame(&Frame::default()), Err(TrajectoryError::Closed)));
    }

    #[test]
    fn compression_beats_raw_f32_on_smooth_data() {
        let frames = smooth_frames(400, 30);
        let bytes = write_all(&frames, 30, Precision::F32, 64);
        let raw_f32 = 400 * (8 + 8 + 24 + 30 * 12);
        assert!(
            bytes.len() < raw_f32,
            "compressed {} bytes, raw f32 {} bytes",
            bytes.len(),
            raw_f32
        );
    }

    #[test]
    fn every_single_byte_corruption_is_detected() {
        let frames = smooth_frames(5, 2);
        let bytes = write_all(&frames, 2, Precision::F64, 2);
        for i in 0..bytes.len() {
            let mut bad = bytes.clone();
            bad[i] ^= 0x10;
            let ok = match TrajectoryReader::open(Cursor::new(bad)) {
                Err(_) => false,
                Ok(mut r) => r.frames().all(|f| f.is_ok()),
            };
            assert!(!ok, "corruption at byte {i} went undetected");
        }
    }

    #[test]
    fn every_truncation_is_rejected() {
        let bytes = write_all(&smooth_frames(5, 2), 2, Precision::F64, 2);
        for n in 0..bytes.len() {
            let ok = match TrajectoryReader::open(Cursor::new(bytes[..n].to_vec())) {
                Err(_) => false,
                Ok(mut r) => r.frames().all(|f| f.is_ok()),
            };
            assert!(!ok, "prefix of {n} bytes accepted");
        }
    }

    #[test]
    fn corrupt_chunk_only_fails_its_own_frames() {
        let frames = smooth_frames(10, 2);
        let bytes = write_all(&frames, 2, Precision::F64, 4);
        let r = TrajectoryReader::open(Cursor::new(bytes.clone())).unwrap();
        let target = r.chunks()[1];
        let mut bad = bytes;
        let pos = (target.record_offset + RECORD_HEADER_LEN + 3) as usize;
        bad[pos] ^= 0xff;
        let mut r = TrajectoryReader::open(Cursor::new(bad)).unwrap();
        assert!(r.frame(0).is_ok());
        assert!(matches!(r.frame(5), Err(TrajectoryError::ChecksumMismatch { .. })));
        assert!(r.frame(9).is_ok());
    }

    #[test]
    fn long_trajectory_streams_with_one_chunk_resident() {
        let n_frames = 5000usize;
        let k = 256u32;
        let frames = smooth_frames(n_frames, 4);
        let bytes = write_all(&frames, 4, Precision::F32, k);
        let mut r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        let n_chunks = r.chunks().len() as u64;
        assert_eq!(n_chunks, (n_frames as u64).div_ceil(u64::from(k)));
        let mut seen = 0u64;
        let mut max_resident = 0usize;
        let mut it = r.frames();
        while let Some(f) = it.next() {
            let f = f.unwrap();
            assert_eq!(f.step, seen * 10);
            seen += 1;
            max_resident = max_resident.max(it.resident_frames());
        }
        assert_eq!(seen, n_frames as u64);
        assert!(max_resident < k as usize, "resident {max_resident} frames");
        assert_eq!(r.stats().chunks_loaded, n_chunks);
    }

    #[test]
    fn random_access_reads_only_the_needed_chunk() {
        let n_frames = 3000usize;
        let frames = smooth_frames(n_frames, 4);
        let bytes = write_all(&frames, 4, Precision::F64, 128);
        let counter = Arc::new(AtomicU64::new(0));
        let src = Counting { inner: Cursor::new(bytes), bytes: Arc::clone(&counter) };
        let mut r = TrajectoryReader::open(src).unwrap();
        let info = r.chunks()[r.chunk_of(2500)];
        counter.store(0, Ordering::SeqCst);
        r.reset_stats();
        let f = r.frame(2500).unwrap();
        assert_eq!(f.step, 25_000);
        assert_eq!(r.stats(), ReadStats { chunks_loaded: 1, bytes_read: info.record_len });
        assert_eq!(counter.load(Ordering::SeqCst), info.record_len);
    }

    #[test]
    fn lazy_source_projects_columns_and_rows() {
        use crate::lazy::{col, LazyFrame};
        let frames = smooth_frames(20, 2);
        let bytes = write_all(&frames, 2, Precision::F64, 6);
        let r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        let src = LazyTrajectory::new(r);
        assert_eq!(src.num_rows(), 20);
        let names = vec!["step".to_string(), "x1".to_string()];
        let df = src.read(5..9, &names).unwrap();
        assert_eq!(df.len(), 4);
        assert_eq!(i64s(df.column("step").unwrap()), &[50, 60, 70, 80]);
        let x1 = f64s(df.column("x1").unwrap());
        assert_eq!(x1[0], frames[5].positions[1][0]);
        assert!(matches!(src.read(0..1, &["nope".into()]), Err(LazyError::UnknownColumn(_))));

        let out = LazyFrame::scan(src).filter(col("time").gt(0.02)).collect().unwrap();
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn lazy_source_matches_eager_random_data() {
        let mut rng = SplitMix64::new(7);
        let frames: Vec<Frame> = (0..30)
            .map(|s| Frame {
                step: s,
                time: rng.uniform(0.0, 1.0),
                box_lengths: [2.0, 2.0, 2.0],
                positions: (0..3)
                    .map(|_| [rng.uniform(0.0, 1.0), rng.uniform(0.0, 1.0), rng.uniform(0.0, 1.0)])
                    .collect(),
            })
            .collect();
        let bytes = write_all(&frames, 3, Precision::F64, 8);
        let r = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
        let src = LazyTrajectory::new(r);
        let all: Vec<String> = src.schema().into_iter().map(|f| f.name).collect();
        let df = src.read(0..30, &all).unwrap();
        for (i, f) in frames.iter().enumerate() {
            assert_eq!(f64s(df.column("time").unwrap())[i], f.time);
            assert_eq!(f64s(df.column("z2").unwrap())[i], f.positions[2][2]);
        }
    }
}
