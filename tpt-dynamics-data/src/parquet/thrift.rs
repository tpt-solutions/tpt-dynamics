//! Minimal Thrift compact-protocol reader, enough for Parquet metadata.

use super::ParquetError;

/// Compact-protocol type codes.
pub(crate) const T_BOOL_TRUE: u8 = 1;
pub(crate) const T_BOOL_FALSE: u8 = 2;
pub(crate) const T_BYTE: u8 = 3;
pub(crate) const T_I16: u8 = 4;
pub(crate) const T_I32: u8 = 5;
pub(crate) const T_I64: u8 = 6;
pub(crate) const T_DOUBLE: u8 = 7;
pub(crate) const T_BINARY: u8 = 8;
pub(crate) const T_LIST: u8 = 9;
pub(crate) const T_SET: u8 = 10;
pub(crate) const T_MAP: u8 = 11;
pub(crate) const T_STRUCT: u8 = 12;

/// Cursor over a Thrift compact-encoded byte slice.
pub(crate) struct Thrift<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Thrift<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }

    /// Current byte offset into the slice.
    pub(crate) fn position(&self) -> usize {
        self.pos
    }

    fn byte(&mut self) -> Result<u8, ParquetError> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or_else(|| ParquetError::Format("truncated metadata".into()))?;
        self.pos += 1;
        Ok(b)
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], ParquetError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.buf.len())
            .ok_or_else(|| ParquetError::Format("truncated metadata".into()))?;
        let out = &self.buf[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    /// Unsigned LEB128 varint.
    pub(crate) fn varint(&mut self) -> Result<u64, ParquetError> {
        let mut out = 0u64;
        for shift in (0..64).step_by(7) {
            let b = self.byte()?;
            out |= u64::from(b & 0x7f) << shift;
            if b & 0x80 == 0 {
                return Ok(out);
            }
        }
        Err(ParquetError::Format("varint too long".into()))
    }

    /// Zigzag-encoded signed varint.
    // The casts reinterpret the bit pattern on purpose; that is what zigzag means.
    #[allow(clippy::cast_possible_wrap)]
    pub(crate) fn zigzag(&mut self) -> Result<i64, ParquetError> {
        let u = self.varint()?;
        Ok(((u >> 1) as i64) ^ -((u & 1) as i64))
    }

    /// Length-prefixed binary value.
    pub(crate) fn binary(&mut self) -> Result<&'a [u8], ParquetError> {
        let n = self.varint()? as usize;
        self.take(n)
    }

    /// Length-prefixed UTF-8 string.
    pub(crate) fn string(&mut self) -> Result<String, ParquetError> {
        let b = self.binary()?;
        String::from_utf8(b.to_vec())
            .map_err(|_| ParquetError::Format("metadata string is not UTF-8".into()))
    }

    /// Reads a field header. Returns `None` at the struct's stop marker.
    pub(crate) fn field(&mut self, last_id: &mut i16) -> Result<Option<(i16, u8)>, ParquetError> {
        let b = self.byte()?;
        if b == 0 {
            return Ok(None);
        }
        let delta = i16::from(b >> 4);
        let ttype = b & 0x0f;
        let id = if delta == 0 { self.zigzag()? as i16 } else { *last_id + delta };
        *last_id = id;
        Ok(Some((id, ttype)))
    }

    /// Reads a list header: element count and element type.
    pub(crate) fn list(&mut self) -> Result<(usize, u8), ParquetError> {
        let b = self.byte()?;
        let size = match b >> 4 {
            15 => self.varint()? as usize,
            n => usize::from(n),
        };
        Ok((size, b & 0x0f))
    }

    /// Reads a struct, calling `f` for each field with (reader, id, type).
    /// `f` must consume the field's value, or return `Ok(false)` to have the
    /// value skipped.
    pub(crate) fn structure<F>(&mut self, mut f: F) -> Result<(), ParquetError>
    where
        F: FnMut(&mut Self, i16, u8) -> Result<bool, ParquetError>,
    {
        let mut last = 0i16;
        while let Some((id, ttype)) = self.field(&mut last)? {
            if !f(self, id, ttype)? {
                self.skip(ttype)?;
            }
        }
        Ok(())
    }

    /// Skips one value of the given compact type.
    pub(crate) fn skip(&mut self, ttype: u8) -> Result<(), ParquetError> {
        match ttype {
            T_BOOL_TRUE | T_BOOL_FALSE => Ok(()),
            T_BYTE => self.take(1).map(|_| ()),
            T_I16 | T_I32 | T_I64 => self.varint().map(|_| ()),
            T_DOUBLE => self.take(8).map(|_| ()),
            T_BINARY => self.binary().map(|_| ()),
            T_LIST | T_SET => {
                let (n, et) = self.list()?;
                for _ in 0..n {
                    self.skip(et)?;
                }
                Ok(())
            }
            T_MAP => {
                let n = self.varint()? as usize;
                if n > 0 {
                    let kv = self.byte()?;
                    for _ in 0..n {
                        self.skip(kv >> 4)?;
                        self.skip(kv & 0x0f)?;
                    }
                }
                Ok(())
            }
            T_STRUCT => self.structure(|_, _, _| Ok(false)),
            other => Err(ParquetError::Format(format!("unknown Thrift type {other}"))),
        }
    }
}
