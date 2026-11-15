//! Memory bound for streaming trajectory reads (todo Phase 6).
//!
//! The reader decodes one chunk at a time, so an iteration over a trajectory
//! of any length holds at most one chunk of decoded frames. This test checks
//! that bound on a trajectory much larger than one chunk. The multi-gigabyte
//! measurement lives in `docs/benchmarks/README.md`.

use std::io::Cursor;

use tpt_dynamics_data::rng::SplitMix64;
use tpt_dynamics_data::trajectory::{
    Frame, FrameSink, Precision, TrajectoryReader, TrajectoryWriter,
};

const ATOMS: usize = 500;
const FRAMES: u64 = 301; // not a multiple of the chunk size, so the last chunk is short
const CHUNK: u32 = 16;

fn make_frames() -> Vec<Frame> {
    let mut rng = SplitMix64::new(77);
    let mut positions: Vec<[f64; 3]> = (0..ATOMS)
        .map(|_| [rng.uniform(0.0, 5.0), rng.uniform(0.0, 5.0), rng.uniform(0.0, 5.0)])
        .collect();
    (0..FRAMES)
        .map(|step| {
            for p in &mut positions {
                for c in p.iter_mut() {
                    *c += rng.uniform(-0.01, 0.01);
                }
            }
            Frame {
                step,
                time: step as f64 * 0.002,
                box_lengths: [5.0; 3],
                positions: positions.clone(),
            }
        })
        .collect()
}

#[test]
fn streaming_iteration_holds_at_most_one_chunk() {
    let frames = make_frames();
    let mut w = TrajectoryWriter::new(Vec::new(), ATOMS, Precision::F64, CHUNK).unwrap();
    for f in &frames {
        w.push_frame(f).unwrap();
    }
    let bytes = w.finish().unwrap();

    let mut reader = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
    assert_eq!(reader.n_frames(), FRAMES);
    let mut it = reader.frames();
    let mut seen = 0u64;
    while let Some(frame) = it.next() {
        let frame = frame.unwrap();
        assert!(
            it.resident_frames() < CHUNK as usize,
            "frame {seen}: {} decoded frames resident, bound is one chunk",
            it.resident_frames()
        );
        assert_eq!(frame, frames[seen as usize], "frame {seen} differs from the written frame");
        seen += 1;
    }
    assert_eq!(seen, FRAMES);
}

#[test]
fn random_access_matches_streaming() {
    let frames = make_frames();
    let mut w = TrajectoryWriter::new(Vec::new(), ATOMS, Precision::F64, CHUNK).unwrap();
    for f in &frames {
        w.push_frame(f).unwrap();
    }
    let bytes = w.finish().unwrap();
    let mut reader = TrajectoryReader::open(Cursor::new(bytes)).unwrap();
    for index in [0u64, 15, 16, 17, 200, FRAMES - 1] {
        assert_eq!(reader.frame(index).unwrap(), frames[index as usize], "frame {index}");
    }
}
