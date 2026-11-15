//! In-process model of distributed histogram filling.
//!
//! Each "node" is a scoped thread that fills its own partial histogram, encodes it
//! with the binary serialiser, and sends only bytes back. The coordinator decodes
//! the partials and merges them in a different order from the one they were
//! produced in. The result must match a single serial fill. No network is used;
//! distributed execution of lazy plans is not implemented.

use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::rng::SplitMix64;

const NEVENTS: usize = 200_000;
const NODES: usize = 4;

fn events(seed: u64) -> Vec<(f64, f64)> {
    let mut rng = SplitMix64::new(seed);
    (0..NEVENTS)
        .map(|_| {
            // Gaussian around 0; the +-6 range leaves a few values in the flows.
            let x = rng.gaussian() * 2.0;
            // Unit weights keep bin sums integral, so the comparison is exact.
            (x, 1.0)
        })
        .collect()
}

fn fresh() -> Histogram1D {
    Histogram1D::uniform(60, -6.0, 6.0).unwrap()
}

fn serial_fill(data: &[(f64, f64)]) -> Histogram1D {
    let mut h = fresh();
    for &(x, w) in data {
        h.fill_weighted(x, w);
    }
    h
}

/// Fill one chunk per scoped thread and return each partial as serialised bytes.
fn worker_bytes(data: &[(f64, f64)]) -> Vec<Vec<u8>> {
    let chunk = data.len().div_ceil(NODES);
    std::thread::scope(|s| {
        let handles: Vec<_> = data
            .chunks(chunk)
            .map(|part| {
                s.spawn(move || {
                    let mut h = fresh();
                    for &(x, w) in part {
                        h.fill_weighted(x, w);
                    }
                    h.to_bytes()
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    })
}

fn assert_same_counts(a: &Histogram1D, b: &Histogram1D) {
    assert_eq!(a.edges(), b.edges());
    assert_eq!(a.contents(), b.contents(), "bin contents differ");
    assert_eq!(a.sum_weights_squared(), b.sum_weights_squared(), "sumw2 differs");
    assert_eq!(a.underflow(), b.underflow());
    assert_eq!(a.overflow(), b.overflow());
    assert_eq!(a.entries(), b.entries());
    assert_eq!(a.nan_fills(), b.nan_fills());
    let (ma, mb) = (a.mean().unwrap(), b.mean().unwrap());
    assert!((ma - mb).abs() < 1e-9, "mean differs: {ma} vs {mb}");
}

#[test]
fn merged_partials_equal_serial_fill_in_any_order() {
    let data = events(42);
    let serial = serial_fill(&data);

    let bytes = worker_bytes(&data);
    assert_eq!(bytes.len(), NODES);

    let partials: Vec<Histogram1D> =
        bytes.iter().map(|b| Histogram1D::from_bytes(b).unwrap()).collect();

    // Merge in several permutations; every order must agree with the serial fill.
    for order in [[0usize, 1, 2, 3], [3, 1, 0, 2], [2, 3, 1, 0]] {
        let mut merged = fresh();
        for &i in &order {
            merged.merge(&partials[i]).unwrap();
        }
        assert_same_counts(&merged, &serial);
    }
}

#[test]
fn tree_reduction_matches_linear_merge() {
    let data = events(7);
    let partials: Vec<Histogram1D> =
        worker_bytes(&data).iter().map(|b| Histogram1D::from_bytes(b).unwrap()).collect();

    // Pairwise (tree) reduction, as a real coordinator would do with many nodes.
    let mut level = partials;
    while level.len() > 1 {
        let mut next = Vec::new();
        let mut it = level.into_iter();
        while let Some(mut a) = it.next() {
            if let Some(b) = it.next() {
                a.merge(&b).unwrap();
            }
            next.push(a);
        }
        level = next;
    }
    assert_same_counts(&level[0], &serial_fill(&data));
}

#[test]
fn serialised_partial_survives_round_trip_bit_exact() {
    let data = events(99);
    let part = serial_fill(&data[..1000]);
    let bytes = part.to_bytes();
    let back = Histogram1D::from_bytes(&bytes).unwrap();
    assert_eq!(back, part);
    assert_eq!(back.to_bytes(), bytes);
}

#[test]
fn corrupted_partial_is_rejected_by_coordinator() {
    let data = events(3);
    let mut bytes = worker_bytes(&data).remove(0);
    let mid = bytes.len() / 2;
    bytes[mid] ^= 0x01;
    assert!(Histogram1D::from_bytes(&bytes).is_err());
}

#[test]
fn mismatched_binning_is_rejected_on_merge() {
    let mut a = fresh();
    let b = Histogram1D::uniform(59, -6.0, 6.0).unwrap();
    assert!(a.merge(&b).is_err());
}
