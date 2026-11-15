//! Property tests for histogram `merge`.
//!
//! For randomly generated event sets with non-integer weights, partitioned at
//! random, the merged result must agree with a single serial fill for every
//! histogram type. `merge` must also be commutative and associative. Weights are
//! not integers here, so floating-point sums are compared within a relative
//! tolerance rather than bit-for-bit.

use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::hist2d::Histogram2D;
use tpt_dynamics_data::profile::Profile1D;
use tpt_dynamics_data::rng::SplitMix64;

const CASES: u64 = 64;
const REL_TOL: f64 = 1e-12;

#[derive(Clone, Copy, Debug)]
struct Event {
    x: f64,
    y: f64,
    w: f64,
}

fn events(rng: &mut SplitMix64, n: usize) -> Vec<Event> {
    (0..n)
        .map(|_| Event {
            // Range [-6, 6] with some values outside, so the flows are exercised.
            x: rng.gaussian() * 3.0,
            y: rng.gaussian() * 3.0,
            w: 0.1 + 2.0 * rng.next_f64(),
        })
        .collect()
}

/// Random cut points: returns `parts + 1` monotone indices into `0..=n`.
fn random_cuts(rng: &mut SplitMix64, n: usize, parts: usize) -> Vec<usize> {
    let mut cuts: Vec<usize> =
        (0..parts - 1).map(|_| (rng.next_u64() % (n as u64 + 1)) as usize).collect();
    cuts.push(0);
    cuts.push(n);
    cuts.sort_unstable();
    cuts
}

fn approx_eq(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= REL_TOL * scale
}

fn assert_close(a: f64, b: f64, what: &str) {
    assert!(approx_eq(a, b), "{what}: {a} vs {b}");
}

fn fill1(data: &[Event]) -> Histogram1D {
    let mut h = Histogram1D::uniform(40, -6.0, 6.0).unwrap();
    for e in data {
        h.fill_weighted(e.x, e.w);
    }
    h
}

fn fill2(data: &[Event]) -> Histogram2D {
    let mut h = Histogram2D::uniform(12, -6.0, 6.0, 12, -6.0, 6.0).unwrap();
    for e in data {
        h.fill_weighted(e.x, e.y, e.w);
    }
    h
}

fn fillp(data: &[Event]) -> Profile1D {
    let mut p = Profile1D::uniform(20, -6.0, 6.0).unwrap();
    for e in data {
        p.fill_weighted(e.x, e.y, e.w);
    }
    p
}

fn assert_h1_equal(a: &Histogram1D, b: &Histogram1D) {
    assert_eq!(a.edges(), b.edges());
    assert_eq!(a.entries(), b.entries(), "entries");
    assert_eq!(a.nan_fills(), b.nan_fills(), "nan_fills");
    for (i, (&ca, &cb)) in a.contents().iter().zip(b.contents()).enumerate() {
        assert_close(ca, cb, &format!("bin {i} content"));
    }
    for (i, (&ea, &eb)) in a.sum_weights_squared().iter().zip(b.sum_weights_squared()).enumerate() {
        assert_close(ea, eb, &format!("bin {i} sumw2"));
    }
    assert_close(a.underflow(), b.underflow(), "underflow");
    assert_close(a.overflow(), b.overflow(), "overflow");
    assert_close(a.integral(), b.integral(), "integral");
}

fn assert_h2_equal(a: &Histogram2D, b: &Histogram2D) {
    assert_eq!(a.x_edges(), b.x_edges());
    assert_eq!(a.y_edges(), b.y_edges());
    assert_eq!(a.entries(), b.entries(), "entries");
    for ix in 0..a.nbins_x() {
        for iy in 0..a.nbins_y() {
            assert_close(a.bin(ix, iy), b.bin(ix, iy), &format!("cell ({ix},{iy})"));
        }
    }
    assert_close(a.integral(), b.integral(), "integral");
}

fn assert_profile_equal(a: &Profile1D, b: &Profile1D) {
    assert_eq!(a.edges(), b.edges());
    for i in 0..a.nbins() {
        let (ca, cb) = (a.cell(i), b.cell(i));
        assert_eq!(ca.entries, cb.entries, "bin {i} entries");
        assert_close(ca.sumw, cb.sumw, &format!("bin {i} sumw"));
        assert_close(ca.sumwy, cb.sumwy, &format!("bin {i} sumwy"));
        assert_close(ca.sumwy2, cb.sumwy2, &format!("bin {i} sumwy2"));
        assert_close(ca.sumw2, cb.sumw2, &format!("bin {i} sumw2"));
    }
}

#[test]
fn merged_partitions_match_serial_fill_for_all_histogram_types() {
    for case in 0..CASES {
        let mut rng = SplitMix64::new(0xA5A5_0000 + case);
        let n = 1 + (rng.next_u64() % 4000) as usize;
        let parts = 2 + (rng.next_u64() % 6) as usize;
        let data = events(&mut rng, n);
        let cuts = random_cuts(&mut rng, n, parts);

        let serial1 = fill1(&data);
        let serial2 = fill2(&data);
        let serialp = fillp(&data);

        let mut m1 = Histogram1D::uniform(40, -6.0, 6.0).unwrap();
        let mut m2 = Histogram2D::uniform(12, -6.0, 6.0, 12, -6.0, 6.0).unwrap();
        let mut mp = Profile1D::uniform(20, -6.0, 6.0).unwrap();
        for w in cuts.windows(2) {
            let chunk = &data[w[0]..w[1]];
            m1.merge(&fill1(chunk)).unwrap();
            m2.merge(&fill2(chunk)).unwrap();
            mp.merge(&fillp(chunk)).unwrap();
        }

        assert_h1_equal(&m1, &serial1);
        assert_h2_equal(&m2, &serial2);
        assert_profile_equal(&mp, &serialp);
    }
}

#[test]
fn merge_is_commutative_and_associative() {
    for case in 0..CASES {
        let mut rng = SplitMix64::new(0x5EED_0000 + case);
        let data = events(&mut rng, 300);
        let cuts = random_cuts(&mut rng, data.len(), 3);
        let (a, b, c) = (&data[cuts[0]..cuts[1]], &data[cuts[1]..cuts[2]], &data[cuts[2]..cuts[3]]);

        // (a + b) + c
        let mut left = fill1(a);
        left.merge(&fill1(b)).unwrap();
        left.merge(&fill1(c)).unwrap();

        // c + (b + a)
        let mut ba = fill1(b);
        ba.merge(&fill1(a)).unwrap();
        let mut right = fill1(c);
        right.merge(&ba).unwrap();

        assert_h1_equal(&left, &right);

        let mut p_left = fillp(a);
        p_left.merge(&fillp(b)).unwrap();
        p_left.merge(&fillp(c)).unwrap();
        let mut p_ba = fillp(b);
        p_ba.merge(&fillp(a)).unwrap();
        let mut p_right = fillp(c);
        p_right.merge(&p_ba).unwrap();
        assert_profile_equal(&p_left, &p_right);

        let mut h2_left = fill2(a);
        h2_left.merge(&fill2(b)).unwrap();
        h2_left.merge(&fill2(c)).unwrap();
        let mut h2_ba = fill2(b);
        h2_ba.merge(&fill2(a)).unwrap();
        let mut h2_right = fill2(c);
        h2_right.merge(&h2_ba).unwrap();
        assert_h2_equal(&h2_left, &h2_right);
    }
}

#[test]
fn merging_an_empty_histogram_is_the_identity() {
    let mut rng = SplitMix64::new(77);
    let data = events(&mut rng, 500);
    let mut h = fill1(&data);
    let before = h.clone();
    h.merge(&Histogram1D::uniform(40, -6.0, 6.0).unwrap()).unwrap();
    assert_h1_equal(&h, &before);
}
