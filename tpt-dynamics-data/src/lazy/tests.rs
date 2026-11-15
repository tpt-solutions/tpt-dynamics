//! Tests for the lazy data frame: correctness against eager Rust, optimiser
//! equivalence and plan shape, memory bounds, and validation errors.

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::*;

/// Deterministic uniform numbers in [0, 1) without an external crate.
fn lcg(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

/// Two-particle event table: kinematics of two objects, a group key and a flag.
fn synthetic(n: usize, seed: u64) -> DataFrame {
    let mut s = seed;
    let mut pt1 = Vec::with_capacity(n);
    let mut eta1 = Vec::with_capacity(n);
    let mut phi1 = Vec::with_capacity(n);
    let mut m1 = Vec::with_capacity(n);
    let mut pt2 = Vec::with_capacity(n);
    let mut eta2 = Vec::with_capacity(n);
    let mut phi2 = Vec::with_capacity(n);
    let mut m2 = Vec::with_capacity(n);
    let mut g = Vec::with_capacity(n);
    let mut flag = Vec::with_capacity(n);
    for _ in 0..n {
        pt1.push(5.0 + 95.0 * lcg(&mut s));
        eta1.push(-2.5 + 5.0 * lcg(&mut s));
        phi1.push(-PI + TAU * lcg(&mut s));
        m1.push(0.5 * lcg(&mut s));
        pt2.push(5.0 + 95.0 * lcg(&mut s));
        eta2.push(-2.5 + 5.0 * lcg(&mut s));
        phi2.push(-PI + TAU * lcg(&mut s));
        m2.push(0.5 * lcg(&mut s));
        g.push((lcg(&mut s) * 5.0) as i64);
        flag.push(lcg(&mut s) < 0.7);
    }
    DataFrame::new(vec![
        ("pt1", Column::F64(pt1)),
        ("eta1", Column::F64(eta1)),
        ("phi1", Column::F64(phi1)),
        ("m1", Column::F64(m1)),
        ("pt2", Column::F64(pt2)),
        ("eta2", Column::F64(eta2)),
        ("phi2", Column::F64(phi2)),
        ("m2", Column::F64(m2)),
        ("g", Column::I64(g)),
        ("flag", Column::Bool(flag)),
    ])
    .unwrap()
}

fn f64_col<'a>(df: &'a DataFrame, name: &str) -> &'a [f64] {
    match df.column(name).unwrap() {
        Column::F64(v) => v,
        other => panic!("{name} is {}", other.dtype()),
    }
}

fn i64_col<'a>(df: &'a DataFrame, name: &str) -> &'a [i64] {
    match df.column(name).unwrap() {
        Column::I64(v) => v,
        other => panic!("{name} is {}", other.dtype()),
    }
}

/// Compares names, row counts and values. Floats must agree to `tol` relative.
fn assert_frames_close(a: &DataFrame, b: &DataFrame, tol: f64) {
    assert_eq!(a.names(), b.names());
    assert_eq!(a.len(), b.len());
    for name in a.names() {
        match (a.column(name).unwrap(), b.column(name).unwrap()) {
            (Column::F64(x), Column::F64(y)) => {
                for (i, (p, q)) in x.iter().zip(y).enumerate() {
                    let scale = p.abs().max(q.abs()).max(1.0);
                    assert!((p - q).abs() <= tol * scale, "{name}[{i}]: {p} vs {q}");
                }
            }
            (x, y) => assert_eq!(x, y, "column {name}"),
        }
    }
}

/// Invariant mass computed from E^2 = p^2 + m^2 and M^2 = m1^2 + m2^2 + 2(E1E2 - p1.p2).
/// Deliberately a different formula from the library kernel.
fn reference_mass(r: [f64; 8]) -> f64 {
    let p = |pt: f64, eta: f64, phi: f64| [pt * phi.cos(), pt * phi.sin(), pt * eta.sinh()];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let p1 = p(r[0], r[1], r[2]);
    let p2 = p(r[4], r[5], r[6]);
    let e1 = (dot(p1, p1) + r[3] * r[3]).sqrt();
    let e2 = (dot(p2, p2) + r[7] * r[7]).sqrt();
    let m_sq = r[3] * r[3] + r[7] * r[7] + 2.0 * (e1 * e2 - dot(p1, p2));
    m_sq.max(0.0).sqrt()
}

/// The analysis used throughout: kinematic filter, mass per event, then
/// per-group aggregates.
fn analysis(df: &DataFrame) -> LazyFrame {
    LazyFrame::scan(df.clone())
        .filter(col("pt1").gt(20.0).and(col("flag")))
        .with_column(
            invariant_mass(
                col("pt1"),
                col("eta1"),
                col("phi1"),
                col("m1"),
                col("pt2"),
                col("eta2"),
                col("phi2"),
                col("m2"),
            )
            .alias("mass"),
        )
        .group_by(["g"])
        .agg([
            AggExpr::count(),
            AggExpr::sum(col("mass")).alias("mass_sum"),
            AggExpr::mean(col("pt1")).alias("pt1_mean"),
            AggExpr::min(col("pt2")).alias("pt2_min"),
            AggExpr::max(col("pt2")).alias("pt2_max"),
            AggExpr::sum(col("g")).alias("g_sum"),
        ])
}

/// Eager reference for [`analysis`], written as plain loops over the columns.
fn eager_analysis(df: &DataFrame) -> DataFrame {
    #[derive(Default)]
    struct Group {
        count: i64,
        mass_sum: f64,
        pt1_sum: f64,
        pt2_min: f64,
        pt2_max: f64,
        g_sum: i64,
    }
    let pt1 = f64_col(df, "pt1");
    let eta1 = f64_col(df, "eta1");
    let phi1 = f64_col(df, "phi1");
    let m1 = f64_col(df, "m1");
    let pt2 = f64_col(df, "pt2");
    let eta2 = f64_col(df, "eta2");
    let phi2 = f64_col(df, "phi2");
    let m2 = f64_col(df, "m2");
    let g = i64_col(df, "g");
    let Column::Bool(flag) = df.column("flag").unwrap() else { panic!("flag must be bool") };
    let mut groups: BTreeMap<i64, Group> = BTreeMap::new();
    for i in 0..df.len() {
        if !(pt1[i] > 20.0 && flag[i]) {
            continue;
        }
        let mass =
            reference_mass([pt1[i], eta1[i], phi1[i], m1[i], pt2[i], eta2[i], phi2[i], m2[i]]);
        let grp = groups.entry(g[i]).or_insert(Group {
            pt2_min: f64::INFINITY,
            pt2_max: f64::NEG_INFINITY,
            ..Group::default()
        });
        grp.count += 1;
        grp.mass_sum += mass;
        grp.pt1_sum += pt1[i];
        grp.pt2_min = grp.pt2_min.min(pt2[i]);
        grp.pt2_max = grp.pt2_max.max(pt2[i]);
        grp.g_sum += g[i];
    }
    let keys: Vec<i64> = groups.keys().copied().collect();
    let vals: Vec<&Group> = groups.values().collect();
    DataFrame::new(vec![
        ("g", Column::I64(keys)),
        ("count", Column::I64(vals.iter().map(|v| v.count).collect())),
        ("mass_sum", Column::F64(vals.iter().map(|v| v.mass_sum).collect())),
        ("pt1_mean", Column::F64(vals.iter().map(|v| v.pt1_sum / v.count as f64).collect())),
        ("pt2_min", Column::F64(vals.iter().map(|v| v.pt2_min).collect())),
        ("pt2_max", Column::F64(vals.iter().map(|v| v.pt2_max).collect())),
        ("g_sum", Column::I64(vals.iter().map(|v| v.g_sum).collect())),
    ])
    .unwrap()
}

#[test]
fn small_pipeline_gives_expected_values() {
    let df = DataFrame::new(vec![
        ("pt", Column::F64(vec![10.0, 30.0, 50.0, 40.0])),
        ("g", Column::I64(vec![1, 2, 2, 1])),
    ])
    .unwrap();
    let out = LazyFrame::scan(df)
        .filter(col("pt").gt(15.0))
        .with_column((col("pt") * lit(2.0)).alias("pt2"))
        .group_by(["g"])
        .agg([
            AggExpr::count().alias("n"),
            AggExpr::sum(col("pt2")).alias("s"),
            AggExpr::mean(col("pt")).alias("m"),
            AggExpr::min(col("pt")).alias("lo"),
            AggExpr::max(col("pt")).alias("hi"),
        ])
        .collect()
        .unwrap();
    assert_eq!(out.names(), ["g", "n", "s", "m", "lo", "hi"]);
    assert_eq!(i64_col(&out, "g"), [1, 2]);
    assert_eq!(i64_col(&out, "n"), [1, 2]);
    assert_eq!(f64_col(&out, "s"), [80.0, 160.0]);
    assert_eq!(f64_col(&out, "m"), [40.0, 40.0]);
    assert_eq!(f64_col(&out, "lo"), [40.0, 30.0]);
    assert_eq!(f64_col(&out, "hi"), [40.0, 50.0]);
}

#[test]
fn select_and_with_column_replace_outputs() {
    let df = DataFrame::new(vec![
        ("a", Column::I64(vec![1, 2, 3])),
        ("b", Column::F64(vec![0.5, 1.5, 2.5])),
    ])
    .unwrap();
    let out = LazyFrame::scan(df.clone())
        .with_column((col("a") + lit(10_i64)).alias("a"))
        .select([col("a"), (col("b") * lit(2.0)).alias("bb")])
        .collect()
        .unwrap();
    assert_eq!(out.names(), ["a", "bb"]);
    assert_eq!(i64_col(&out, "a"), [11, 12, 13]);
    assert_eq!(f64_col(&out, "bb"), [1.0, 3.0, 5.0]);
}

#[test]
fn lazy_matches_eager_reference_on_one_million_rows() {
    let df = synthetic(1_000_000, 42);
    let lazy = analysis(&df).collect().unwrap();
    let eager = eager_analysis(&df);
    assert_frames_close(&lazy, &eager, 1e-6);
    // Group sizes are integers and must match exactly.
    assert_eq!(i64_col(&lazy, "count"), i64_col(&eager, "count"));

    // Sequential execution with the same chunking is bit-identical.
    let seq = analysis(&df)
        .with_options(CollectOptions { parallel: false, ..CollectOptions::default() })
        .collect()
        .unwrap();
    assert_frames_close(&lazy, &seq, 0.0);
}

#[test]
fn optimised_and_unoptimised_plans_agree() {
    let df = synthetic(20_000, 7);
    let queries: Vec<LazyFrame> = vec![
        analysis(&df),
        LazyFrame::scan(df.clone())
            .with_column((col("pt1") * lit(2.0)).alias("x"))
            .filter(col("g").gt(1_i64))
            .group_by(["g"])
            .agg([AggExpr::sum(col("x")).alias("xs"), AggExpr::count()]),
        LazyFrame::scan(df.clone())
            .select([col("pt1"), col("g")])
            .filter(col("pt1").gt(50.0))
            .select([(col("pt1") / lit(2.0)).alias("half"), col("g")])
            .group_by(["g"])
            .agg([AggExpr::max(col("half")).alias("hmax")]),
        LazyFrame::scan(df.clone())
            .group_by(["g"])
            .agg([AggExpr::count(), AggExpr::mean(col("pt2")).alias("pm")])
            .filter(col("g").ge(2_i64)),
        LazyFrame::scan(df.clone())
            .with_column(col("pt1").alias("g"))
            .filter(col("flag"))
            .select([col("g").alias("v")]),
    ];
    for q in queries {
        let plan = q.clone();
        let opt = q.with_chunk_size(3_000).collect().unwrap();
        let raw = plan.with_chunk_size(3_000).collect_unoptimized().unwrap();
        assert_frames_close(&opt, &raw, 0.0);
    }
}

#[test]
fn result_is_independent_of_chunk_size_and_threads() {
    let df = synthetic(5_000, 3);
    let base = analysis(&df)
        .with_options(CollectOptions {
            chunk_size: usize::MAX,
            parallel: false,
            ..CollectOptions::default()
        })
        .collect()
        .unwrap();
    for chunk in [1_usize, 7, 1_000] {
        for parallel in [false, true] {
            let out = analysis(&df)
                .with_options(CollectOptions { chunk_size: chunk, window: 3, parallel })
                .collect()
                .unwrap();
            // Integer columns are exact; float sums differ only by rounding.
            assert_eq!(i64_col(&out, "count"), i64_col(&base, "count"));
            assert_eq!(i64_col(&out, "g_sum"), i64_col(&base, "g_sum"));
            assert_frames_close(&out, &base, 1e-9);
        }
    }
}

#[test]
fn explain_shows_predicate_pushdown_and_projection_pushdown() {
    let df = synthetic(100, 1);
    // The filter reads only `pt1`, so it sinks below the with-column and the
    // unused with-column is dropped. The scan then reads just pt1 and g.
    let lf = LazyFrame::scan(df)
        .with_column((col("m1") * lit(2.0)).alias("m2x"))
        .filter(col("pt1").gt(20.0))
        .group_by(["g"])
        .agg([AggExpr::sum(col("pt1")).alias("pt_sum")]);
    let raw = lf.explain_unoptimized();
    assert!(raw.contains("WithColumn (m1 * 2) AS m2x"), "{raw}");
    assert!(raw.contains("Scan [*]"), "{raw}");

    let lines: Vec<String> = lf.explain().lines().map(|l| l.trim().to_string()).collect();
    assert_eq!(
        lines,
        ["Aggregate by=[g] aggs=[pt_sum=sum(pt1)]", "Filter (pt1 > 20)", "Scan [pt1, g] rows=100",]
    );
}

#[test]
fn filter_on_group_key_moves_below_aggregate() {
    let df = synthetic(100, 2);
    let lf = LazyFrame::scan(df).group_by(["g"]).agg([AggExpr::count()]).filter(col("g").gt(1_i64));
    let lines: Vec<String> = lf.explain().lines().map(|l| l.trim().to_string()).collect();
    assert_eq!(
        lines,
        ["Aggregate by=[g] aggs=[count=count()]", "Filter (g > 1)", "Scan [g] rows=100",]
    );
}

#[test]
fn filter_on_computed_column_stays_above_select() {
    let df = synthetic(100, 4);
    let lf =
        LazyFrame::scan(df).select([(col("pt1") * lit(2.0)).alias("x")]).filter(col("x").gt(1.0));
    let lines: Vec<String> = lf.explain().lines().map(|l| l.trim().to_string()).collect();
    assert_eq!(lines, ["Filter (x > 1)", "Select [(pt1 * 2) AS x]", "Scan [pt1] rows=100",]);
}

#[test]
fn filter_through_select_of_identity_columns() {
    let df = synthetic(100, 5);
    let lf =
        LazyFrame::scan(df).select([col("pt1"), col("g").alias("grp")]).filter(col("pt1").gt(10.0));
    let lines: Vec<String> = lf.explain().lines().map(|l| l.trim().to_string()).collect();
    assert_eq!(lines, ["Select [pt1, g AS grp]", "Filter (pt1 > 10)", "Scan [pt1, g] rows=100",]);
}

#[test]
fn physics_functions_match_known_values() {
    // Back-to-back massless objects: invariant mass is 2 * pt.
    let df = DataFrame::new(vec![
        ("pt1", Column::F64(vec![50.0, 10.0])),
        ("eta1", Column::F64(vec![0.0, 1.0])),
        ("phi1", Column::F64(vec![0.0, 0.5])),
        ("m1", Column::F64(vec![0.0, 0.105])),
        ("pt2", Column::F64(vec![50.0, 20.0])),
        ("eta2", Column::F64(vec![0.0, -0.5])),
        ("phi2", Column::F64(vec![PI, 2.0])),
        ("m2", Column::F64(vec![0.0, 0.0])),
        ("eta3", Column::F64(vec![0.0, 3.0])),
        ("phi3", Column::F64(vec![0.0, 6.0])),
        ("px", Column::F64(vec![3.0, 0.0])),
        ("py", Column::F64(vec![4.0, 0.0])),
    ])
    .unwrap();
    let out = LazyFrame::scan(df.clone())
        .select([
            invariant_mass(
                col("pt1"),
                col("eta1"),
                col("phi1"),
                col("m1"),
                col("pt2"),
                col("eta2"),
                col("phi2"),
                col("m2"),
            )
            .alias("mll"),
            delta_r(col("eta1"), col("phi1"), col("eta3"), col("phi3")).alias("dr"),
            transverse_momentum(col("px"), col("py")).alias("pt"),
        ])
        .collect()
        .unwrap();
    let mll = f64_col(&out, "mll");
    assert!((mll[0] - 100.0).abs() < 1e-9, "{}", mll[0]);
    let expect = reference_mass([10.0, 1.0, 0.5, 0.105, 20.0, -0.5, 2.0, 0.0]);
    assert!((mll[1] - expect).abs() < 1e-9, "{} vs {expect}", mll[1]);

    let dr = f64_col(&out, "dr");
    assert!((dr[0]).abs() < 1e-12);
    // Row 1: deta = 1 - 3 = -2, dphi = 0.5 - 6 wraps into [-pi, pi).
    let expect_dr = (4.0_f64 + (0.5 - 6.0 + TAU).powi(2)).sqrt();
    assert!((dr[1] - expect_dr).abs() < 1e-12, "{}", dr[1]);

    let pt = f64_col(&out, "pt");
    assert!((pt[0] - 5.0).abs() < 1e-12);
    assert!(pt[1].abs() < 1e-12);
}

#[test]
fn memory_is_bounded_by_the_window() {
    /// A source that produces chunks on demand. Each chunk holds a lease that
    /// counts how many chunks are alive at once.
    struct Counting {
        rows: usize,
        live: Arc<AtomicUsize>,
        peak: Arc<AtomicUsize>,
    }

    struct Alive {
        live: Arc<AtomicUsize>,
    }

    impl Alive {
        fn new(live: &Arc<AtomicUsize>, peak: &Arc<AtomicUsize>) -> Self {
            let now = live.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            Self { live: Arc::clone(live) }
        }
    }

    impl Drop for Alive {
        fn drop(&mut self) {
            self.live.fetch_sub(1, Ordering::SeqCst);
        }
    }

    impl Source for Counting {
        fn schema(&self) -> Schema {
            vec![
                Field { name: "x".into(), dtype: DType::F64 },
                Field { name: "g".into(), dtype: DType::I64 },
            ]
        }

        fn num_rows(&self) -> usize {
            self.rows
        }

        fn read(&self, rows: Range<usize>, projection: &[String]) -> Result<DataFrame, LazyError> {
            let end = rows.end.min(self.rows);
            let start = rows.start.min(end);
            let mut cols = Vec::new();
            for name in projection {
                match name.as_str() {
                    "x" => cols.push((
                        "x",
                        Column::F64((start..end).map(|i| (i % 100) as f64 - 50.0).collect()),
                    )),
                    "g" => cols.push((
                        "g",
                        Column::I64((start..end).map(|i| i64::from((i % 4) as u8)).collect()),
                    )),
                    other => return Err(LazyError::UnknownColumn(other.into())),
                }
            }
            Ok(DataFrame::with_rows(end - start, cols)?
                .with_lease(Alive::new(&self.live, &self.peak)))
        }
    }

    let rows = 100_000;
    let chunk = 1_000;
    let window = 4;
    let expected_pass = (0..rows).filter(|i| (i % 100) as f64 - 50.0 > 0.0).count();
    for parallel in [false, true] {
        let live = Arc::new(AtomicUsize::new(0));
        let peak = Arc::new(AtomicUsize::new(0));
        let source = Counting { rows, live: Arc::clone(&live), peak: Arc::clone(&peak) };
        let out = LazyFrame::scan(source)
            .filter(col("x").gt(0.0))
            .group_by(["g"])
            .agg([AggExpr::count()])
            .with_options(CollectOptions { chunk_size: chunk, window, parallel })
            .collect()
            .unwrap();
        let total: usize = i64_col(&out, "count").iter().map(|&c| c as usize).sum();
        assert_eq!(total, expected_pass);
        let peak = peak.load(Ordering::SeqCst);
        assert!((1..=window).contains(&peak), "peak {peak} chunks alive, window {window}");
        assert_eq!(live.load(Ordering::SeqCst), 0, "all chunks released");
    }
}

#[test]
fn empty_source_gives_typed_empty_result() {
    let df = DataFrame::new(vec![("g", Column::I64(vec![])), ("x", Column::F64(vec![]))]).unwrap();
    let out = LazyFrame::scan(df)
        .filter(col("x").gt(0.0))
        .group_by(["g"])
        .agg([AggExpr::count(), AggExpr::sum(col("x")).alias("xs")])
        .collect()
        .unwrap();
    assert_eq!(out.len(), 0);
    assert_eq!(out.names(), ["g", "count", "xs"]);
    assert_eq!(out.schema()[2].dtype, DType::F64);
    assert_eq!(out.schema()[1].dtype, DType::I64);
}

#[test]
fn validation_reports_clear_errors() {
    let df = DataFrame::new(vec![
        ("x", Column::F64(vec![1.0, 2.0])),
        ("g", Column::I64(vec![1, 2])),
        ("f", Column::F64(vec![1.0, 2.0])),
    ])
    .unwrap();
    let scan = || LazyFrame::scan(df.clone());

    let err = scan().filter(col("nope").gt(0.0)).collect().unwrap_err();
    assert_eq!(err, LazyError::UnknownColumn("nope".into()));

    let err = scan().filter(col("x") + lit(1.0)).collect().unwrap_err();
    assert!(matches!(err, LazyError::Type(_)), "{err}");

    let err = scan().group_by(["f"]).agg([AggExpr::count()]).collect().unwrap_err();
    assert!(matches!(err, LazyError::Type(_)), "{err}");

    let err = scan().select([col("x") * lit(2.0)]).collect().unwrap_err();
    assert!(matches!(err, LazyError::Invalid(_)), "{err}");

    let err = scan().group_by(["g"]).agg([AggExpr::count().alias("g")]).collect().unwrap_err();
    assert!(matches!(err, LazyError::DuplicateColumn(_)), "{err}");

    let err = scan().group_by(Vec::<String>::new()).agg([AggExpr::count()]).collect().unwrap_err();
    assert!(matches!(err, LazyError::Invalid(_)), "{err}");

    let err = scan().group_by(["g"]).agg([AggExpr::sum(col("g").gt(0_i64))]).collect().unwrap_err();
    assert!(matches!(err, LazyError::Type(_)), "{err}");

    let err = scan().with_chunk_size(0).collect().unwrap_err();
    assert!(matches!(err, LazyError::Invalid(_)), "{err}");
}

#[test]
fn ragged_columns_are_rejected() {
    let err = DataFrame::new(vec![("a", Column::I64(vec![1, 2])), ("b", Column::I64(vec![1]))])
        .unwrap_err();
    assert!(matches!(err, LazyError::Invalid(_)), "{err}");
}

#[test]
fn constant_folding_collapses_literal_arithmetic() {
    use super::expr::fold_constants;
    // (2 * 3) + 1 folds to 7.0; the column reference is untouched.
    let folded = fold_constants((lit(2.0) * lit(3.0)) + lit(1.0));
    assert_eq!(folded, lit(7.0));
    let mixed = fold_constants(col("pt1") * (lit(2_i64) + lit(3_i64)));
    assert_eq!(mixed, col("pt1") * lit(5_i64));
    let neg = fold_constants(-lit(4.0));
    assert_eq!(neg, lit(-4.0));
}

#[test]
fn constant_folding_leaves_type_errors_for_runtime() {
    use super::expr::fold_constants;
    // Bool + number is a type error; folding must not hide it.
    let bad = lit(true) + lit(1.0);
    assert_eq!(fold_constants(bad.clone()), bad);
}

#[test]
fn folded_plans_match_unfolded_results_bit_for_bit() {
    let df = synthetic(2_000, 9);
    let lf = LazyFrame::scan(df.clone())
        .with_column((col("pt1") * (lit(2.0) * lit(0.5) + lit(0.25))).alias("scaled"))
        .filter(col("scaled").gt(lit(10.0) / lit(4.0)))
        .select([col("scaled"), col("g")]);
    let explained = lf.explain();
    let opt = lf.clone().collect().unwrap();
    let raw = lf.collect_unoptimized().unwrap();
    assert_frames_close(&opt, &raw, 0.0);
    assert!(explained.contains("2.5"), "{explained}");
}
