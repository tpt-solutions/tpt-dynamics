"""Basic behaviour tests for the tpt_dynamics Python bindings."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import tpt_dynamics as td

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def _lattice_gro(n_side: int = 3, spacing: float = 2.0) -> str:
    """A simple cubic lattice of Ar-like atoms in a cubic box, as GRO text."""
    n = n_side**3
    box = n_side * spacing
    lines = ["lattice test", f"{n:5d}"]
    serial = 0
    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                serial += 1
                x, y, z = (ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing
                lines.append(f"{1:5d}{'AR':<5s}{'AR':>5s}{serial:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append(f"{box:10.5f}{box:10.5f}{box:10.5f}")
    return "\n".join(lines) + "\n"


MD_CONFIG = """
seed = 7
units = "reduced"
dt = 0.001
nsteps = 20
temperature = 0.7
integrator = "velocity-verlet"

[neighbour]
cutoff = 2.5
skin = 0.3

[lennard_jones]
epsilon = 1.0
sigma = 1.0
modifier = "shift"
"""


# --------------------------------------------------------------------------
# Module
# --------------------------------------------------------------------------


def test_version_is_exposed() -> None:
    assert isinstance(td.version(), str)
    assert td.__version__ == td.version()


# --------------------------------------------------------------------------
# Histogram1D
# --------------------------------------------------------------------------


def test_histogram1d_uniform_fill_contents_and_errors() -> None:
    h = td.Histogram1D(4, 0.0, 4.0)
    h.fill(np.array([0.5, 0.5, 1.5, 3.9, -1.0, 4.0, np.nan]))

    np.testing.assert_allclose(h.contents(), [2.0, 1.0, 0.0, 1.0])
    np.testing.assert_allclose(h.errors(), [math.sqrt(2), 1.0, 0.0, 1.0])
    np.testing.assert_allclose(h.edges(), [0.0, 1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(h.bin_centers(), [0.5, 1.5, 2.5, 3.5])
    assert h.underflow == pytest.approx(1.0)
    assert h.overflow == pytest.approx(1.0)
    assert h.nan_fills == 1
    assert h.entries == 6
    assert h.integral() == pytest.approx(4.0)
    assert h.bin(0) == pytest.approx(2.0)
    assert h.bin_error(0) == pytest.approx(math.sqrt(2))


def test_histogram1d_weighted_fill() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    h.fill(np.array([0.5, 0.5, 1.5]), np.array([2.0, 3.0, 0.5]))

    np.testing.assert_allclose(h.contents(), [5.0, 0.5])
    np.testing.assert_allclose(h.errors(), [math.sqrt(13.0), 0.5])


def test_histogram1d_weight_length_mismatch_raises() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    with pytest.raises(ValueError):
        h.fill(np.array([0.5, 1.5]), np.array([1.0]))


def test_histogram1d_variable_binning() -> None:
    h = td.Histogram1D.variable([0.0, 1.0, 10.0])
    h.fill(np.array([0.5, 5.0, 5.0]))
    np.testing.assert_allclose(h.contents(), [1.0, 2.0])


def test_histogram1d_invalid_binning_raises() -> None:
    with pytest.raises(ValueError):
        td.Histogram1D(0, 0.0, 1.0)
    with pytest.raises(ValueError):
        td.Histogram1D.variable([1.0, 0.5])


def test_histogram1d_mean_and_quantile() -> None:
    h = td.Histogram1D(10, 0.0, 10.0)
    h.fill(np.array([2.5, 2.5, 7.5]))
    assert h.mean() == pytest.approx((2.5 + 2.5 + 7.5) / 3.0, rel=0.1)
    assert h.quantile(0.5) is not None
    assert td.Histogram1D(3, 0.0, 1.0).mean() is None


def test_histogram1d_merge_is_additive() -> None:
    a = td.Histogram1D(2, 0.0, 2.0)
    b = td.Histogram1D(2, 0.0, 2.0)
    a.fill(np.array([0.5]))
    b.fill(np.array([0.5, 1.5]))
    a.merge(b)
    np.testing.assert_allclose(a.contents(), [2.0, 1.0])
    assert a.entries == 3


def test_histogram1d_merge_binning_mismatch_raises() -> None:
    a = td.Histogram1D(2, 0.0, 2.0)
    b = td.Histogram1D(3, 0.0, 2.0)
    with pytest.raises(ValueError):
        a.merge(b)


def test_histogram1d_bin_index_out_of_range_raises() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    with pytest.raises(IndexError):
        h.bin(2)


def test_histogram1d_to_dict_and_copy() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    h.fill(np.array([0.5]))
    d = h.to_dict()
    assert set(d) >= {"edges", "contents", "errors", "underflow", "overflow", "entries"}
    np.testing.assert_allclose(d["contents"], [1.0, 0.0])

    c = h.copy()
    c.fill(np.array([1.5]))
    assert h.integral() == pytest.approx(1.0)
    assert c.integral() == pytest.approx(2.0)


# --------------------------------------------------------------------------
# Histogram2D
# --------------------------------------------------------------------------


def test_histogram2d_fill_and_projections() -> None:
    h = td.Histogram2D(2, 0.0, 2.0, 3, 0.0, 3.0)
    h.fill(np.array([0.5, 0.5, 1.5]), np.array([0.5, 2.5, 1.5]))

    contents = h.contents()
    assert contents.shape == (2, 3)
    assert contents[0, 0] == pytest.approx(1.0)
    assert contents[0, 2] == pytest.approx(1.0)
    assert contents[1, 1] == pytest.approx(1.0)
    assert h.integral() == pytest.approx(3.0)

    np.testing.assert_allclose(h.project_x().contents(), [2.0, 1.0])
    np.testing.assert_allclose(h.project_y().contents(), [1.0, 1.0, 1.0])


def test_histogram2d_merge_and_index_errors() -> None:
    a = td.Histogram2D(2, 0.0, 2.0, 2, 0.0, 2.0)
    b = td.Histogram2D(2, 0.0, 2.0, 2, 0.0, 2.0)
    a.fill(np.array([0.5]), np.array([0.5]))
    b.fill(np.array([0.5]), np.array([0.5]))
    a.merge(b)
    assert a.bin(0, 0) == pytest.approx(2.0)
    with pytest.raises(IndexError):
        a.bin(2, 0)
    with pytest.raises(ValueError):
        a.merge(td.Histogram2D(3, 0.0, 2.0, 2, 0.0, 2.0))


# --------------------------------------------------------------------------
# Lazy frames
# --------------------------------------------------------------------------


def test_lazy_scan_dict_filter_with_column_group_by_agg() -> None:
    pt = np.array([10.0, 30.0, 50.0, 25.0])
    g = np.array([1, 2, 2, 1], dtype=np.int64)
    lf = (
        td.LazyFrame.scan_dict({"pt": pt, "g": g})
        .filter(td.col("pt") > 20.0)
        .with_column((td.col("pt") * 2.0).alias("pt2"))
        .group_by(["g"])
        .agg([td.Agg.sum(td.col("pt")).alias("pt_sum"), td.Agg.count().alias("n")])
    )
    out = lf.collect()

    order = np.argsort(out["g"])
    np.testing.assert_array_equal(out["g"][order], [1, 2])
    np.testing.assert_allclose(out["pt_sum"][order], [25.0, 80.0])
    np.testing.assert_array_equal(out["n"][order], [1, 2])


def test_lazy_select_and_expression_result_matches_numpy() -> None:
    px = np.array([3.0, 0.0, -1.0])
    py = np.array([4.0, 2.0, 1.0])
    out = (
        td.LazyFrame.scan_dict({"px": px, "py": py})
        .select([td.transverse_momentum(td.col("px"), td.col("py")).alias("pt")])
        .collect()
    )
    assert list(out) == ["pt"]
    np.testing.assert_allclose(out["pt"], np.hypot(px, py))


def test_lazy_bool_filter_and_chunk_size() -> None:
    x = np.arange(10, dtype=np.float64)
    out = (
        td.LazyFrame.scan_dict({"x": x})
        .with_chunk_size(3)
        .filter((td.col("x") >= 2.0) & (td.col("x") < 5.0))
        .collect()
    )
    np.testing.assert_allclose(out["x"], [2.0, 3.0, 4.0])


def test_lazy_unknown_column_raises() -> None:
    lf = td.LazyFrame.scan_dict({"x": np.zeros(3)})
    with pytest.raises(ValueError):
        lf.filter(td.col("missing") > 0.0).collect()


def test_lazy_rejects_non_numpy_columns() -> None:
    with pytest.raises(TypeError):
        td.LazyFrame.scan_dict({"x": [1.0, 2.0]})


def test_lazy_scan_csv(tmp_path: Path) -> None:
    path = tmp_path / "events.csv"
    path.write_text("pt,n,isolated\n10.0,1,true\n30.0,2,false\n50.0,2,true\n", encoding="utf-8")
    out = (
        td.LazyFrame.scan_csv(str(path), dtypes={"n": "i64", "isolated": "bool"})
        .filter(td.col("isolated").eq_to(True))
        .collect()
    )
    np.testing.assert_allclose(out["pt"], [10.0, 50.0])
    assert out["n"].dtype == np.int64
    assert out["isolated"].dtype == np.bool_


def test_lazy_scan_csv_bad_value_reports_line(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("x\n1.0\nnot-a-number\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 3"):
        td.LazyFrame.scan_csv(str(path)).collect()


def test_expr_repr_and_explain_are_strings() -> None:
    e = td.col("pt") > 1.0
    assert "pt" in repr(e)
    assert isinstance(td.LazyFrame.scan_dict({"pt": np.ones(2)}).explain(), str)


# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------


def test_validate_config_accepts_and_rejects() -> None:
    td.validate_config(MD_CONFIG)
    with pytest.raises(ValueError):
        td.validate_config("nsteps = 10\nunknown_key = 1\n")


def test_run_simulation_returns_energy_arrays() -> None:
    out = td.run_simulation(MD_CONFIG, _lattice_gro(), format="gro")

    expected_keys = {"step", "time", "kinetic", "potential", "total", "temperature"}
    assert set(out) == expected_keys
    np.testing.assert_array_equal(out["step"], np.arange(21, dtype=np.uint64))
    for key in expected_keys - {"step"}:
        assert out[key].shape == (21,), key
        assert np.all(np.isfinite(out[key])), key
    np.testing.assert_allclose(out["total"], out["kinetic"] + out["potential"], rtol=1e-9, atol=1e-9)
    assert out["time"][-1] == pytest.approx(20 * 0.001)


def test_run_simulation_is_deterministic_and_sampling_works() -> None:
    gro = _lattice_gro()
    a = td.run_simulation(MD_CONFIG, gro)
    b = td.run_simulation(MD_CONFIG, gro)
    np.testing.assert_array_equal(a["total"], b["total"])

    sampled = td.run_simulation(MD_CONFIG, gro, sample_every=5)
    np.testing.assert_array_equal(sampled["step"], [0, 5, 10, 15, 20])


def test_run_simulation_errors() -> None:
    with pytest.raises(ValueError):
        td.run_simulation(MD_CONFIG, _lattice_gro(), format="xyz")
    with pytest.raises(ValueError):
        td.run_simulation(MD_CONFIG, _lattice_gro(), sample_every=0)
    with pytest.raises(ValueError):
        td.run_simulation("dt = -1.0\n", _lattice_gro())
