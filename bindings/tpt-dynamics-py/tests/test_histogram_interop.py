"""Hist-style numpy interop for Histogram1D and Histogram2D."""

from __future__ import annotations

import numpy as np
import pytest

import tpt_dynamics as td

# --------------------------------------------------------------------------
# Axes
# --------------------------------------------------------------------------


def test_axis_1d_exposes_edges_centers_and_length() -> None:
    h = td.Histogram1D.variable([0.0, 1.0, 3.0, 7.0])
    (axis,) = h.axes
    assert isinstance(axis, td.Axis)
    np.testing.assert_allclose(axis.edges, [0.0, 1.0, 3.0, 7.0])
    np.testing.assert_allclose(axis.centers, [0.5, 2.0, 5.0])
    assert len(axis) == 3
    assert axis.edges.dtype == np.float64
    assert "nbins=3" in repr(axis)


def test_axes_2d_are_x_then_y() -> None:
    h = td.Histogram2D.variable([0.0, 2.0, 4.0], [-1.0, 0.0, 1.0, 2.0])
    xa, ya = h.axes
    np.testing.assert_allclose(xa.edges, [0.0, 2.0, 4.0])
    np.testing.assert_allclose(ya.edges, [-1.0, 0.0, 1.0, 2.0])
    assert len(xa) == 2
    assert len(ya) == 3


def test_axis_edges_are_copies() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    h.axes[0].edges[0] = 99.0
    assert h.edges()[0] == pytest.approx(0.0)


# --------------------------------------------------------------------------
# values / variances / to_numpy
# --------------------------------------------------------------------------


def test_1d_values_variances_and_to_numpy() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    h.fill(np.array([0.5, 0.5, 1.5]), weights=np.array([2.0, 3.0, 4.0]))

    np.testing.assert_allclose(h.values(), [5.0, 4.0])
    np.testing.assert_allclose(h.variances(), [4.0 + 9.0, 16.0])
    np.testing.assert_allclose(h.errors() ** 2, h.variances())

    values, edges = h.to_numpy()
    np.testing.assert_allclose(values, h.values())
    np.testing.assert_allclose(edges, [0.0, 1.0, 2.0])
    assert values.dtype == np.float64


def test_1d_values_are_copies() -> None:
    h = td.Histogram1D(2, 0.0, 2.0)
    h.fill(np.array([0.5]))
    h.values()[0] = 100.0
    assert h.bin(0) == pytest.approx(1.0)


def test_2d_values_variances_and_to_numpy_shapes() -> None:
    h = td.Histogram2D(2, 0.0, 2.0, 3, 0.0, 3.0)
    h.fill(np.array([0.5, 0.5, 1.5]), np.array([0.5, 0.5, 2.5]), np.array([1.0, 2.0, 3.0]))

    values = h.values()
    assert values.shape == (2, 3)
    assert values[0, 0] == pytest.approx(3.0)
    assert values[1, 2] == pytest.approx(3.0)
    assert h.variances()[0, 0] == pytest.approx(1.0 + 4.0)
    assert h.variances()[1, 2] == pytest.approx(9.0)

    v, xe, ye = h.to_numpy()
    assert v.shape == (2, 3)
    np.testing.assert_allclose(xe, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(ye, [0.0, 1.0, 2.0, 3.0])


# --------------------------------------------------------------------------
# from_numpy
# --------------------------------------------------------------------------


def test_1d_round_trip_through_to_numpy() -> None:
    h = td.Histogram1D.variable([0.0, 1.0, 2.0, 4.0])
    h.fill(np.array([0.2, 1.5, 1.7, 3.0]), weights=np.array([1.0, 0.5, 2.0, 1.5]))

    values, edges = h.to_numpy()
    copy = td.Histogram1D.from_numpy(values, edges, variances=h.variances())
    np.testing.assert_allclose(copy.values(), h.values())
    np.testing.assert_allclose(copy.variances(), h.variances())
    np.testing.assert_allclose(copy.edges(), h.edges())
    assert copy.nbins() == h.nbins()


def test_1d_from_numpy_defaults_variances_to_values() -> None:
    values = np.array([4.0, 9.0])
    h = td.Histogram1D.from_numpy(values, np.array([0.0, 1.0, 2.0]))
    np.testing.assert_allclose(h.variances(), values)
    np.testing.assert_allclose(h.errors(), [2.0, 3.0])


def test_1d_imported_histogram_has_no_moments_and_accepts_fills() -> None:
    h = td.Histogram1D.from_numpy(np.array([1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
    assert h.entries == 0
    assert h.underflow == 0.0
    assert h.mean() is None
    h.fill(np.array([0.5]))
    assert h.bin(0) == pytest.approx(2.0)
    assert h.entries == 1


@pytest.mark.parametrize(
    ("values", "edges", "variances"),
    [
        ([1.0, 2.0], [0.0, 1.0], None),  # length mismatch
        ([1.0], [0.0, 1.0], [-1.0]),  # negative variance
        ([np.nan], [0.0, 1.0], None),  # non-finite content
        ([1.0], [1.0, 0.5], None),  # decreasing edges
        ([1.0], [0.0], None),  # no bins
    ],
)
def test_1d_from_numpy_rejects_bad_input(
    values: list[float], edges: list[float], variances: list[float] | None
) -> None:
    with pytest.raises(ValueError):
        td.Histogram1D.from_numpy(
            np.asarray(values, dtype=np.float64),
            np.asarray(edges, dtype=np.float64),
            None if variances is None else np.asarray(variances, dtype=np.float64),
        )


def test_2d_round_trip_through_to_numpy() -> None:
    h = td.Histogram2D.variable([0.0, 1.0, 3.0], [0.0, 2.0, 4.0])
    h.fill(np.array([0.5, 2.0, 2.5]), np.array([0.5, 3.0, 3.0]), np.array([1.0, 2.0, 0.5]))

    values, xe, ye = h.to_numpy()
    copy = td.Histogram2D.from_numpy(values, xe, ye, variances=h.variances())
    np.testing.assert_allclose(copy.values(), h.values())
    np.testing.assert_allclose(copy.variances(), h.variances())
    np.testing.assert_allclose(copy.x_edges(), h.x_edges())
    np.testing.assert_allclose(copy.y_edges(), h.y_edges())
    assert copy.bin(1, 1) == pytest.approx(h.bin(1, 1))


def test_2d_from_numpy_defaults_variances_and_checks_shape() -> None:
    values = np.arange(6, dtype=np.float64).reshape(2, 3)
    h = td.Histogram2D.from_numpy(values, np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0, 3.0]))
    np.testing.assert_allclose(h.variances(), values)
    assert h.bin(1, 2) == pytest.approx(5.0)

    with pytest.raises(ValueError, match="shape"):
        td.Histogram2D.from_numpy(
            values.T, np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0, 3.0])
        )


def test_2d_from_numpy_rejects_negative_variance() -> None:
    with pytest.raises(ValueError):
        td.Histogram2D.from_numpy(
            np.ones((1, 1)),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            variances=np.array([[-1.0]]),
        )
