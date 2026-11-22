"""Type stubs for the native tpt-dynamics extension module."""

from collections.abc import Mapping, Sequence
from typing import Final, Literal, Union

import numpy as np
from numpy.typing import NDArray

__version__: Final[str]

_Operand = Union["Expr", int, float, bool]
_ColumnArray = NDArray[np.float64] | NDArray[np.int64] | NDArray[np.bool_]

def version() -> str:
    """Returns the crate version."""
    ...

class Expr:
    """An expression over the columns of a frame. Build with `col` and `lit`.

    Comparison methods are used instead of `==` and `!=` so that Python
    equality on expression objects keeps its usual meaning.
    """

    def alias(self, name: str) -> Expr: ...
    def gt(self, other: _Operand) -> Expr: ...
    def ge(self, other: _Operand) -> Expr: ...
    def lt(self, other: _Operand) -> Expr: ...
    def le(self, other: _Operand) -> Expr: ...
    def eq_to(self, other: _Operand) -> Expr: ...
    def ne_to(self, other: _Operand) -> Expr: ...
    def __gt__(self, other: _Operand) -> Expr: ...
    def __ge__(self, other: _Operand) -> Expr: ...
    def __lt__(self, other: _Operand) -> Expr: ...
    def __le__(self, other: _Operand) -> Expr: ...
    def __add__(self, other: _Operand) -> Expr: ...
    def __radd__(self, other: _Operand) -> Expr: ...
    def __sub__(self, other: _Operand) -> Expr: ...
    def __rsub__(self, other: _Operand) -> Expr: ...
    def __mul__(self, other: _Operand) -> Expr: ...
    def __rmul__(self, other: _Operand) -> Expr: ...
    def __truediv__(self, other: _Operand) -> Expr: ...
    def __rtruediv__(self, other: _Operand) -> Expr: ...
    def __and__(self, other: Expr) -> Expr: ...
    def __or__(self, other: Expr) -> Expr: ...
    def __neg__(self) -> Expr: ...
    def __invert__(self) -> Expr: ...
    def __repr__(self) -> str: ...

class Agg:
    """An aggregate for `LazyFrame.group_by(...).agg(...)`."""

    @staticmethod
    def count() -> Agg: ...
    @staticmethod
    def sum(expr: _Operand) -> Agg: ...
    @staticmethod
    def mean(expr: _Operand) -> Agg: ...
    @staticmethod
    def min(expr: _Operand) -> Agg: ...
    @staticmethod
    def max(expr: _Operand) -> Agg: ...
    def alias(self, name: str) -> Agg: ...

def col(name: str) -> Expr:
    """References a column by name."""
    ...

def lit(value: int | float | bool) -> Expr:
    """A constant broadcast to every row."""
    ...

def invariant_mass(
    pt1: _Operand,
    eta1: _Operand,
    phi1: _Operand,
    m1: _Operand,
    pt2: _Operand,
    eta2: _Operand,
    phi2: _Operand,
    m2: _Operand,
) -> Expr:
    """Invariant mass of two particles, each given as (pt, eta, phi, m)."""
    ...

def delta_r(eta1: _Operand, phi1: _Operand, eta2: _Operand, phi2: _Operand) -> Expr:
    """Angular distance sqrt(d_eta^2 + d_phi^2), with d_phi wrapped to [-pi, pi)."""
    ...

def transverse_momentum(px: _Operand, py: _Operand) -> Expr:
    """Transverse momentum sqrt(px^2 + py^2)."""
    ...

class LazyFrame:
    """A deferred pipeline over a data source. Nothing runs until `collect()`."""

    @staticmethod
    def scan_dict(columns: Mapping[str, _ColumnArray]) -> LazyFrame:
        """Scans in-memory 1D numpy arrays (float64, int64 or bool) of equal length."""
        ...

    @staticmethod
    def scan_csv(
        path: str,
        dtypes: Mapping[str, Literal["f64", "i64", "bool"]] | None = None,
        delimiter: str = ",",
    ) -> LazyFrame:
        """Scans a CSV file with a header line.

        The file is parsed in full when the frame is created. Columns not listed
        in `dtypes` are float64.
        """
        ...

    def filter(self, predicate: Expr) -> LazyFrame: ...
    def select(self, exprs: Sequence[Expr]) -> LazyFrame: ...
    def with_column(self, expr: Expr) -> LazyFrame: ...
    def with_chunk_size(self, chunk_size: int) -> LazyFrame: ...
    def group_by(self, keys: Sequence[str]) -> GroupBy: ...
    def explain(self) -> str: ...
    def collect(self) -> dict[str, _ColumnArray]:
        """Runs the pipeline and returns a dict of numpy arrays, one per column."""
        ...

class GroupBy:
    """Pending group-by of a `LazyFrame`."""

    def agg(self, aggs: Sequence[Agg]) -> LazyFrame: ...

class Axis:
    """One histogram axis: its bin edges.

    Returned by `Histogram1D.axes` and `Histogram2D.axes`.
    """

    @property
    def edges(self) -> NDArray[np.float64]:
        """Bin edges, shape (nbins + 1,)."""
        ...

    @property
    def centers(self) -> NDArray[np.float64]:
        """Bin centres, shape (nbins,)."""
        ...

    def __len__(self) -> int:
        """Number of regular bins."""
        ...

    def __repr__(self) -> str: ...

class Histogram1D:
    """A 1D histogram with weighted fills, under/overflow and error tracking.

    Hist-style accessors: `axes`, `values()`, `variances()`, `to_numpy()` and the
    `from_numpy()` constructor. `to_numpy()` returns `(values, edges)`.
    """

    def __init__(self, nbins: int, lo: float, hi: float) -> None:
        """Creates `nbins` uniform bins over [lo, hi)."""
        ...

    @staticmethod
    def variable(edges: Sequence[float]) -> Histogram1D:
        """Creates a histogram with explicit, strictly increasing bin edges."""
        ...

    @staticmethod
    def from_numpy(
        values: NDArray[np.float64],
        edges: NDArray[np.float64],
        variances: NDArray[np.float64] | None = None,
    ) -> Histogram1D:
        """Builds a histogram from bin values and edges.

        `values` has length `len(edges) - 1`. `variances` defaults to `values`, the
        Poisson variance of unweighted counts. Flow and fill counts start at zero.
        Raises ValueError for bad edges, length mismatch, or negative or non-finite
        values.
        """
        ...

    @property
    def axes(self) -> tuple[Axis]:
        """The histogram's axis, as a one-element tuple."""
        ...

    def fill(self, x: NDArray[np.float64], weights: NDArray[np.float64] | None = None) -> None:
        """Fills from a 1D float64 array, with optional per-entry weights."""
        ...

    def nbins(self) -> int: ...
    def edges(self) -> NDArray[np.float64]: ...
    def values(self) -> NDArray[np.float64]:
        """Bin contents (sum of weights), shape (nbins,). Same as `contents()`."""
        ...

    def variances(self) -> NDArray[np.float64]:
        """Per-bin variance, sum of w^2, shape (nbins,)."""
        ...

    def to_numpy(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Returns `(values, edges)`."""
        ...

    def bin_centers(self) -> NDArray[np.float64]: ...
    def contents(self) -> NDArray[np.float64]: ...
    def errors(self) -> NDArray[np.float64]:
        """Per-bin statistical errors, sqrt(sum of w^2)."""
        ...

    def bin(self, i: int) -> float: ...
    def bin_error(self, i: int) -> float: ...
    @property
    def underflow(self) -> float: ...
    @property
    def overflow(self) -> float: ...
    @property
    def entries(self) -> int: ...
    @property
    def nan_fills(self) -> int: ...
    def integral(self) -> float: ...
    def mean(self) -> float | None: ...
    def std_dev(self) -> float | None: ...
    def quantile(self, p: float) -> float | None: ...
    def merge(self, other: Histogram1D) -> None:
        """Adds `other` (identical binning required) into this histogram."""
        ...

    def copy(self) -> Histogram1D: ...
    def to_dict(self) -> dict[str, NDArray[np.float64] | float | int]: ...
    def __repr__(self) -> str: ...

class Histogram2D:
    """A 2D histogram with weighted fills. Contents are indexed [ix, iy]."""

    def __init__(
        self, nx: int, xlo: float, xhi: float, ny: int, ylo: float, yhi: float
    ) -> None:
        """Creates `nx * ny` uniform bins over [xlo, xhi) x [ylo, yhi)."""
        ...

    @staticmethod
    def variable(xedges: Sequence[float], yedges: Sequence[float]) -> Histogram2D: ...
    @staticmethod
    def from_numpy(
        values: NDArray[np.float64],
        xedges: NDArray[np.float64],
        yedges: NDArray[np.float64],
        variances: NDArray[np.float64] | None = None,
    ) -> Histogram2D:
        """Builds a histogram from a (nbins_x, nbins_y) array and two edge arrays.

        `variances` defaults to `values`. Raises ValueError if the shape does not
        match the edges, or if a value is non-finite or a variance is negative.
        """
        ...

    @property
    def axes(self) -> tuple[Axis, Axis]:
        """The `(x, y)` axes."""
        ...

    def fill(
        self,
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        weights: NDArray[np.float64] | None = None,
    ) -> None:
        """Fills from two 1D float64 arrays of equal length, with optional weights."""
        ...

    def values(self) -> NDArray[np.float64]:
        """Bin contents, shape (nbins_x, nbins_y). Same as `contents()`."""
        ...

    def variances(self) -> NDArray[np.float64]:
        """Per-bin variance, sum of w^2, shape (nbins_x, nbins_y)."""
        ...

    def to_numpy(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Returns `(values, x_edges, y_edges)`."""
        ...

    def nbins_x(self) -> int: ...
    def nbins_y(self) -> int: ...
    def bin(self, ix: int, iy: int) -> float: ...
    def bin_error(self, ix: int, iy: int) -> float: ...
    def contents(self) -> NDArray[np.float64]:
        """Bin contents, shape (nbins_x, nbins_y)."""
        ...

    def errors(self) -> NDArray[np.float64]:
        """Bin errors, shape (nbins_x, nbins_y)."""
        ...

    def x_edges(self) -> NDArray[np.float64]: ...
    def y_edges(self) -> NDArray[np.float64]: ...
    def integral(self) -> float: ...
    @property
    def entries(self) -> int: ...
    @property
    def nan_fills(self) -> int: ...
    def project_x(self) -> Histogram1D:
        """Projection onto the x axis (sums over y)."""
        ...

    def project_y(self) -> Histogram1D:
        """Projection onto the y axis (sums over x)."""
        ...

    def merge(self, other: Histogram2D) -> None: ...
    def copy(self) -> Histogram2D: ...
    def to_dict(self) -> dict[str, NDArray[np.float64] | int]: ...
    def __repr__(self) -> str: ...

def validate_config(config: str) -> None:
    """Parses and validates a TOML run configuration. Raises ValueError if invalid."""
    ...

def run_simulation(
    config: str,
    structure: str,
    *,
    format: Literal["gro", "pdb"] = "gro",
    sample_every: int = 1,
    minimize: bool | None = None,
    checkpoint_restart: str | None = None,
    checkpoint_write: str | None = None,
) -> dict[str, NDArray[np.float64] | NDArray[np.uint64]]:
    """Runs a simulation from a TOML config and GRO/PDB text.

    Returns arrays `step` (uint64), `time`, `kinetic`, `potential`, `total`
    and `temperature` (float64), sampled at step 0 and every `sample_every` steps.

    The keyword arguments override the TOML settings for this call:
    `minimize` sets `[minimize] enabled`; `checkpoint_restart` restores the
    state from that file before the run (sampling then starts at the
    checkpoint's step); `checkpoint_write` writes the final state to that file.

    The `[output]` table is honoured as in `tptdyn simulate`: `energy_file` and
    `trajectory_file` are written at their configured intervals (paths are
    relative to the working directory). The returned arrays are sampled
    separately and are not read back from those files.
    """
    ...
