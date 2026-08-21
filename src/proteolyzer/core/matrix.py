"""Matrix transformation utilities.

Provide performant, well-documented routines to perform common matrix
operations required by the analysis pipeline, including normalization,
scaling and basic imputations.

Example
    >>> import numpy as np
    >>> builder = MatrixBuilder(processed_data)  # doctest: +SKIP
    >>> builder.matrix_generation(
    ...     "Ms1.Area", index=["Precursor.Id"], columns=["Run"]
    ... ).normalize_matrix(within_groups=["Run"], agg_func=np.nansum).matrix
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .logging import Logged
from .models import Report


@dataclass(frozen=True)
class Missingness:
    """How much of a matrix is absent, and from where.

    Gaps and zeros are counted apart because they mean different things: an NA
    is a measurement that was never made (missing at random), a zero is one
    that came back empty (missing not at random), and imputation has to treat
    them differently.
    """

    #: Percentage of cells that are NA.
    mar: float
    #: Percentage of cells that are exactly zero.
    mnar: float
    #: Fraction of each column that is absent, counting both gaps and zeros.
    per_column: pd.Series

    def sparse_columns(self, threshold: float = 0.75) -> list:
        """Columns more than `threshold` absent, as candidates for dropping."""
        if not 0 <= threshold <= 1:
            raise ValueError(f"threshold must be between 0 and 1, got {threshold}")
        return list(self.per_column[self.per_column > threshold].index)


class MatrixBuilder(Logged):
    """Pivots long-form processed data into a quantitative matrix."""

    __slots__ = ("data", "matrix")

    #: Set by :meth:`matrix_generation`.
    matrix: pd.DataFrame | None

    def __init__(self, data: Report | pd.DataFrame):
        """Takes a :class:`Report` or a plain frame."""
        self.data = data.frame if isinstance(data, Report) else data
        self.matrix = None

    def missingness(
        self, matrix: pd.DataFrame | None = None, warning_threshold: float = 0.75
    ) -> Missingness:
        """Measure how much of the matrix is absent, logging the headline numbers.

        Defaults to the generated matrix. Returns the numbers as well as
        logging them, so a caller can act on them rather than read them.
        """
        if matrix is None:
            matrix = self._require_matrix()

        cells = matrix.shape[0] * matrix.shape[1]
        if not cells:
            self.logger.warning("Matrix is empty; nothing to measure.")
            return Missingness(np.nan, np.nan, pd.Series(dtype="float64"))

        absent = matrix.isna() | (matrix == 0)
        report = Missingness(
            mar=round((matrix.isna().to_numpy().sum() / cells) * 100, 2),
            mnar=round(((matrix == 0).to_numpy().sum() / cells) * 100, 2),
            per_column=absent.mean(),
        )
        self.logger.info(f"Data has {report.mar}% Missing At Random (MAR)")
        self.logger.info(f"Data has {report.mnar}% Missing Not At Random (MNAR)")

        try:
            sparse = report.sparse_columns(warning_threshold)
        except ValueError as exc:
            self.logger.error(f"{exc}. Ignoring further missingness checks...")
            return report

        if sparse:
            self.logger.warning(
                f"Missingness over {round(warning_threshold * 100)}% for "
                f"{sparse}, recommend dropping these."
            )
        return report

    def matrix_generation(
        self, values: str, index: list[str], columns: list[str]
    ) -> MatrixBuilder:
        try:
            self.matrix = self.data.pivot(index=index, columns=columns, values=values)
        except ValueError as exc:
            # pivot rejects duplicate index/column pairs but does not say how
            # many there are. Counting them is another pass over the frame, so
            # it happens only once the pivot has actually failed.
            duplicates = self.data.duplicated(subset=[*index, *columns], keep=False)
            if not duplicates.any():
                raise
            raise ValueError(
                f"{duplicates.sum()} Duplicate combinations of '{index}' and "
                f"'{columns}' found. Matrix cannot be created."
            ) from exc

        self.missingness(self.matrix)

        return self

    def _require_matrix(self) -> pd.DataFrame:
        if self.matrix is None:
            raise ValueError("No matrix yet. Call matrix_generation() first.")
        return self.matrix

    def normalize_matrix(
        self, within_groups: list[str], agg_func: Callable, replace_zeros: bool = True
    ) -> MatrixBuilder:
        """Divide each value by an aggregate of its row within each column group.

        Parameters
        ----------
        within_groups
            Column-index level(s) defining the groups aggregated over.
        agg_func
            Row-wise aggregation, called as ``agg_func(block, axis=1)`` (e.g.
            ``np.nansum``, ``np.nanmedian``).
        replace_zeros
            Treat zeros as missing before normalizing.
        """
        matrix = self._require_matrix()
        matrix_norm = matrix.replace(0, np.nan) if replace_zeros else matrix

        # Nullable integer/float dtypes hold pd.NA, which numpy cannot reduce
        # over, so normalize on a plain float64 copy. The copy also has to be
        # writeable: to_numpy can hand back a read-only view of the frame.
        norm_data = matrix_norm.to_numpy(dtype="float64", na_value=np.nan, copy=True)

        cols_df = matrix_norm.columns.to_frame(index=False)
        grouped = cols_df.groupby(within_groups, sort=False).indices

        for group in grouped.values():
            sub_data = norm_data[:, group]
            row_sums = agg_func(sub_data, axis=1)
            row_sums = np.where(row_sums == 0, np.nan, row_sums)
            norm_data[:, group] = sub_data / row_sums[:, None]

        self.matrix = pd.DataFrame(
            norm_data, index=matrix.index, columns=matrix.columns
        )
        return self
