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

import numpy as np
import pandas as pd

from .logging import Logged
from .models import Report


class MatrixBuilder(Logged):
    """Pivots long-form processed data into a quantitative matrix."""

    __slots__ = ("data", "matrix")

    #: Set by :meth:`matrix_generation`.
    matrix: pd.DataFrame | None

    def __init__(self, data: Report | pd.DataFrame):
        """Takes a :class:`Report` or a plain frame."""
        self.data = data.frame if isinstance(data, Report) else data
        self.matrix = None

    def missingness_check(
        self, matrix: pd.DataFrame, warning_threshold: float = 0.75
    ) -> None:
        total_vals = matrix.shape[0] * matrix.shape[1]

        MAR = round((matrix.isna().sum().sum() / total_vals) * 100, 2)
        MNAR = round((matrix[matrix == 0].notna().sum().sum() / total_vals) * 100, 2)
        self.logger.info(f"Data has {MAR}% Missing At Random (MAR)")
        self.logger.info(f"Data has {MNAR}% Missing Not At Random (MNAR)")

        if not 0 <= warning_threshold <= 1:
            self.logger.error(
                "warning_threshold must be between 0 and 1. Ignoring further "
                "missingness checks..."
            )
            return

        replaced_mnar = matrix.replace(0, np.nan)
        recommend_drop = (
            replaced_mnar.apply(lambda x: x.isna().sum() / len(x)) > warning_threshold
        )
        recommend_drop = list(recommend_drop[recommend_drop == True].index)
        if len(recommend_drop) > 0:
            self.logger.warning(
                f"Missingness over {round(warning_threshold * 100)}% for "
                f"{recommend_drop}, recommend dropping these."
            )

    def matrix_generation(
        self, values: str, index: list[str], columns: list[str]
    ) -> MatrixBuilder:
        duplicate_count = self.data.duplicated(
            subset=[*index, *columns], keep=False
        ).sum()

        if duplicate_count > 0:
            raise ValueError(
                f"{duplicate_count} Duplicate combinations of '{index}' and "
                f"'{columns}' found. Matrix cannot be created."
            )

        self.matrix = self.data.pivot(index=index, columns=columns, values=values)
        self.missingness_check(self.matrix)

        return self

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
        if self.matrix is None:
            raise ValueError("No matrix to normalize. Call matrix_generation() first.")

        matrix_norm = self.matrix.replace(0, np.nan) if replace_zeros else self.matrix

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
            norm_data, index=self.matrix.index, columns=self.matrix.columns
        )
        return self
