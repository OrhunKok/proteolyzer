"""Matrix transformation utilities.

Provide performant, well-documented routines to perform common matrix
operations required by the analysis pipeline, including normalization,
scaling and basic imputations.

Example
    >>> from proteolyzer.transformers.matrix import normalize_rows
    >>> norm = normalize_rows(matrix, method="zscore")
"""

import pandas as pd
import numpy as np
from typing import Callable
from proteolyzer.utils.logging import MetaLogging
from proteolyzer.utils.models import ProcessedData


class MatrixBuilder(metaclass=MetaLogging):
    __slots__ = ("processed_data", "data", "matrix", "logger")

    def __init__(self, processed_data: ProcessedData):
        self.data = processed_data

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
                "warning_threshold must be between 0 and 1. Ignoring further missingness checks..."
            )
        else:
            replaced_mnar = matrix.replace(0, np.nan)
            recommend_drop = (
                replaced_mnar.apply(lambda x: x.isna().sum() / len(x))
                > warning_threshold
            )
            recommend_drop = list(recommend_drop[recommend_drop == True].index)
            if len(recommend_drop) > 0:
                self.logger.warning(
                    f"Missingness over {round(warning_threshold * 100)}% for {recommend_drop}, recommend dropping these."
                )

    def matrix_generation(
        self, values: str, index: list[str], columns: list[str]
    ) -> pd.DataFrame:
        duplicate_count = self.data.duplicated(
            subset=[*index, *columns], keep=False
        ).sum()

        if duplicate_count > 0:
            raise ValueError(
                f"{duplicate_count} Duplicate combinations of '{index}' and '{columns}' found. Matrix cannot be created."
            )
        else:
            self.matrix = self.data.pivot(index=index, columns=columns, values=values)

        self.missingness_check(self.matrix)

        return self

    def normalize_matrix(
        self, within_groups: list[str], agg_func: Callable, replace_zeros: bool = True
    ) -> pd.DataFrame:
        matrix_norm = (
            self.matrix.replace(0, np.nan) if replace_zeros else self.matrix.copy()
        )

        norm_data = matrix_norm.values.copy()

        cols_df = matrix_norm.columns.to_frame(index=False)
        grouped = cols_df.groupby(within_groups, sort=False).indices

        for _, group in grouped.items():
            sub_data = norm_data[:, group]
            row_sums = agg_func(sub_data, axis=1)
            row_sums[row_sums == 0] = np.nan
            norm_data[:, group] = sub_data / row_sums[:, None]

        self.matrix = pd.DataFrame(
            norm_data, index=self.matrix.index, columns=self.matrix.columns
        )
        return self
