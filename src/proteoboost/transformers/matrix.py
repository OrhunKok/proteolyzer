import pandas as pd
import numpy as np
from proteoboost.utils.logging import MetaLogging
from proteoboost.utils.models import ProcessedData


class MatrixBuilder(metaclass=MetaLogging):
    __slots__ = ("processed_data", "data", "logger")

    def __init__(self, processed_data: ProcessedData):
        self.data = processed_data.data

    def _missingness_check(
        self, matrix: pd.DataFrame, warning_threshold: float = 0.75
    ) -> None:
        total_vals = matrix.shape[0] * matrix.shape[1]

        MAR = round((matrix.isna().sum().sum() / total_vals) * 100)
        MNAR = round((1 - (matrix[matrix == 0].isna().sum().sum() / total_vals)) * 100)
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
            recommend_drop = list(recommend_drop[recommend_drop is True].index)
            if len(recommend_drop) > 0:
                self.logger.warning(
                    f"Missingness over {round(warning_threshold * 100)}% for {recommend_drop}, recommend dropping these."
                )

    def matrix_generation(self, values: str, index: str, columns: str) -> pd.DataFrame:
        duplicate_count = self.data.duplicated(
            subset=[index, columns], keep=False
        ).sum()

        if duplicate_count > 0:
            raise ValueError(
                f"{duplicate_count} Duplicate combinations of '{index}' and '{columns}' found. Matrix cannot be created."
            )
        else:
            matrix = self.data.pivot(index=index, columns=columns, values=values)

        self._missingness_check(matrix)

        return matrix
