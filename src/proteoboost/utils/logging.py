import logging
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(asctime)s - %(funcName)s - %(name)s - %(message)s",
)

class MetaLogging(type):
    def __new__(cls, name, bases, dct):

        def _memory_check(self, data : pd.DataFrame) -> None:
            """
            Calculates and logs the memory usage of the provided data.

            Args:
                data: The pandas DataFrame to check memory usage for.
            """
            try:
                summed_bytes = sum(data.memory_usage(deep=True))
                total_mib = round(summed_bytes / (1024 ** 2), 1)
                self.logger.info(f"Total memory usage of data: {total_mib} MiB")
            except ImportError:
                self.logger.error("Pandas library is required for memory check.")
            except AttributeError:
                self.logger.error("The provided data does not have a memory_usage method. Make sure the data is a Pandas DataFrame")
            except Exception as e:
                self.logger.error(f"An unexpected error occurred during memory check: {e}")

        dct['_memory_check'] = _memory_check
        log_class = super().__new__(cls, name, bases, dct)
        log_class.logger = logging.getLogger(name)

        return log_class
    
