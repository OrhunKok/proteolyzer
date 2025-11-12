import numpy as np

def cv(data, min_datapoints : int = 3):
    data = np.asarray(data)
    if data.size < min_datapoints:
        return np.nan
    mean = np.mean(data)
    std_dev = np.std(data, ddof=1)
    if mean == 0:
        return np.nan
    return std_dev / mean