[![Pylint](https://github.com/OrhunKok/proteoboost/actions/workflows/pylint.yml/badge.svg)](https://github.com/OrhunKok/proteoboost/actions/workflows/pylint.yml)

# ProteoBoost

ProteoBoost is a Python package designed for processing, analyzing, and visualizing proteomics data. It provides utilities for data loading, processing, transformation, and visualization, making it easier to work with complex proteomics datasets.

## Features

- **Data Loading**: Load data from various formats, including Parquet and Excel.
- **Data Processing**: Clean and structure raw proteomics data for downstream analysis.
- **Matrix Transformation**: Build matrices for statistical and machine learning workflows.
- **Visualization**: Generate publication-quality plots for data exploration and presentation.
- **Custom Logging**: Integrated logging for better traceability and debugging.

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd proteoboost
   ```

2. Install the package and its dependencies:
   ```bash
   pip install .
   ```

3. Alternatively, install the package in editable mode for development:
   ```bash
   pip install -e .
   ```

## Dependencies

ProteoBoost requires the following Python packages:

- `pandas` (for data manipulation)
- `numpy` (for numerical operations)
- `pyarrow` (for handling Parquet files)
- `matplotlib` (for plotting)
- `scienceplots` (for plot theming)
- `seaborn` (for advanced visualizations)

## Contributing

Contributions are welcome! If you'd like to contribute, please fork the repository, make your changes, and submit a pull request. Ensure that your code passes all tests and follows the project's coding standards.

## License

This project is licensed under the MIT License.