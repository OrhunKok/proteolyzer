---
sidebar_label: base
title: proteolyzer.plots.base
---

## plt

## inspect

## contextlib

## scienceplots

## MetaLogging

## PlotBase Objects

```python
class PlotBase(metaclass=MetaLogging)
```

PlotBase is a base class for creating and managing plots with customizable themes.

**Attributes**:

- `theme` _str_ - The theme to be applied to the plots. Defaults to &#x27;science&#x27;.
- `ax` _matplotlib.axes.Axes_ - The current Axes object for the plot.

**Methods**:

- `__init__(theme` - str = &#x27;science&#x27;):
  Initializes the PlotBase object with a specified theme.
  filter_kwargs(func, kwargs):
  Filters the provided keyword arguments to include only those accepted by the given function.
  __getattr__(name):
  Dynamically retrieves attributes from the current Axes object (`ax`) if available.
  plot_theme():
  A context manager that applies the specified plot theme during plotting.
  plot(plot_func, plot_kws):
  Generates a plot using the specified plotting function and keyword arguments.
  save(*args, **kwargs):
  Saves the current figure to a file and closes it.
  show():
  Displays the current plot using `plt.show()`.

#### \_\_init\_\_

```python
def __init__(theme: str = "science")
```

#### filter\_kwargs

```python
def filter_kwargs(func, kwargs)
```

Filters a dictionary of keyword arguments to include only those accepted by a given function.

**Arguments**:

- `func` _Callable_ - The function whose accepted arguments will be used for filtering.
- `kwargs` _dict_ - A dictionary of keyword arguments to filter.

**Returns**:

- `dict` - A dictionary containing only the keyword arguments that are accepted by the function.

#### \_\_getattr\_\_

```python
def __getattr__(name)
```

#### plot\_theme

```python
@contextlib.contextmanager
def plot_theme()
```

A context manager that sets a specified plot theme.

#### plot

```python
def plot(plot_func, plot_kws)
```

Generates a plot using the specified plotting function and keyword arguments.

This method applies the current plot theme, filters the provided keyword arguments
to ensure compatibility with the plotting function, and then calls the function
to create the plot.

**Arguments**:

- `plot_func` _callable_ - The plotting function to be used for generating the plot.
- `plot_kws` _dict_ - A dictionary of keyword arguments to be passed to the plotting
  function. The dictionary should include a key &#x27;kwargs&#x27; containing additional
  arguments to be unpacked.
  

**Returns**:

- `matplotlib.axes.Axes` - The Axes object of the generated plot.

#### save

```python
def save(*args, **kwargs)
```

Saves the current figure to a file.

Parameters
----------
*args : tuple
    Positional arguments to be passed to `matplotlib.figure.Figure.savefig`.
**kwargs : dict
    Keyword arguments to be passed to `matplotlib.figure.Figure.savefig`.

Behavior
--------
- If the `ax` attribute is set, the associated figure is saved using the
  provided arguments and then closed to free up resources.
- If the `ax` attribute is not set, a message is printed indicating that
  the figure could not be saved.

#### show

```python
def show()
```

Display the plot if it has been generated.

This method checks if the `ax` attribute is not `None`. If a plot
exists, it will render the plot using `matplotlib.pyplot.show()`.
Otherwise, it will print a message indicating that no plot has
been generated yet.

