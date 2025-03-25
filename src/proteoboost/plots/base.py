import matplotlib.pyplot as plt
import inspect
import contextlib
import os
import scienceplots # type: ignore
from ..utils.logging import MetaLogging

class PlotBase(metaclass = MetaLogging):
    def __init__(self, theme : str = 'science'):
        self.theme = theme

    def filter_kwargs(self, func, kwargs):
        """Filters kwargs based on the function's accepted arguments."""
        func_args = inspect.signature(func).parameters
        ignore_keys = {'self', '__class__'}

        return {k: v for k, v in kwargs.items() if k in func_args and k not in ignore_keys}

    def __getattr__(self, name):
        if self.ax and hasattr(self.ax, name):
            return getattr(self.ax, name)
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    @contextlib.contextmanager
    def plot_theme(self):
        """A context manager that sets a specified plot theme."""
        try:
            with plt.style.context(self.theme):
                yield
        except ValueError as e:
            print(f"Error applying theme '{self.theme}': {e}")
            yield 
        except ModuleNotFoundError:
            print("scienceplots is not installed.")
            yield

    def plot(self, plot_func, plot_kws):
        with self.plot_theme():
            plot_kws = self.filter_kwargs(plot_func, plot_kws)
            kwargs = plot_kws.pop('kwargs')
            self.ax = plot_func(**plot_kws, **kwargs)
        return self.ax
    
    def save(self, *args, **kwargs):
        """Saves the current figure"""
        if self.ax:
            self.ax.figure.savefig(*args, **kwargs)
            plt.close(self.ax.figure)
        else:
            print("Save path or Figure not available. Figure not saved.")

    def show(self):
        """Shows the plot using plt.show()."""
        if self.ax is not None:
            plt.show()
        else:
            print("No plot has been generated yet.")