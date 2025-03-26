import matplotlib.pyplot as plt
import inspect
import contextlib
import scienceplots  # type: ignore
from ..utils.logging import MetaLogging


class PlotBase(metaclass=MetaLogging):
    """
    PlotBase is a base class for creating and managing plots with customizable themes.
    Attributes:
        theme (str): The theme to be applied to the plots. Defaults to 'science'.
        ax (matplotlib.axes.Axes): The current Axes object for the plot.
    Methods:
        __init__(theme: str = 'science'):
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
    """

    def __init__(self, theme: str = "science"):
        self.theme = theme

    def filter_kwargs(self, func, kwargs):
        """
        Filters a dictionary of keyword arguments to include only those accepted by a given function.
        Args:
            func (Callable): The function whose accepted arguments will be used for filtering.
            kwargs (dict): A dictionary of keyword arguments to filter.
        Returns:
            dict: A dictionary containing only the keyword arguments that are accepted by the function.
        """
        func_args = inspect.signature(func).parameters
        ignore_keys = {"self", "__class__"}

        return {
            k: v for k, v in kwargs.items() if k in func_args and k not in ignore_keys
        }

    def __getattr__(self, name):
        if self.ax and hasattr(self.ax, name):
            return getattr(self.ax, name)
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

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
        """
        Generates a plot using the specified plotting function and keyword arguments.

        This method applies the current plot theme, filters the provided keyword arguments
        to ensure compatibility with the plotting function, and then calls the function
        to create the plot.

        Args:
            plot_func (callable): The plotting function to be used for generating the plot.
            plot_kws (dict): A dictionary of keyword arguments to be passed to the plotting
                function. The dictionary should include a key 'kwargs' containing additional
                arguments to be unpacked.

        Returns:
            matplotlib.axes.Axes: The Axes object of the generated plot.
        """
        with self.plot_theme():
            plot_kws = self.filter_kwargs(plot_func, plot_kws)
            kwargs = plot_kws.pop("kwargs")
            self.ax = plot_func(**plot_kws, **kwargs)
        return self.ax

    def save(self, *args, **kwargs):
        """
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
        """
        if self.ax:
            self.ax.figure.savefig(*args, **kwargs)
            plt.close(self.ax.figure)
        else:
            print("Save path or Figure not available. Figure not saved.")

    def show(self):
        """
        Display the plot if it has been generated.

        This method checks if the `ax` attribute is not `None`. If a plot
        exists, it will render the plot using `matplotlib.pyplot.show()`.
        Otherwise, it will print a message indicating that no plot has
        been generated yet.
        """
        if self.ax is not None:
            plt.show()
        else:
            print("No plot has been generated yet.")
