import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Literal
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from .base import PlotBase


class RelPlot(PlotBase):
    def __init__(self, **kwargs):
        super().__init__()

    def _default_xaxis(self) -> None:
        max_val = max(np.abs(self.data[self.x]))
        xlim = round(max_val + (0.05 * max_val), 1)
        xlim = (-xlim, xlim)
        self.logger.info(f"Setting x-axis limits to {xlim}.")
        return self.ax.set_xlim(xlim)


    def label_points(self, ax, data: pd.DataFrame, label: str, signif_filter : bool = True, 
                     size: float = 4, max_label: int = 50, ha='center', va='center', **kwargs) -> None:
        """Labels points on the plot using the provided label column, with optional condition and styling."""

        # Ensure labels are within current axis limits
        ax_ylim = ax.get_ylim()
        ax_xlim = ax.get_xlim()
        data = data[(data[self.y] <= ax_ylim[1]) & (data[self.y] >= ax_ylim[0])]
        data = data[(data[self.x] <= ax_xlim[1]) & (data[self.x] >= ax_xlim[0])]

        # Filter for significance
        if signif_filter:
            data = data[data['Significance'] == True]

        data_length = len(data)
        if data_length > max_label:
            self.logger.warning(f"Too many points to label ({data_length}). Labelling 50 most extreme points...")
            sorted_vals = data[self.x].sort_values()
            bottom_25 = sorted_vals[:25].index
            top_25 = sorted_vals[-25:].index
            label_idx = bottom_25.union(top_25)
            data = data.loc[label_idx]
        texts = [ax.text(x=row[self.x], y=row[self.y], s=row[label], size=size, ha=ha, va=va, **kwargs) for _, row in data.iterrows()]
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax, min_arrow_len=0)



class VolcanoPlot(RelPlot):
    def __init__(self, data: pd.DataFrame, x: str, y: str, hue: Literal['Regulation', 'Significance'] = 'Regulation', label : str = None, signif: float = 0.05,
                 **kwargs):
        """
        - data (DataFrame): The input data for plotting.
        - x (str): The column name for the x-axis values.
        - y (str): The column name for the y-axis values (e.g., p-values).
        - signif (float, optional): The threshold for significance.
        - **kwargs: Additional keyword arguments passed to the parent class and plotting functions.
        """
        super().__init__()
        self.orig_data = data
        self.x = x
        self.y = y
        self.hue = hue
        self.label = label
        self.signif = -np.log10(signif)
        self.delta_box_color = kwargs.get('delta_box_color', 'lightgrey')
        self.delta_text_size = kwargs.get('delta_text_size', 6)

        if self.hue == 'Significance':
            self.hue_order = [False, True]
        elif self.hue == 'Regulation':
            self.hue_order = ['notsig', 'up', 'down']

        self.data = self._prepare_data()

        plot_kws = locals()
        plot_kws['data'] = self.data
        plot_kws['hue_order'] = self.hue_order

        self.ax = self.plot(sns.scatterplot, plot_kws)
        self._default_xaxis()
        self._add_threshold_lines()
        self._add_delta_count_box()
        self.ax.set_ylabel('-Log10 (FDR)')
                                                         
    def _prepare_data(self) -> pd.DataFrame:
        data = self.orig_data.copy()
        data[self.y] = -np.log10(data[self.y])
        data['Significance'] = data[self.y] > self.signif
        data['Regulation'] = 'notsig'
        data.loc[(data['Significance'] == True) & (data[self.x] > 0), 'Regulation'] = 'up'
        data.loc[(data['Significance'] == True) & (data[self.x] < 0), 'Regulation'] = 'down'
        return data

    def _add_threshold_lines(self) -> None:
        self.ax.axhline(y=self.signif, color='black', linestyle='--', linewidth=0.5, dashes=(5, 15))
        self.ax.axvline(x=0, color='black', linewidth=0.5)

    def _add_delta_count_box(self, box_position=(0.725, 0.75), box_width=0.25, box_height=0.2) -> None:
        greater_delta = sum(self.data['Regulation'] == 'up')
        less_delta = sum(self.data['Regulation'] == 'down')
        greater_text = f'{self.x} $>$ 0: {greater_delta}'
        less_text = f'{self.x} $<$ 0: {less_delta}'
        text_y_greater = box_position[1] + box_height * 0.75
        text_y_less = box_position[1] + box_height * 0.25
        text_x = box_position[0] + 0.025
        self.ax.add_patch(Rectangle(box_position, box_width, box_height,
                                    facecolor=self.delta_box_color, alpha=0.75, transform=self.ax.transAxes))
        self.ax.text(x=text_x, y=text_y_greater, s=greater_text, ha='left', va='center',
                     transform=self.ax.transAxes, fontsize=self.delta_text_size)
        self.ax.text(x=text_x, y=text_y_less, s=less_text, ha='left', va='center',
                     transform=self.ax.transAxes, fontsize=self.delta_text_size)