"""Relational plotting utilities (scatter, pairwise, joint plots).

This module implements higher-level relational plots used for exploratory
data analysis in proteolyzer. Plots accept data frames or arrays and
return Matplotlib axes or Seaborn objects.

Example
    >>> from proteolyzer.plots.relational import scatter_plot
    >>> ax = scatter_plot(df, x="intensity", y="ratio")

"""

import numpy as np
import seaborn as sns
import pandas as pd
from typing import Literal, List
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from .base import PlotBase


class RelPlot(PlotBase):
    def __init__(self, **kwargs):
        super().__init__()

    def _symmetric_xaxis(self) -> None:
        max_val = max(np.abs(self.data[self.x]))
        xlim = round(max_val + (0.05 * max_val), 1)
        xlim = (-xlim, xlim)
        self.logger.info(f"Setting x-axis limits to {xlim}.")
        return self.ax.set_xlim(xlim)

    def label_points(
        self,
        signif_filter: bool = True,
        size: float = 4,
        max_label: int = 50,
        ha="center",
        va="center",
        **kwargs,
    ) -> None:
        # Ensure labels are within current axis limits
        ax_ylim = self.ax.get_ylim()
        ax_xlim = self.ax.get_xlim()

        data = self.data.copy()
        data = data[(data[self.y] <= ax_ylim[1]) & (data[self.y] >= ax_ylim[0])]
        data = data[(data[self.x] <= ax_xlim[1]) & (data[self.x] >= ax_xlim[0])]

        # Filter for significance
        if signif_filter:
            data = data[data["Significance"] == True]

        data_length = len(data)
        if data_length > max_label:
            self.logger.warning(
                f"Too many points to label ({data_length}). Labelling {max_label} most extreme points..."
            )
            # Ensure max_label is even
            label_count = max_label // 2

            sorted_vals = data[self.x].sort_values()
            bottom = sorted_vals[:label_count].index
            top = sorted_vals[-label_count:].index
            label_idx = bottom.union(top)
            data = data.loc[label_idx]
        texts = [
            self.ax.text(
                x=row[self.x],
                y=row[self.y],
                s=row[self.label],
                size=size,
                ha=ha,
                va=va,
                **kwargs,
            )
            for _, row in data.iterrows()
        ]
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="-", color="k", lw=0.5),
            ax=self.ax,
            min_arrow_len=0,
            clip_on=True,
        )

    def add_data_point_count(
        self,
        x_pos: float = 0.7625,
        y_pos: float = 0.825,
        fontsize: int = 6,
        ha: str = "left",
        va: str = "top",
        color="black",
        **kwargs,
    ):
        text = f"$\\textit{{n}}$  = {len(self.ax.collections[0].get_offsets())}"
        self.ax.text(
            x_pos,
            y_pos,
            text,
            ha=ha,
            va=va,
            color=color,
            fontsize=fontsize,
            transform=self.ax.transAxes,
            **kwargs,
        )


class VolcanoPlot(RelPlot):
    def __init__(
        self,
        data: pd.DataFrame,
        x: str,
        y: str,
        hue: Literal["Regulation", "Significance", str] = "Regulation",
        hue_order: List = None,
        label: str = None,
        signif: float = 0.05,
        **kwargs,
    ):
        super().__init__()
        self.orig_data = data
        self.x = x
        self.y = y
        self.hue = hue
        self.hue_order = hue_order
        self.label = label
        self.signif = -np.log10(signif)
        self.delta_text_size = kwargs.get("delta_text_size", 6)

        if self.hue == "Significance":
            hue_order = [False, True]
        elif self.hue == "Regulation":
            hue_order = ["notsig", "up", "down"]

        self.data = self._prepare_data()
        plot_kws = locals()
        plot_kws["data"] = self.data
        plot_kws["hue_order"] = hue_order

        self.ax = self.plot(sns.scatterplot, plot_kws)
        # self._symmetric_xaxis()
        self._add_threshold_lines()
        self._add_delta_count_box()
        self.ax.set_ylim(bottom=0)

    def _prepare_data(self) -> pd.DataFrame:
        data = self.orig_data.copy()
        if data[self.y].max() <= 1:
            self.logger.info(f"-Log10 normalizing {self.y}")
            data[self.y] = -np.log10(data[self.y])
            data["Significance"] = data[self.y] > self.signif
        elif data[self.y].max() > 1:
            self.logger.info(
                f"{self.y} already -Log10 transformed. Skipping transform..."
            )
            data["Significance"] = data[self.y] > self.signif

        data["Regulation"] = "notsig"
        data.loc[(data["Significance"] == True) & (data[self.x] > 0), "Regulation"] = (
            "up"
        )
        data.loc[(data["Significance"] == True) & (data[self.x] < 0), "Regulation"] = (
            "down"
        )
        return data

    def _add_threshold_lines(self) -> None:
        self.ax.axhline(
            y=self.signif, color="black", linestyle="--", linewidth=0.5, dashes=(1, 5)
        )
        self.ax.axvline(x=0, color="black", linewidth=0.5)

    def _add_delta_count_box(
        self, box_position=(0.675, 0.85), box_width=0.3, box_height=0.125
    ) -> None:
        greater_delta = sum(self.data["Regulation"] == "up")
        less_delta = sum(self.data["Regulation"] == "down")
        count = len(self.data)

        max_digits = max(len(str(greater_delta)), len(str(less_delta)), len(str(count)))

        greater_text_part = f"{self.x}$>$0: "
        greater_num_part = f"{greater_delta:>{max_digits}}"

        less_text_part = f"{self.x}$<$0: "
        less_num_part = f"{less_delta:>{max_digits}}"

        # Calculate vertical positions with spacing
        text_y_greater = box_position[1] + box_height * 0.9
        text_y_less = box_position[1] + box_height * 0.35

        text_x_text = box_position[0] + ((box_position[1] - box_position[0]) * 0.1)
        text_x_num = box_position[0] + box_width * 0.9

        self.ax.add_patch(
            Rectangle(
                box_position,
                box_width,
                box_height,
                facecolor="none",
                edgecolor="black",
                linewidth=0.5,
                alpha=0.75,
                transform=self.ax.transAxes,
            )
        )

        self.ax.text(
            x=text_x_text,
            y=text_y_greater,
            s=greater_text_part,
            ha="left",
            va="top",
            transform=self.ax.transAxes,
            fontsize=self.delta_text_size,
        )
        self.ax.text(
            x=text_x_num,
            y=text_y_greater,
            s=greater_num_part,
            ha="right",
            va="top",
            transform=self.ax.transAxes,
            fontsize=self.delta_text_size,
        )

        self.ax.text(
            x=text_x_text,
            y=text_y_less,
            s=less_text_part,
            ha="left",
            va="top",
            transform=self.ax.transAxes,
            fontsize=self.delta_text_size,
        )
        self.ax.text(
            x=text_x_num,
            y=text_y_less,
            s=less_num_part,
            ha="right",
            va="top",
            transform=self.ax.transAxes,
            fontsize=self.delta_text_size,
        )
