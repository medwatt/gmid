# imports <<<
from typing import List, Optional, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import EngFormatter
# >>>


class Plotter:
    def __init__(
        self,
        fig_size: Tuple[int, int] = (8, 4),
        line_width: float = 1.5,
        grid_color: str = "0.9",
    ) -> None:
        self.fig_size = fig_size
        self.line_width = line_width
        self.grid_color = grid_color

    # single axis plot <<<
    def create_figure(
        self,
        title: str = "",
        x_label: str = "",
        y_label: str = "",
        x_lim: Optional[Tuple[float, float]] = None,
        y_lim: Optional[Tuple[float, float]] = None,
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
    ) -> Tuple[plt.Figure, plt.Axes]:

        fig, ax = plt.subplots(figsize=self.fig_size, tight_layout=True)
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(True, which="both", ls="--", color=self.grid_color)

        if x_lim is not None:
            ax.set_xlim(*x_lim)

        if y_lim is not None:
            ax.set_ylim(*y_lim)

        if x_scale:
            ax.set_xscale(x_scale)

        if y_scale:
            ax.set_yscale(y_scale)

        if x_eng_format:
            ax.xaxis.set_major_formatter(EngFormatter(unit=""))

        if y_eng_format:
            ax.yaxis.set_major_formatter(EngFormatter(unit=""))

        return fig, ax

    # >>>

    # two axis plot <<<
    def create_figure_with_twin(
        self,
        title: str = "",
        x_label: str = "",
        y_label: str = "",
        y2_label: str = "",
        x_lim: Optional[Tuple[float, float]] = None,
        y_lim: Optional[Tuple[float, float]] = None,
        y2_lim: Optional[Tuple[float, float]] = None,
        x_scale: str = "",
        y_scale: str = "",
        y2_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        y2_eng_format: bool = False,
    ) -> Tuple[plt.Figure, plt.Axes, plt.Axes]:

        fig, ax = plt.subplots(figsize=self.fig_size, tight_layout=True)
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax2 = ax.twinx()
        ax2.set_ylabel(y2_label)
        ax.grid(True, which="both", ls="--", color=self.grid_color)

        if x_lim is not None:
            ax.set_xlim(*x_lim)

        if y_lim is not None:
            ax.set_ylim(*y_lim)

        if x_scale:
            ax.set_xscale(x_scale)

        if y_scale:
            ax.set_yscale(y_scale)

        if x_eng_format:
            ax.xaxis.set_major_formatter(EngFormatter(unit=""))

        if y_eng_format:
            ax.yaxis.set_major_formatter(EngFormatter(unit=""))

        if y2_lim is not None:
            ax2.set_ylim(*y2_lim)

        if y2_scale:
            ax2.set_yscale(y2_scale)

        if y2_eng_format:
            ax2.yaxis.set_major_formatter(EngFormatter(unit=""))

        return fig, ax, ax2
    # >>>

    # plot data <<<
    def plot_data(
        self,
        ax: plt.Axes,
        x: Union[np.ndarray, List[np.ndarray]],
        y: Union[np.ndarray, List[np.ndarray]],
        line_style: Optional[str] = "solid",
        legend: Optional[Union[np.ndarray, List[Union[float, str]]]] = None,
        legend_title: Optional[str] = None,
        legend_placement: Optional[str] = None,
        bbox_to_anchor: Optional[Tuple[float, float]] = None,
        save_fig: str = "",
        end_plotting: bool = True,
    ) -> None:

        # Plot data when x and y are numpy arrays.
        if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
            if x.ndim > 1 and y.ndim > 1:
                ax.plot(x.T, y.T, lw=self.line_width, ls=line_style, picker=True)
            else:
                ax.plot(x, y, lw=self.line_width, ls=line_style, picker=True)

        # Plot data when x and y are lists or tuples.
        elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            for x_item, y_item in zip(x, y):
                if isinstance(x_item, np.ndarray) and isinstance(y_item, np.ndarray):
                    if x_item.ndim > 1 and y_item.ndim > 1:
                        ax.plot(
                            x_item.T,
                            y_item.T,
                            lw=self.line_width,
                            ls=line_style,
                            picker=True,
                        )
                    else:
                        ax.plot(
                            x_item,
                            y_item,
                            lw=self.line_width,
                            ls=line_style,
                            picker=True,
                        )
                else:
                    ax.plot(
                        np.asarray(x_item),
                        np.asarray(y_item),
                        lw=self.line_width,
                        ls=line_style,
                        picker=True,
                    )
        else:
            # Fallback: convert inputs to arrays.
            ax.plot(
                np.asarray(x),
                np.asarray(y),
                lw=self.line_width,
                ls=line_style,
                picker=True,
            )

        if end_plotting:
            if legend is not None:
                # Determine if legend items are numeric.
                if isinstance(legend[0], (int, float, np.number)):
                    formatter = EngFormatter(unit="")
                    formatted_legend = [formatter(val) for val in legend]
                else:
                    formatted_legend = legend

                if legend_placement == "bottom":
                    anchor = (0.5, -0.25) if bbox_to_anchor is None else bbox_to_anchor
                    leg = ax.legend(
                        formatted_legend,
                        loc="upper center",
                        bbox_to_anchor=anchor,
                        ncol=len(formatted_legend),
                        title=legend_title,
                    )

                elif legend_placement == "top":
                    anchor = (0.5, 1.05) if bbox_to_anchor is None else bbox_to_anchor
                    leg = ax.legend(
                        formatted_legend,
                        loc="lower center",
                        bbox_to_anchor=anchor,
                        ncol=len(formatted_legend),
                        title=legend_title,
                    )
                else:
                    anchor = (1, 0.5) if bbox_to_anchor is None else bbox_to_anchor
                    leg = ax.legend(
                        formatted_legend,
                        loc="center left",
                        bbox_to_anchor=anchor,
                        title=legend_title,
                    )

                leg.get_title().set_fontsize("large")

            if save_fig:
                ax.figure.savefig(save_fig, bbox_inches="tight")
            plt.show()
    # >>>
