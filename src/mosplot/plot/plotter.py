# imports <<<
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
# >>>


class Plotter:
    def __init__(self, fig_size=(8, 4), line_width=1.5, grid_color="0.9"):
        self.fig_size = fig_size
        self.line_width = line_width
        self.grid_color = grid_color
        self.dots = []
        self.annotations = []

    def create_figure(
        self,
        title="",
        x_label="",
        y_label="",
        x_lim=None,
        y_lim=None,
        x_scale="",
        y_scale="",
        x_eng_format=False,
        y_eng_format=False,
    ):
        fig, ax = plt.subplots(figsize=self.fig_size, tight_layout=True)
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", color=self.grid_color)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if x_lim:
            ax.set_xlim(*x_lim)
        if y_lim:
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

    def plot_data(self, ax, x, y, legend=None, save_fig=""):
        if isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and x.ndim == y.ndim:
            ax.plot(x.T, y.T, lw=self.line_width, picker=True)
        elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            for x_, y_ in zip(x, y):
                if x_.ndim == 1 and x_.shape[0] != y_.shape[0]:
                    ax.plot(x_, y_.T, lw=self.line_width, picker=True)
                else:
                    ax.plot(x_, y_, lw=self.line_width, picker=True)
        elif x.ndim == 1:
            if x.shape[0] != y.shape[0]:
                ax.plot(x, y.T, lw=self.line_width, picker=True)
            else:
                ax.plot(x, y, lw=self.line_width, picker=True)
        elif y.ndim == 1:
            if y.shape[0] != x.shape[0]:
                ax.plot(x.T, y, lw=self.line_width, picker=True)
            else:
                ax.plot(x, y, lw=self.line_width, picker=True)
        if legend:
            ax.legend(legend, loc="center left", bbox_to_anchor=(1, 0.5))
        if save_fig:
            ax.figure.savefig(save_fig, bbox_inches="tight")

        plt.show()
