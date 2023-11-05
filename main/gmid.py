# -----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
# -----------------------------------------------------------------------------#
from functools import partial
import numpy as np
from scipy.interpolate import interpn
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from matplotlib.widgets import Cursor

# don't warn user about bad divisions
np.seterr(divide="ignore", invalid="ignore")

# plot settings {{{
FIG_SIZE = (8, 4)
LINE_WIDTH = 1.5
GRID_COLOR = "0.9"

PLOT_SETTINGS = {
    "FIG_SIZE": FIG_SIZE,
    "LINE_WIDTH": LINE_WIDTH,
    "GRID_COLOR": GRID_COLOR,
}


def set_plot_settings(var_name, new_value):
    global FIG_SIZE, LINE_WIDTH, GRID_COLOR
    if var_name in PLOT_SETTINGS:
        globals()[var_name] = new_value


# }}}

# matplotlib interraction {{{
dots = []
annotations = []


def on_canvas_click(event, fig, ax):
    if fig.canvas.toolbar.mode != "":
        return

    x = event.xdata
    y = event.ydata

    if x is None or y is None:
        return

    print(f"x={x}, y={y}")

    (dot,) = ax.plot(x, y, "ro")
    dots.append(dot)

    formatter = EngFormatter()
    x_eng = formatter(x)
    y_eng = formatter(y)
    annotation = ax.annotate(
        f"({x_eng}, {y_eng})",
        (x, y),
        textcoords="offset points",
        xytext=(0, 10),
        ha="center",
    )
    annotations.append(annotation)

    fig.canvas.draw()


def clear_annotations_and_dots(fig):
    for dot in dots:
        dot.remove()
    for annotation in annotations:
        annotation.remove()
    fig.canvas.draw()
    dots.clear()
    annotations.clear()


# }}}


def load_lookup_table(path: str):
    return np.load(path, allow_pickle=True).tolist()


class GMID:
    # init method {{{
    def __init__(
        self, lookup_table: dict, *, mos: str, vsb=None, vgs=None, vds=None, slice_independent=None
    ):
        # extract from table
        self.parameters = lookup_table["parameter_names"]
        self.w = lookup_table["w"]
        self.l = lookup_table["l"]
        self.lengths = self.l  # just an alias

        self.mos = mos
        self.lookup_table = lookup_table[mos]
        self.lookup_table["l"] = self.lengths
        self.lookup_table["w"] = self.w
        self.slice_independent = slice_independent

        self.__parse_source_list(locals())
        self.vsb, self.vgs, self.vds = self.voltage_sources.values()

        # define commonly-used expressions to avoid typing them every time
        self.__common_expressions()
        self.__common_plot_methods()

        # extract parameters by varying the independent source
        self.__extract_parameters_by_independent_source()

    # }}}

    # private function: parse argument list to determine the independent source{{{
    def __parse_source_list(self, args: dict):
        sources = {"vsb": None, "vgs": None, "vds": None}
        independent_variable = list(sources.keys())
        fixed_variables = []

        for k, v in args.items():
            if k in sources.keys() and v is not None:
                sources[k] = v
                fixed_variables.append(k)
                independent_variable.remove(k)

        self.independent_variable = independent_variable[0]
        self.independent_source = self.lookup_table[self.independent_variable]

        sources[self.independent_variable] = self.independent_source
        self.voltage_sources = sources
        self.fixed_variables = fixed_variables

    # }}}

    # private function: parse argument list of lookup function {{{
    def __parse_arg_list(self, arg: dict):
        variables = {"length": None, "vsb": None, "vgs": None, "vds": None}
        primary = arg["primary"]
        primary_variable = ""
        secondary_variable = ""
        for k, v in arg.items():
            if k in variables.keys():
                if isinstance(v, tuple):
                    if k == primary:
                        primary_variable = k
                    else:
                        secondary_variable = k
                    variables[k] = np.arange(*v)
                else:
                    variables[k] = v
        if not primary_variable and secondary_variable:
            primary_variable, secondary_variable = secondary_variable, ""
        return variables, primary_variable, secondary_variable

    # }}}

    # private function: calculate parameter from expression {{{
    def __calculate_from_expression(
        self,
        expression: dict,
        table: dict,
        filter_by_rows: np.ndarray = np.array([]),
    ):
        var_list = []
        for v in expression["variables"]:
            var = table[v]
            if filter_by_rows is not None and filter_by_rows.size > 0 and var.size > 1:
                var_list.append((np.take(var, filter_by_rows, 0)))
            else:
                var_list.append(var)
        if "function" in expression:
            values = expression["function"](*var_list)
        else:
            values = var_list[0]
        try:
            return values, expression["label"]
        except KeyError:
            return values, ""

    # }}}

    # private function: extract parameters by varying the independent source {{{
    def __extract_parameters_by_independent_source(self):
        title_label = []
        lookup_table = self.lookup_table
        slice_idx = {
            "length": slice(None),
            "vsb": slice(None),
            "vgs": slice(None),
            "vds": slice(None),
        }

        # find closest match to the dependent sources in the lookup table
        for s in self.fixed_variables:
            value = self.voltage_sources[s]
            slice_idx[s] = (np.abs(lookup_table[s] - value)).argmin()
            title_label.append(f"V_{{ \\mathrm{{ { (s[1:]).upper() } }} }}")
            title_label.append(value)

        # TODO-SLICING: I'm not that happy with the way I implemented this!
        if self.slice_independent:
            var = self.independent_variable
            start = self.slice_independent[0]
            end = self.slice_independent[1]
            slice_idx[var] = (lookup_table[var] >= start) & (lookup_table[var] <= end)

        # extract parameters based on the slice
        extracted_table = {}
        for p in self.parameters:
            extracted_table[p] = lookup_table[p][tuple(slice_idx.values())]

        # TODO-SLICING
        if not self.slice_independent:
            extracted_table[self.independent_variable] = np.tile(
                lookup_table[self.independent_variable], (len(self.lengths), 1)
            )
        else:
            extracted_table[self.independent_variable] = np.tile(
                lookup_table[self.independent_variable][slice_idx[var]], (len(self.lengths), 1)
            )

        extracted_table["w"] = np.array(self.w)
        extracted_table["l"] = self.l
        extracted_table[
            "title"
        ] = f"{lookup_table['model_name']}, " + "$%s=%.2f$, $%s=%.2f$" % tuple(title_label)

        legend_formatter = EngFormatter(unit="m")
        extracted_table["label"] = [legend_formatter.format_eng(sw) for sw in self.l]

        if self.independent_variable == "vds":
            self.independent_source_expression = self.vds_expression
        elif self.independent_variable == "vgs":
            self.independent_source_expression = self.vgs_expression
        elif self.independent_variable == "vsb":
            self.independent_source_expression = self.vsb_expression

        self.extracted_table = extracted_table
        # }}}

    # private function: plot settings {{{
    def __plot_settings(
        self,
        fig,
        ax,
        y: np.ndarray,
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        x_label: str = "",
        y_label: str = "",
        legend: list = [],
        title: str = "",
        save_fig: str = "",
    ):
        if x_scale:
            ax.set_xscale(x_scale)

        if y_scale:
            ax.set_yscale(y_scale)
        else:
            # TODO: this might not always give expected result
            if np.max(y) / np.min(y) > 1000:
                ax.set_yscale("log")

        formatter0 = EngFormatter(unit="")
        if y_eng_format:
            ax.yaxis.set_major_formatter(formatter0)
        if x_eng_format:
            ax.xaxis.set_major_formatter(formatter0)

        if x_limit:
            ax.set_xlim(*x_limit)

        if y_limit:
            ax.set_ylim(*y_limit)

        if legend:
            ax.legend(legend, loc="center left", bbox_to_anchor=(1, 0.5))

        ax.set_title(title)
        ax.grid(True, which="both", ls="--", color=GRID_COLOR)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        if save_fig:
            fig.savefig(save_fig)

    # }}}

    # private function: look for a value using a calculated expression {{{
    def __lookup_by(
        self,
        length: float,
        independent_expression: dict,
        val: float,
        dependent_expression: dict,
    ):
        values, _ = self.__calculate_from_expression(independent_expression, self.extracted_table)
        g = values[(np.abs(self.extracted_table["l"] - length)).argmin()]
        if dependent_expression:
            values, _ = self.__calculate_from_expression(dependent_expression, self.extracted_table)
        else:
            values, _ = self.__calculate_from_expression(dependent_expression, self.extracted_table)
        point = np.array([length, val])
        return interpn((self.lookup_table["l"], g), values, point)

    # }}}

    # quick plot {{{
    def quick_plot(
        self,
        x: np.ndarray | list | tuple,
        y: np.ndarray | list | tuple,
        *,
        x_label: str = "",
        y_label: str = "",
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        legend: list = [],
        save_fig: str = "",
    ):
        """
        Make quick plots. As a reminder, when `x` and `y` are of size m x n, pass
        them to this function as x.T and y.T
        """
        fig, ax = plt.subplots(1, 1, figsize=FIG_SIZE, tight_layout=True)
        fig.canvas.mpl_connect("button_press_event", lambda event: on_canvas_click(event, fig, ax))
        fig.canvas.mpl_connect(
            "key_press_event",
            lambda event: clear_annotations_and_dots(fig) if event.key == "d" else None,
        )

        if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
            ax.plot(x, y, lw=LINE_WIDTH, picker=True)

        if isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            for x_, y_ in zip(x, y):
                ax.plot(x_, y_, lw=LINE_WIDTH, picker=True)
            y = y[0]

        self.__plot_settings(
            fig,
            ax,
            y,
            x_label=x_label,
            y_label=y_label,
            x_limit=x_limit,
            y_limit=y_limit,
            y_scale=y_scale,
            x_scale=x_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format,
            legend=legend,
            save_fig=save_fig,
        )

    # }}}

    # function: plot by expression {{{
    def plot_by_expression(
        self,
        *,
        x_axis: dict,
        y_axis: dict,
        lengths: tuple = (),
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        save_fig: str = "",
        return_result: bool = False,
    ):
        extracted_table = self.extracted_table

        # TODO: this does not always work due to floating-point errors
        lengths_idx = None
        if lengths:
            lengths_idx = np.where(np.in1d(extracted_table["l"], np.array(lengths)))[0]

        if lengths_idx is not None and lengths_idx.size > 0:
            legend = [extracted_table["label"][i] for i in lengths_idx]
        else:
            legend = extracted_table["label"]

        x, x_label = self.__calculate_from_expression(x_axis, extracted_table, lengths_idx)
        y, y_label = self.__calculate_from_expression(y_axis, extracted_table, lengths_idx)

        fig, ax = plt.subplots(1, 1, figsize=(8, 4), tight_layout=True)
        fig.canvas.mpl_connect("button_press_event", lambda event: on_canvas_click(event, fig, ax))
        fig.canvas.mpl_connect(
            "key_press_event",
            lambda event: clear_annotations_and_dots(fig) if event.key == "d" else None,
        )

        ax.plot(x.T, y.T, lw=LINE_WIDTH, picker=True)

        self.__plot_settings(
            fig,
            ax,
            y,
            x_limit,
            y_limit,
            x_scale,
            y_scale,
            x_eng_format,
            y_eng_format,
            x_label,
            y_label,
            legend,
            extracted_table["title"],
            save_fig,
        )

        plt.show()

        if return_result:
            return x.T, y.T

    # }}}

    # function: plot by sweep {{{
    def plot_by_sweep(
        self,
        *,
        length: float | tuple,
        vsb: float | tuple,
        vgs: float | tuple,
        vds: float | tuple,
        primary: str = "",
        x_axis_expression: dict | np.ndarray,
        y_axis_expression: dict | np.ndarray,
        x_label: str = "",
        y_label: str = "",
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        save_fig: str = "",
        return_result: bool = False,
    ):
        _, primary_variable, secondary_variable = self.__parse_arg_list(locals())

        if primary_variable not in x_axis_expression["variables"]:
            x = self.lookup(
                length=length,
                vsb=vsb,
                vgs=vgs,
                vds=vds,
                expression=x_axis_expression,
                primary=primary,
            )
        else:
            x = np.arange(*(locals()[primary_variable]))

        if primary_variable not in y_axis_expression["variables"]:
            y = self.lookup(
                length=length,
                vsb=vsb,
                vgs=vgs,
                vds=vds,
                expression=y_axis_expression,
                primary=primary,
            )
        else:
            y = np.arange(*(locals()[primary_variable]))

        if x.ndim < y.ndim:
            x = np.tile(x, y.shape[0]).reshape(y.shape[0], -1)
        elif y.ndim < x.ndim:
            y = np.tile(y, x.shape[0]).reshape(x.shape[0], -1)

        if x.ndim > 1:
            x = x.T
            y = y.T

        if secondary_variable:
            legend_formatter = EngFormatter(unit="m")
            legend = [
                legend_formatter.format_eng(sw) for sw in np.arange(*locals()[secondary_variable])
            ]
        else:
            legend = []

        if not x_label:
            try:
                x_label = x_axis_expression["label"]
            except KeyError:
                x_label = ""
            except TypeError:
                x_label = x_axis_expression

        if not y_label:
            try:
                y_label = y_axis_expression["label"]
            except KeyError:
                y_label = ""
            except TypeError:
                y_label = y_axis_expression

        self.quick_plot(
            x,
            y,
            x_label=x_label,
            y_label=y_label,
            x_limit=x_limit,
            y_limit=y_limit,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format,
            legend=legend,
            save_fig=save_fig,
        )
        if return_result:
            return (x, y)

    # }}}

    # function: lookup a parameter from the table usign values used in defining the table {{{
    def lookup(
        self,
        *,
        length: float | np.ndarray,
        vsb: float | np.ndarray,
        vgs: float | np.ndarray,
        vds: float | np.ndarray,
        expression: str | dict,
        primary: str = "",
    ):
        """
        Return the interpolated value calculated by `expression` for  given
        `length`, `vsb`, `vgs`, and `vds` variables.
        Any two of these variables can take a range in the form (start, stop, step).
        When two variables take a range, the arguement `primary` must be used to
        indicate which of the two is the primary variable.
        --------
        Example:
        --------
            nmos.lookup(
                length = 180e-9,
                vsb = 0,
                vgs = (0.1, 1.8, 0.01),
                vds = (0.4, 1.4+0.5, 0.5),
                expression = nmos.gmid_expression,
                primary = "vgs"
            )
        """
        l = self.lengths
        s = self.lookup_table["vsb"]
        g = self.lookup_table["vgs"]
        d = self.lookup_table["vds"]
        points = (l, s, d, g)

        values, _ = self.__calculate_from_expression(expression, self.lookup_table)

        # parse argument list
        variables, primary_variable, secondary_variable = self.__parse_arg_list(locals())

        if not (primary_variable or secondary_variable):
            point = np.array(list(variables.values()))
            return interpn(points, values, point)
        else:
            n = variables[primary_variable].size
            if secondary_variable:
                m = variables[secondary_variable].size
                secondary_range = variables[secondary_variable]
            else:
                m = 1
                secondary_range = np.array([1])  # put any value to make the for loop start
            result = np.zeros(shape=(m, n))
            for i, s in enumerate(secondary_range):
                stack_points = []
                for k, v in variables.items():
                    if k == primary_variable:
                        stack_points.append(v)
                    elif k == secondary_variable:  # this condition will never be meet when m=1
                        stack_points.append(np.repeat(s, n))
                    else:
                        stack_points.append(np.repeat(v, n))
                point = np.column_stack(tuple(stack_points))
                result[i] = interpn(points, values, point)
            if secondary_variable:
                return result
            else:
                return result[0]

    # }}}

    # function: lookup a parameter by gmid from the extracted table {{{
    def lookup_by_gmid(self, length: float | tuple, gmid: float, expression: dict):
        """
        Return the interpolated value calculated by `expression` for  given
        `length` and `gmid` variables.
        `length` can take a range in the form (start, stop, step)
        --------
        Example:
        --------
            nmos.lookup_by_gmid(
                length=450e-9,
                gmid=15,
                expression=nmos.gain_expression
            )
        """
        return self.lookup_by_expression(
            length=length,
            look_by_expression=self.gmid_expression,
            look_by_value=gmid,
            look_for_expression=expression,
        )

    # }}}

    # function: lookup a parameter by an expression from the extracted table {{{
    def lookup_by_expression(
        self,
        length: float | tuple,
        look_by_expression: dict,
        look_by_value,
        look_for_expression: dict,
    ):
        if isinstance(length, tuple):
            length = np.arange(*length)
            result = np.empty(shape=(length.size,))
            for i, v in enumerate(length):
                result[i] = self.__lookup_by(
                    v, look_by_expression, look_by_value, look_for_expression
                )
            return result
        elif isinstance(length, float):
            return self.__lookup_by(length, look_by_expression, look_by_value, look_for_expression)

    # }}}

    # private function: define commonly-used expressions {{{
    def __common_expressions(self):
        # create attributes for parameters from the lookup table
        LABEL_TABLE = {
            "vsb": ["V_{\\mathrm{SB}}", "V"],
            "vgs": ["V_{\\mathrm{GS}}", "V"],
            "vds": ["V_{\\mathrm{DS}}", "V"],
            "id": ["I_{D}", "A"],
            "vth": ["V_{\\mathrm{TH}}", "V"],
            "vdsat": ["V_{\\mathrm{DS_{\\mathrm{SAT}}}}", "V"],
            "gm": ["g_{m}", "S"],
            "gmbs": ["g_{\\mathrm{mbs}}", "S"],
            "gds": ["g_{\\mathrm{ds}}", "S"],
            "cgg": ["c_{\\mathrm{gg}}", "F"],
            "cgs": ["c_{\\mathrm{gs}}", "F"],
            "cbg": ["c_{\\mathrm{bg}}", "F"],
            "cgd": ["c_{\\mathrm{gd}}", "F"],
            "cdd": ["c_{\\mathrm{dd}}", "F"],
        }
        for parameter, (label, unit) in LABEL_TABLE.items():
            if parameter in self.parameters or parameter in ["vsb", "vgs", "vds"]:
                setattr(
                    self,
                    f"{parameter}_expression",
                    {"variables": [parameter], "label": f"${label}\ ({unit})$"},
                )

        self.gmid_expression = {
            "variables": ["gm", "id"],
            "function": lambda x, y: x / y,
            "label": "$g_m/I_D (S/A)$",
        }
        self.vstar_expression = {
            "variables": ["gm", "id"],
            "function": lambda x, y: (2 * y) / x,
            "label": "$V^{\\star} (V)$",
        }
        self.gain_expression = {
            "variables": ["gm", "gds"],
            "function": lambda x, y: x / y,
            "label": "$g_{m}/g_{\\mathrm{ds}}$",
        }
        self.current_density_expression = {
            "variables": ["id", "w"],
            "function": lambda x, y: x / y,
            "label": "$I_{D}/W (A/m)$",
        }
        self.transist_frequency_expression = {
            "variables": ["gm", "cgg"],
            "function": lambda x, y: x / (2 * np.pi * y),
            "label": "$f_{T} (Hz)$",
        }
        self.early_voltage_expression = {
            "variables": ["id", "gds"],
            "function": lambda x, y: x / y,
            "label": "$V_{A} (V)$",
        }

    # }}}

    # create commonly-used plot functions {{{
    def __common_plot_methods(self):
        PLOT_METHODS = {
            "current_density_plot": [self.gmid_expression, self.current_density_expression],
            "gain_plot": [self.gmid_expression, self.gain_expression],
            "transit_frequency_plot": [self.gmid_expression, self.transist_frequency_expression],
            "early_voltage_plot": [self.gmid_expression, self.early_voltage_expression],
        }

        # This is not ideal since I need to keep track of the signature of `plot_by_expression`.
        # I can use `partial` from `functools` here, but it doesn't hide the bound variables, which I don't like.
        def create_plot_method(self, x_axis, y_axis):
            def plot_method(
                self,
                lengths: tuple = (),
                x_limit: tuple = (),
                y_limit: tuple = (),
                x_scale: str = "",
                y_scale: str = "",
                x_eng_format: bool = False,
                y_eng_format: bool = False,
                save_fig: str = "",
                return_result: bool = False,
            ):
                return self.plot_by_expression(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    lengths=lengths,
                    x_scale=x_scale,
                    y_scale=y_scale,
                    x_limit=x_limit,
                    y_limit=y_limit,
                    x_eng_format=x_eng_format,
                    y_eng_format=y_eng_format,
                    save_fig=save_fig,
                    return_result=return_result,
                )

            return plot_method

        for method_name, (x, y) in PLOT_METHODS.items():
            setattr(GMID, method_name, create_plot_method(self, x_axis=x, y_axis=y))


# }}}

# vim:fdm=marker
