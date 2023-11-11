# -----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
# -----------------------------------------------------------------------------#

import numpy as np
from scipy.interpolate import interpn, griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from matplotlib.widgets import Cursor

# don't warn user about bad divisions
np.seterr(divide="ignore", invalid="ignore")

################################################################################
#                               Helper Functions                               #
################################################################################
# load the generated lookup table
def load_lookup_table(path: str):
    return np.load(path, allow_pickle=True).tolist()

# tile one array so that it has the same shape as the other
def tile_arrays(A, B):
    if A.ndim == 1 and B.ndim == 2:
        if A.shape[0] == B.shape[0]:
            return np.tile(A, (B.shape[1], 1)).T, B
        elif A.shape[0] == B.shape[1]:
            return np.tile(A, (B.shape[0], 1)), B
    elif B.ndim == 1 and A.ndim == 2:
        if B.shape[0] == A.shape[0]:
            return A, np.tile(B, (A.shape[1], 1)).T
        elif B.shape[0] == A.shape[1]:
            return A, np.tile(B, (A.shape[0], 1))
    return A, B

################################################################################
#                    Override Plot Settings During Runtime                     #
################################################################################
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

################################################################################
#                         Matplotlib Plot Interraction                         #
################################################################################
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

################################################################################
#                                     GMID                                     #
################################################################################
class LoadMosfet:
    def __init__( self, *, lookup_table, mos, lengths=None, vsb=None, vgs=None, vds=None, primary=None):
        """
        Initialize a mosfet object.
        Two of `lengths, vsb, vgs, vds` must be fixed at any time.

        Args:
            lookup_table (dict): dictionary of mosfet parameters
            mos (str): type of mosfet: "nmos" or "pmos"
            lengths (float, list, ndarray): length(s) of the mosfet
            vsb (float, tuple): source-body voltage: tuple of the form (start, stop, step)
            vgs (float, tuple): gate-source voltage: tuple of the form (start, stop, step)
            vds (float, tuple): drain-source voltage: tuple of the form (start, stop, step)
            primary (str): name of the primary sweep source

        Example:
            nmos = LoadMosfet(lookup_table=lookup_table, mos="nmos", vsb=0.0, vds=0.5, vgs=(0.3, 1))
        """
        # extract data from lookup table
        self.mos = mos
        self.lookup_table = lookup_table
        self.width = lookup_table[mos]["width"]
        self.lengths = lookup_table[mos]["lengths"]
        self.parameters = lookup_table[mos]["parameter_names"]

        # extract a 2d table of the parameters
        self.secondary_variable_idx, self.filtered_variables, self.extracted_table = \
        self.extract_2d_table(lookup_table=self.lookup_table[self.mos], primary=primary, lengths=lengths, vsb=vsb, vgs=vgs, vds=vds)
        self.lengths, self.vsb, self.vgs, self.vds = self.filtered_variables

        # define commonly-used expressions to avoid typing them every time
        self.__common_expressions()
        self.__common_plot_methods()

    def extract_2d_table(self, *, lookup_table, parameters=None, lengths=None, vsb=None, vgs=None, vds=None, primary=None):
        """
        Filter the lookup table based

        Args:
            lookup_table (dict): dictionary of parameters of one of the mosfets
            lengths (float, list, ndarray): length(s) of the mosfet
            vsb (float, tuple): source-body voltage: tuple of the form (start, stop, step)
            vgs (float, tuple): gate-source voltage: tuple of the form (start, stop, step)
            vds (float, tuple): drain-source voltage: tuple of the form (start, stop, step)
            primary (str): name of the primary sweep source

        Returns:
            secondary_idx: index of the secondary sweep variable in `lengths, vsb, vgs, vds`
            filter_values: filtered values of `lengths, vsb, vgs, vds`
            extracted_table: filtered values of parameters
        """
        # check that at least two parameters are provided
        params = [lengths is not None, vsb is not None, vgs is not None, vds is not None]
        if sum(params) < 2:
            raise ValueError("Please provide at least two parameters.")

        def get_indices(var, target):
            data = lookup_table[var]
            if isinstance(target, tuple): # when `vsb, vgs, vds` is a range
                start, end = target[:2]
                indices = np.where((data >= start) & (data <= end))[0]
                if len(target) == 3:  # if it contains a step
                    step = int(target[2] / (data[1] - data[0]))
                    indices = indices[::step]
            # lengths must be handled separately since they are provided as list of values
            elif var == "lengths" and isinstance(target, (list, np.ndarray)):
                # filter by lengths
                mask = np.isin(self.lookup_table["lengths"], np.array(target))
                indices = np.nonzero(mask)[0]
                indices = np.array(indices, dtype=int)
            else:
                # by eliminating all float variables, we will be left with one variable
                # that is not a float, and that will be the secondary variable
                variables[var] = True
                index = (np.abs(data - target)).argmin()
                return np.array([index]), data[index]

            return indices, data[indices]

        secondary_idx = None
        variables = {"lengths": False, "vsb": False, "vgs": False, "vds": False}
        if primary:
            variables[primary] = True

        indices_and_values = {
            "lengths": get_indices("lengths", lengths) if lengths is not None else (slice(None), lookup_table["lengths"]),
            "vsb": get_indices("vsb", vsb) if vsb is not None else (slice(None), lookup_table["vsb"]),
            "vgs": get_indices("vgs", vgs) if vgs is not None else (slice(None), lookup_table["vgs"]),
            "vds": get_indices("vds", vds) if vds is not None else (slice(None), lookup_table["vds"]),
        }

        slice_indices = []
        filter_values = []
        for idx, item in enumerate(variables.keys()):
            slice_indices.append(indices_and_values[item][0])
            filter_values.append(indices_and_values[item][1])
        slice_indices = tuple(slice_indices)

        def slice_me(a, slices):
            x = a[slices[0], :, :, :]
            x = x[:, slices[1], :, :]
            x = x[:, :, slices[2], :]
            x = x[:, :, :, slices[3]]
            return x

        # extract the data based on the indices
        extracted_table = {}
        if not parameters:
            parameters = lookup_table["parameter_names"]

        for p in parameters:
            if p in lookup_table:
                x = np.squeeze(slice_me(lookup_table[p], slice_indices))
                if x.ndim > 1 and x.shape[0] > x.shape[1]:
                    extracted_table[p] = x.T
                else:
                    extracted_table[p] = x

        one_key = next(iter(extracted_table))
        extracted_table["width"] = np.array(self.width)
        extracted_table["lengths"], _ = tile_arrays(filter_values[0], extracted_table[one_key])
        extracted_table["vsb"], _ = tile_arrays(filter_values[1], extracted_table[one_key])
        extracted_table["vgs"], _ = tile_arrays(filter_values[2], extracted_table[one_key])
        extracted_table["vds"], _ = tile_arrays(filter_values[3], extracted_table[one_key])

        if primary and secondary_idx:
            secondary_idx = list(variables.values()).index(False)

        return secondary_idx, filter_values, extracted_table

    ################################################################################
    #                               Plotting Methods                               #
    ################################################################################
    def __generate_plot_labels(self):
        variables_labels = ["lengths", "vsb", "vgs", "vds"]
        model_name = self.lookup_table[self.mos]["model_name"]
        title_label = []

        for i, v in enumerate(self.filtered_variables):
            if not isinstance(v, np.ndarray):
                label = variables_labels[i]
                title_label.append(f"V_{{ \\mathrm{{ { (label[1:]).upper() } }} }}")
                title_label.append(v)

        self.plot_labels = {}
        self.plot_labels["title"] = f"{model_name}, " + "$%s=%.2f$, $%s=%.2f$" % tuple(title_label)

        legend_formatter = EngFormatter(unit="m")
        self.plot_labels["lengths"] = [legend_formatter.format_eng(sw) for sw in self.lengths]

    def __plot_settings(
        self,
        y: np.ndarray,
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        x_label: str = "",
        y_label: str = "",
        title: str = "",
        save_fig: str = "",
    ):
        fig, ax = plt.subplots(1, 1, figsize=FIG_SIZE, tight_layout=True)
        fig.canvas.mpl_connect(
            "button_press_event",
            lambda event: on_canvas_click(event, fig, ax),
        )
        fig.canvas.mpl_connect(
            "key_press_event",
            lambda event: clear_annotations_and_dots(fig) if event.key == "d" else None,
        )

        ax.set_title(title)
        ax.grid(True, which="both", ls="--", color=GRID_COLOR)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        if x_limit:
            ax.set_xlim(*x_limit)

        if y_limit:
            ax.set_ylim(*y_limit)

        if x_scale:
            ax.set_xscale(x_scale)

        if y_scale:
            ax.set_yscale(y_scale)
        else:
            # TODO: this might not always give expected result
            if np.max(y) / np.min(y) > 1000:
                ax.set_yscale("log")

        # set engineering format if specified
        if y_eng_format:
            ax.yaxis.set_major_formatter(EngFormatter(unit=""))
        if x_eng_format:
            ax.xaxis.set_major_formatter(EngFormatter(unit=""))

        return fig, ax

    def __plot(self, x, y, fig, ax, legend, save_fig):
        if isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and x.ndim == y.ndim:
            ax.plot(x.T, y.T, lw=LINE_WIDTH, picker=True)

        elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            for x_, y_ in zip(x, y):
                if x_.ndim == 1 and x_.shape[0] != y_.shape[0]:
                    ax.plot(x_, y_.T, lw=LINE_WIDTH, picker=True)
                else:
                    ax.plot(x_, y_, lw=LINE_WIDTH, picker=True)

        elif x.ndim == 1:
            if x.shape[0] != y.shape[0]:
                ax.plot(x, y.T, lw=LINE_WIDTH, picker=True)
            else:
                ax.plot(x, y, lw=LINE_WIDTH, picker=True)

        elif y.ndim == 1:
            if y.shape[0] != x.shape[0]:
                ax.plot(x.T, y, lw=LINE_WIDTH, picker=True)
            else:
                ax.plot(x, y, lw=LINE_WIDTH, picker=True)


        if legend:
            ax.legend(legend, loc="center left", bbox_to_anchor=(1, 0.5))

        if save_fig:
            ax.figure.savefig(save_fig, bbox_inches="tight")

    def plot_by_expression(
        self,
        *,
        x_expression: dict,
        y_expression: dict,
        lengths: tuple = (),
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        title: str = None,
        save_fig: str = "",
        return_result: bool = False,
    ):
        extracted_table = self.extracted_table

        # plot labels
        self.__generate_plot_labels()
        if title is not None:
            self.plot_labels["title"] = title

        # filter by lengths
        mask = np.isin(self.lengths, np.array(lengths))
        indices = np.nonzero(mask)[0]
        length_indices = np.array(indices, dtype=int)

        if length_indices.size > 0:
            legend = [self.plot_labels["lengths"][i] for i in length_indices]
        else:
            legend = self.plot_labels["lengths"]

        x, x_label = self.__calculate_from_expression(x_expression, extracted_table, length_indices)
        y, y_label = self.__calculate_from_expression(y_expression, extracted_table, length_indices)

        fig, ax = self.__plot_settings(
            y,
            x_limit,
            y_limit,
            x_scale,
            y_scale,
            x_eng_format,
            y_eng_format,
            x_label,
            y_label,
            self.plot_labels["title"],
            save_fig,
        )

        self.__plot(x, y, fig, ax, legend, save_fig)

        if return_result:
            return x, y

    def plot_by_sweep(
        self,
        *,
        lengths,
        vsb,
        vgs,
        vds,
        primary,
        x_expression,
        y_expression,
        title: str = "",
        x_limit: tuple = (),
        y_limit: tuple = (),
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        save_fig: str = "",
        return_result: bool = False,
    ):
        secondary_variable_idx, filtered_variables, extracted_table = self.extract_2d_table(
            lookup_table=self.lookup_table[self.mos], lengths=lengths, vsb=vsb, vgs=vgs, vds=vds, primary=primary
        )

        x, x_label = self.__calculate_from_expression(x_expression, extracted_table)
        y, y_label = self.__calculate_from_expression(y_expression, extracted_table)

        fig, ax = self.__plot_settings(
            y,
            x_limit,
            y_limit,
            x_scale,
            y_scale,
            x_eng_format,
            y_eng_format,
            x_label,
            y_label,
            title,
            save_fig,
        )

        if secondary_variable_idx:
            legend_formatter = EngFormatter(unit="")
            legend = [legend_formatter.format_eng(sw) for sw in filtered_variables[secondary_variable_idx]]
        else:
            legend = None

        self.__plot(x, y, fig, ax, legend, save_fig)

        if return_result:
            return x, y

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
        title: str = None,
        save_fig: str = "",
    ):
        """
        Make quick plots. As a reminder, when `x` and `y` are of size m x n, pass
        them to this function as x.T and y.T
        """
        fig, ax = self.__plot_settings(
            y,
            x_limit,
            y_limit,
            x_scale,
            y_scale,
            x_eng_format,
            y_eng_format,
            x_label,
            y_label,
            title,
            save_fig,
        )

        self.__plot(x, y, fig, ax, legend, save_fig)

    # }}}

    ################################################################################
    #                             Expression Handling                              #
    ################################################################################
    def __calculate_from_expression(
        self,
        expression: dict,
        table: dict,
        filter_by_rows: np.ndarray = np.array([]),
    ):
        if isinstance(expression, dict):
            var_list = []
            for v in expression["variables"]:
                var = table[v]
                if filter_by_rows.size > 0 and var.ndim > 1:
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
        else:
            return expression, None

    def __common_expressions(self):
        # create attributes for parameters from the lookup table
        LABEL_TABLE = {
            "lengths": ["\\mathrm{Length}", "m"],
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
            if parameter in self.parameters or parameter in ["lengths", "vsb", "vgs", "vds"]:
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
            "variables": ["id", "width"],
            "function": lambda x, y: x / y,
            "label": "$I_{D}/W (A/m)$",
        }
        self.transist_frequency_expression = {
            "variables": ["gm", "cgg"],
            "function": lambda x, y: x / (2 * np.pi * y),
            "label": "$f_{T} (\\mathrm{Hz})$",
        }
        self.early_voltage_expression = {
            "variables": ["id", "gds"],
            "function": lambda x, y: x / y,
            "label": "$V_{A} (V)$",
        }
        self.rds_expression = {
            "variables": ["gds"],
            "function": lambda x: 1 / x,
            "label": "$r_{\\mathrm{ds}} (\\Omega)$",
        }

    def __common_plot_methods(self):
        PLOT_METHODS = {
            "current_density_plot": [self.gmid_expression, self.current_density_expression],
            "gain_plot": [self.gmid_expression, self.gain_expression],
            "transit_frequency_plot": [self.gmid_expression, self.transist_frequency_expression],
            "early_voltage_plot": [self.gmid_expression, self.early_voltage_expression],
        }

        # This is not ideal since I need to keep track of the signature of `plot_by_expression`.
        # I can use `partial` from `functools` here, but it doesn't hide the bound variables, which I don't like.
        def create_plot_method(self, x_expression, y_expression):
            def plot_method(
                self,
                lengths: tuple = (),
                x_limit: tuple = (),
                y_limit: tuple = (),
                x_scale: str = "",
                y_scale: str = "",
                x_eng_format: bool = False,
                y_eng_format: bool = False,
                title: str = None,
                save_fig: str = "",
                return_result: bool = False,
            ):
                return self.plot_by_expression(
                    x_expression=x_expression,
                    y_expression=y_expression,
                    lengths=lengths,
                    x_scale=x_scale,
                    y_scale=y_scale,
                    x_limit=x_limit,
                    y_limit=y_limit,
                    x_eng_format=x_eng_format,
                    y_eng_format=y_eng_format,
                    title=title,
                    save_fig=save_fig,
                    return_result=return_result,
                )

            return plot_method

        for method_name, (x, y) in PLOT_METHODS.items():
            setattr(LoadMosfet, method_name, create_plot_method(self, x_expression=x, y_expression=y))

    ################################################################################
    #                                Lookup Methods                                #
    ################################################################################
    def interpolate(self, *, x_expression, x_value, y_expression, y_value, z_expression):
        """
        Given (1) a value from x_expression,
              (2) a value from y_expression,
        find value of z_expression using interpolation.

        Args:
            x_expression (dict): expression of how to calculate the points on the x-axis
            x_value (float, dict): value(s) inside the domain of x_expression
            y_expression (dict): expression of how to calculate the points on the y-axis
            y_value (float, dict): value(s) inside the domain of y_expression
            z_expression (dict): expression of how to calculate the value you're looking for

        Returns:
            value of expression you're looking for

        Example:
            x = nmos.interpolate(
                x_expression=nmos.vgs_expression,
                x_value=0.65,
                y_expression=nmos.gmid_expression,
                y_value=15,
                z_expression=nmos.lengths_expression,
            )
        """
        x_array, _ = self.__calculate_from_expression(x_expression, self.extracted_table)
        y_aray, _ = self.__calculate_from_expression(y_expression, self.extracted_table)
        z_array, _ = self.__calculate_from_expression(z_expression, self.extracted_table)

        points = np.column_stack((x_array.ravel(), y_aray.ravel()))

        if isinstance(x_value, (tuple, np.ndarray)) and isinstance(y_value, (int, float)):
            if isinstance(x_value, tuple):
                x = np.arange(*x_value)
            else:
                x = x_value
            evaluate_at = np.column_stack((x, np.full(x.shape, y_value)))
        elif isinstance(y_value, (tuple, np.ndarray)) and isinstance(x_value, (int, float)):
            if isinstance(y_value, tuple):
                y = np.arange(*y_value)
            else:
                y = y_value
            evaluate_at = np.column_stack((np.full(y.shape, x_value), y))
        elif isinstance(x_value, tuple) and isinstance(y_value, tuple):
            x = np.arange(*x_value)
            y = np.arange(*y_value)
            X, Y = np.meshgrid(x, y)
            evaluate_at = np.dstack((X, Y)).transpose(1, 0, 2)
        else:
            evaluate_at = np.array([x_value, y_value])

        z_value = griddata(points, z_array.ravel(), evaluate_at, method='cubic', rescale=True)
        return z_value


    def lookup_expression_from_table(self, *, lengths, vsb, vgs, vds, primary, expression):
        """
        Calculate a parameter using the entire table.
        No interpolation is used.

        Args:
            lengths (float, list, ndarray): length(s) of the mosfet
            vsb (float, tuple): source-body voltage, tuple of the form (start, stop, step)
            vgs (float, tuple): gate-source voltage, tuple of the form (start, stop, step)
            vds (float, tuple): drain-source voltage, tuple of the form (start, stop, step)
            primary (str): name of the primary sweep source: "lengths", "vsb", "vgs", or "vds"
            expression(dict): expression of how to calculate the value you're looking for

        Example:
            x = nmos.lookup_expression_from_table(
                lengths=100e-9,
                vsb=0,
                vds=(0.0, 1, 0.01),
                vgs=(0.0, 1.01, 0.2),
                primary="vds",
                expression=nmos.current_density_expression,
            )
        """
        parameters = expression["variables"].copy()
        remove_from_parameters = ["width", "length"]
        for item in remove_from_parameters:
            if item in parameters:
                parameters.remove(item)
        _, _, extracted_table = self.extract_2d_table(lookup_table=self.lookup_table[self.mos], parameters=parameters, lengths=lengths, vsb=vsb, vgs=vgs, vds=vds, primary=primary)
        x, _ = self.__calculate_from_expression(expression, extracted_table)
        return x
