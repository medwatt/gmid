# imports <<<
import numpy as np
from scipy.interpolate import griddata
from .helpers import extract_2d_table
from .plotter import Plotter
# >>>


class Mosfet:
    # init <<<
    def __init__(
        self,
        lookup_table: dict,
        mos: str,
        lengths=None,
        vsb=None,
        vgs=None,
        vds=None,
        primary=None,
    ):
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
        self.mos = mos
        self.lookup_table = lookup_table[mos]
        self.width = self.lookup_table["width"]
        self.lengths_all = self.lookup_table["lengths"]
        self.parameters = self.lookup_table["parameter_names"]

        self.secondary_idx, self.filtered_variables, self.extracted_table = (
            extract_2d_table(
                lookup_table=self.lookup_table,
                width=self.width,
                lengths=lengths,
                vsb=vsb,
                vgs=vgs,
                vds=vds,
                primary=primary,
            )
        )
        self.lengths, self.vsb, self.vgs, self.vds = self.filtered_variables

        # define commonly-used expressions/methods dynamically
        self._setup_common_expressions()
        self._settup_common_plot_methods()

        # initalize plotter
        self.plotter = Plotter()

    # >>>

    # setup common expressions <<<
    def _setup_common_expressions(self):
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
            if parameter in self.parameters or parameter in [
                "lengths",
                "vsb",
                "vgs",
                "vds",
            ]:
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

    # >>>

    # plot by expression <<<
    def plot_by_expression(
        self,
        *,
        x_expression: dict,
        y_expression: dict,
        lengths=(),
        x_limit=(),
        y_limit=(),
        x_scale="",
        y_scale="",
        x_eng_format=False,
        y_eng_format=False,
        title=None,
        save_fig="",
        return_result=False,
    ):
        title_label = []
        variables_labels = ["lengths", "vsb", "vgs", "vds"]
        model_name = self.lookup_table.get("model_name", "")
        for i, v in enumerate(self.filtered_variables):
            if not isinstance(v, np.ndarray):
                label = variables_labels[i]
                title_label.extend([f"V_{{\\mathrm{{{label[1:].upper()}}}}}", v])
        plot_title = (
            f"{model_name}, " + "$%s=%.2f$, $%s=%.2f$" % tuple(title_label)
            if title is None
            else title
        )
        indices = (
            np.nonzero(np.isin(self.lengths, np.array(lengths)))[0]
            if lengths
            else np.arange(len(self.lengths))
        )
        legend = (
            [f"{length}" for length in np.array(self.lengths)[indices]]
            if indices.size > 0
            else [f"{length}" for length in self.lengths]
        )
        x, x_label = self._calculate_from_expression(
            x_expression, self.extracted_table, indices
        )
        y, y_label = self._calculate_from_expression(
            y_expression, self.extracted_table, indices
        )
        fig, ax = self.plotter.create_figure(
            title=plot_title,
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit if x_limit else None,
            y_lim=y_limit if y_limit else None,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format,
        )
        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)
        if return_result:
            return x, y

    # >>>

    # plot by sweep <<<
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
        title="",
        x_limit=(),
        y_limit=(),
        x_scale="",
        y_scale="",
        x_eng_format=False,
        y_eng_format=False,
        save_fig="",
        return_result=False,
    ):
        secondary_idx, filtered_variables, extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            lengths=lengths,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )

        x, x_label = self._calculate_from_expression(x_expression, extracted_table)
        y, y_label = self._calculate_from_expression(y_expression, extracted_table)
        fig, ax = self.plotter.create_figure(
            title=title,
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit if x_limit else None,
            y_lim=y_limit if y_limit else None,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format,
        )
        legend = None

        if secondary_idx is not None:
            legend = [str(sw) for sw in filtered_variables[secondary_idx]]
        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)
        if return_result:
            return x, y

    # >>>

    # quick plot <<<
    def quick_plot(
        self,
        x,
        y,
        x_label="",
        y_label="",
        x_limit=(),
        y_limit=(),
        x_scale="",
        y_scale="",
        x_eng_format=False,
        y_eng_format=False,
        legend=None,
        title=None,
        save_fig="",
    ):
        """
        Make quick plots. As a reminder, when `x` and `y` are of size m x n, pass
        them to this function as x.T and y.T
        """
        fig, ax = self.plotter.create_figure(
            title=title if title else "",
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit if x_limit else None,
            y_lim=y_limit if y_limit else None,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format,
        )
        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)

    # >>>

    # common plot methods <<<
    def _settup_common_plot_methods(self):
        PLOT_METHODS = {
            "current_density_plot": [
                self.gmid_expression,
                self.current_density_expression,
            ],
            "gain_plot": [self.gmid_expression, self.gain_expression],
            "transit_frequency_plot": [
                self.gmid_expression,
                self.transist_frequency_expression,
            ],
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
                title: str = "",
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
            setattr(
                Mosfet,
                method_name,
                create_plot_method(self, x_expression=x, y_expression=y),
            )

    # >>>

    # calculate from expression <<<
    def _calculate_from_expression(
        self, expression: dict, table: dict, filter_by_rows=None
    ):
        if filter_by_rows is None:
            filter_by_rows = np.array([])
        if isinstance(expression, dict):
            var_list = []
            for v in expression["variables"]:
                var = table[v]
                if filter_by_rows.size > 0 and var.ndim > 1:
                    var_list.append(np.take(var, filter_by_rows, 0))
                else:
                    var_list.append(var)
            values = expression.get("function", lambda x: x)(*var_list)
            return values, expression.get("label", "")
        return expression, None

    # >>>

    # interpolate <<<
    def interpolate(
        self,
        x_expression: dict,
        x_value,
        y_expression: dict,
        y_value,
        z_expression: dict,
    ):
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
        x_array, _ = self._calculate_from_expression(x_expression, self.extracted_table)
        y_array, _ = self._calculate_from_expression(y_expression, self.extracted_table)
        z_array, _ = self._calculate_from_expression(z_expression, self.extracted_table)
        points = np.column_stack((x_array.ravel(), y_array.ravel()))
        if isinstance(x_value, (tuple, np.ndarray)) and isinstance(
            y_value, (int, float)
        ):
            x_vals = np.arange(*x_value) if isinstance(x_value, tuple) else x_value
            eval_points = np.column_stack((x_vals, np.full(np.shape(x_vals), y_value)))
        elif isinstance(y_value, (tuple, np.ndarray)) and isinstance(
            x_value, (int, float)
        ):
            y_vals = np.arange(*y_value) if isinstance(y_value, tuple) else y_value
            eval_points = np.column_stack((np.full(np.shape(y_vals), x_value), y_vals))
        elif isinstance(x_value, tuple) and isinstance(y_value, tuple):
            x_vals = np.arange(*x_value)
            y_vals = np.arange(*y_value)
            X, Y = np.meshgrid(x_vals, y_vals)
            eval_points = np.dstack((X, Y)).transpose(1, 0, 2)
        else:
            eval_points = np.array([x_value, y_value])
        return griddata(
            points, z_array.ravel(), eval_points, method="cubic", rescale=True
        )

    # >>>

    # lookup expression from table <<<
    def lookup_expression_from_table(
        self, lengths, vsb, vgs, vds, primary, expression: dict
    ):
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
        params = expression["variables"].copy()
        for item in ["width", "length"]:
            if item in params:
                params.remove(item)
        _, _, extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            lengths=lengths,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )
        x, _ = self._calculate_from_expression(expression, extracted_table)
        return x

    # >>>
