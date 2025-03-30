# imports <<<
import numpy as np
from scipy.interpolate import griddata
from .helpers import extract_2d_table
from .plotter import Plotter
from scipy.spatial import cKDTree
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
        z_expression,  # can be a dict or list of dicts
    ):
        """
        Interpolate using the given expressions. If z_expression is a list,
        returns a list of results.

        Args:
            x_expression (dict): Expression for x-axis points.
            x_value (float, tuple, or np.ndarray): Value(s) within x_expression domain.
            y_expression (dict): Expression for y-axis points.
            y_value (float, tuple, or np.ndarray): Value(s) within y_expression domain.
            z_expression (dict or list of dict): Expression(s) for the output value.

        Returns:
            Interpolated value or list of interpolated values.
        """
        x_array, _ = self._calculate_from_expression(x_expression, self.extracted_table)
        y_array, _ = self._calculate_from_expression(y_expression, self.extracted_table)
        points = np.column_stack((x_array.ravel(), y_array.ravel()))

        if isinstance(x_value, (tuple, np.ndarray)) and isinstance(y_value, (int, float)):
            x_vals = np.arange(*x_value) if isinstance(x_value, tuple) else x_value
            eval_points = np.column_stack((x_vals, np.full(np.shape(x_vals), y_value)))
        elif isinstance(y_value, (tuple, np.ndarray)) and isinstance(x_value, (int, float)):
            y_vals = np.arange(*y_value) if isinstance(y_value, tuple) else y_value
            eval_points = np.column_stack((np.full(np.shape(y_vals), x_value), y_vals))
        elif isinstance(x_value, tuple) and isinstance(y_value, tuple):
            x_vals = np.arange(*x_value)
            y_vals = np.arange(*y_value)
            X, Y = np.meshgrid(x_vals, y_vals)
            eval_points = np.dstack((X, Y)).transpose(1, 0, 2)
        else:
            eval_points = np.array([x_value, y_value])

        single_input = False
        if not isinstance(z_expression, list):
            z_expressions = [z_expression]
            single_input = True
        else:
            z_expressions = z_expression

        results = []
        for expr in z_expressions:
            z_array, _ = self._calculate_from_expression(expr, self.extracted_table)
            res = griddata(points, z_array.ravel(), eval_points, method="cubic", rescale=True)
            results.append(res)
        return results[0] if single_input else results
    # >>>

    # fast interpolate <<<
    def fast_interpolate(
        self,
        x_expression: dict,
        x_value,
        y_expression: dict,
        y_value,
        z_expression,
    ):
        """
        Fast interpolation using cKDTree and inverse distance weighting.

        Args:
            x_expression (dict): Expression for x-axis points.
            x_value (float, tuple, or np.ndarray): Value(s) within x_expression domain.
            y_expression (dict): Expression for y-axis points.
            y_value (float, tuple, or np.ndarray): Value(s) within y_expression domain.
            z_expression (dict or list of dict): Expression(s) for the output value.

        Returns:
            Interpolated value or list of interpolated values.
        """
        # Compute original arrays and build points
        x_array, _ = self._calculate_from_expression(x_expression, self.extracted_table)
        y_array, _ = self._calculate_from_expression(y_expression, self.extracted_table)
        points = np.column_stack((x_array.ravel(), y_array.ravel()))

        # Scale coordinates
        x_min, x_max = x_array.min(), x_array.max()
        y_min, y_max = y_array.min(), y_array.max()
        scaled_points = np.column_stack(((points[:, 0] - x_min) / (x_max - x_min),
                                         (points[:, 1] - y_min) / (y_max - y_min)))
        tree = cKDTree(scaled_points)

        # Build evaluation points
        if isinstance(x_value, (tuple, np.ndarray)) and isinstance(y_value, (int, float)):
            x_vals = np.arange(*x_value) if isinstance(x_value, tuple) else x_value
            eval_points = np.column_stack((x_vals, np.full(np.shape(x_vals), y_value)))
        elif isinstance(y_value, (tuple, np.ndarray)) and isinstance(x_value, (int, float)):
            y_vals = np.arange(*y_value) if isinstance(y_value, tuple) else y_value
            eval_points = np.column_stack((np.full(np.shape(y_vals), x_value), y_vals))
        elif isinstance(x_value, tuple) and isinstance(y_value, tuple):
            x_vals = np.arange(*x_value)
            y_vals = np.arange(*y_value)
            X, Y = np.meshgrid(x_vals, y_vals)
            eval_points = np.dstack((X, Y)).transpose(1, 0, 2)
        else:
            eval_points = np.array([x_value, y_value])

        # Prepare evaluation points for tree query
        if eval_points.ndim == 1:
            eval_points = eval_points.reshape(1, 2)
            out_shape = None
        elif eval_points.ndim == 2:
            out_shape = None
        elif eval_points.ndim > 2 and eval_points.shape[-1] == 2:
            out_shape = eval_points.shape[:-1]
            eval_points = eval_points.reshape(-1, 2)
        else:
            out_shape = None

        # Scale evaluation points
        scaled_eval_points = np.column_stack(((eval_points[:, 0] - x_min) / (x_max - x_min),
                                              (eval_points[:, 1] - y_min) / (y_max - y_min)))

        # Prepare z expressions
        single_input = False
        if not isinstance(z_expression, list):
            z_expressions = [z_expression]
            single_input = True
        else:
            z_expressions = z_expression

        results = []
        k = min(8, len(scaled_points))  # Number of neighbors
        eps = 1e-12

        for expr in z_expressions:
            z_array, _ = self._calculate_from_expression(expr, self.extracted_table)
            z_flat = z_array.ravel()
            dist, idx = tree.query(scaled_eval_points, k=k)

            # Ensure 2D arrays when k == 1
            if k == 1:
                dist = dist[:, None]
                idx = idx[:, None]

            res = np.empty(len(scaled_eval_points))
            zero_mask = dist[:, 0] < eps
            res[zero_mask] = z_flat[idx[zero_mask, 0]]
            if np.any(~zero_mask):
                d = dist[~zero_mask]
                i = idx[~zero_mask]
                weights = 1 / (d ** 2)
                res[~zero_mask] = np.sum(weights * z_flat[i], axis=1) / np.sum(weights, axis=1)

            if out_shape is not None:
                res = res.reshape(out_shape)
            results.append(res)

        return results[0] if single_input else results
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
