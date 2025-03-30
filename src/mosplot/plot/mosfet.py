import numpy as np
from typing import Tuple, Optional, Union, List
from .helpers import extract_2d_table, evaluate_expression
from .plot import Plotter
from .expressions import Expression
from .interpolation import GridInterpolator, KDTreeInterpolator

class Mosfet:
    """
    Initialize a mos object.
    Two of `lengths, vsb, vgs, vds` must be fixed at any time.

    Args:
        lookup_table: dictionary of mosfet parameters
        mos: model name of the mos transistor
        lengths: length(s) of the mosfet
        vsb: source-body voltage
        vgs: gate-source voltage
        vds: drain-source voltage
        primary: name of the primary sweep source

    Example:
        nmos = Mosfet(lookup_table=lookup_table, mos="nch_lvt", vsb=0.0, vds=0.4, vgs=(0.01, 1.19))
        pmos = Mosfet(lookup_table=lookup_table, mos="pch_lvt", vsb=0.0, vds=-0.4, vgs=(-1.19, -0.01))
    """
    def __init__(
        self,
        lookup_table: dict,
        mos: str,
        lengths: Optional[Union[float, List[float], np.ndarray]] = None,
        vsb: Optional[Union[float, Tuple[float, float, float]]] = None,
        vgs: Optional[Union[float, Tuple[float, float, float]]] = None,
        vds: Optional[Union[float, Tuple[float, float, float]]] = None,
        primary: Optional[str] = None
    ) -> None:
        self.mos = mos
        self.lookup_table = lookup_table[mos]
        self.width = self.lookup_table["width"]
        self.lengths_all = self.lookup_table["lengths"]
        self.parameters = self.lookup_table["parameter_names"]

        self.secondary_idx, self.filtered_variables, self.extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            lengths=lengths,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )
        self.lengths = self.filtered_variables[0]
        self.vsb = self.filtered_variables[1]
        self.vgs = self.filtered_variables[2]
        self.vds = self.filtered_variables[3]

        self.plotter = Plotter()

        # Initialize basic expressions.
        self.lengths_expression = Expression(variables=["lengths"], label ="$\\mathrm{Length}\\ (m)$")
        self.vsb_expression     = Expression(variables=["vsb"],     label ="$V_{\\mathrm{SB}}\\ (V)$")
        self.vgs_expression     = Expression(variables=["vgs"],     label ="$V_{\\mathrm{GS}}\\ (V)$")
        self.vds_expression     = Expression(variables=["vds"],     label ="$V_{\\mathrm{DS}}\\ (V)$")
        self.id_expression      = Expression(variables=["id"],      label ="$I_{D}\\ (A)$")
        self.vth_expression     = Expression(variables=["vth"],     label ="$V_{\\mathrm{TH}}\\ (V)$")
        self.vdsat_expression   = Expression(variables=["vdsat"],   label ="$V_{\\mathrm{DS_{\\mathrm{SAT}}}}\\ (V)$")
        self.gm_expression      = Expression(variables=["gm"],      label ="$g_{m}\\ (S)$")
        self.gmbs_expression    = Expression(variables=["gmbs"],    label ="$g_{\\mathrm{mbs}}\\ (S)$")
        self.gds_expression     = Expression(variables=["gds"],     label ="$g_{\\mathrm{ds}}\\ (S)$")
        self.cgg_expression     = Expression(variables=["cgg"],     label ="$c_{\\mathrm{gg}}\\ (F)$")
        self.cgs_expression     = Expression(variables=["cgs"],     label ="$c_{\\mathrm{gs}}\\ (F)$")
        self.cbg_expression     = Expression(variables=["cbg"],     label ="$c_{\\mathrm{bg}}\\ (F)$")
        self.cgd_expression     = Expression(variables=["cgd"],     label ="$c_{\\mathrm{gd}}\\ (F)$")
        self.cdd_expression     = Expression(variables=["cdd"],     label ="$c_{\\mathrm{dd}}\\ (F)$")

        # Initialize computed expressions.
        self.gmid_expression = Expression(
            variables=["gm", "id"],
            function=lambda x, y: x / y,
            label="$g_m/I_D\\ (S/A)$"
        )
        self.vstar_expression = Expression(
            variables=["gm", "id"],
            function=lambda x, y: (2 * y) / x,
            label="$V^{\\star}\\ (V)$"
        )
        self.gain_expression = Expression(
            variables=["gm", "gds"],
            function=lambda x, y: x / y,
            label="$g_{m}/g_{\\mathrm{ds}}$"
        )
        self.current_density_expression = Expression(
            variables=["id", "width"],
            function=lambda x, y: x / y,
            label="$I_{D}/W\\ (A/m)$"
        )
        self.transist_frequency_expression = Expression(
            variables=["gm", "cgg"],
            function=lambda x, y: x / (2 * np.pi * y),
            label="$f_{T}\\ (Hz)$"
        )
        self.early_voltage_expression = Expression(
            variables=["id", "gds"],
            function=lambda x, y: x / y,
            label="$V_{A}\\ (V)$"
        )
        self.rds_expression = Expression(
            variables=["gds"],
            function=lambda x: 1 / x,
            label="$r_{\\mathrm{ds}}\\ (\\Omega)$"
        )

    def _calculate_from_expression(
        self,
        expression: Expression,
        filter_by_rows: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, str]:
        """
        Calculates values from a given expression over the extracted lookup table.

        Args:
            expression: The Expression object defining the computation.
            filter_by_rows: Optional array of row indices to filter the data.

        Returns:
            A tuple containing the computed numpy array and the expression's label.
        """
        return evaluate_expression(expression, self.extracted_table, filter_by_rows)

    def plot_by_expression(
        self,
        x_expression: Expression,
        y_expression: Expression,
        lengths: Optional[Union[float, List[float], np.ndarray]] = None,
        x_limit: Optional[Tuple[float, float]] = None,
        y_limit: Optional[Tuple[float, float]] = None,
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        title: Optional[str] = None,
        save_fig: str = "",
        return_result: bool = False
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Plots data computed from expressions for the x and y axes.

        Args:
            x_expression: Expression for computing x-axis values.
            y_expression: Expression for computing y-axis values.
            lengths: Optional length filter.
            x_limit: Optional limits for the x-axis.
            y_limit: Optional limits for the y-axis.
            x_scale: Optional x-axis scale type (e.g., 'log').
            y_scale: Optional y-axis scale type.
            x_eng_format: If True, format x-axis labels in engineering units.
            y_eng_format: If True, format y-axis labels in engineering units.
            title: Optional title for the plot.
            save_fig: Optional filename to save the figure.
            return_result: If True, return the computed x and y arrays.

        Returns:
            A tuple (x, y) of numpy arrays if return_result is True; otherwise, None.
        """
        title_label: List[str] = []
        variables_labels = ["lengths", "vsb", "vgs", "vds"]
        model_name = self.lookup_table.get("model_name", "")

        for i, v in enumerate(self.filtered_variables):
            if not isinstance(v, np.ndarray):
                label = variables_labels[i]
                title_label.extend([f"V_{{\\mathrm{{{label[1:].upper()}}}}}", str(v)])

        plot_title = title if title is not None else (
            f"{model_name}, " + "$%s=%s$, $%s=%s$" % tuple(title_label)
        )

        indices = (
            np.nonzero(np.isin(self.lengths, np.array(lengths)))[0]
            if lengths is not None
            else np.arange(len(self.lengths))
        )

        legend = [f"{length}" for length in np.array(self.lengths)[indices]]

        x, x_label = self._calculate_from_expression(x_expression, indices)
        y, y_label = self._calculate_from_expression(y_expression, indices)

        fig, ax = self.plotter.create_figure(
            title=plot_title,
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit,
            y_lim=y_limit,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format
        )

        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)

        if return_result:
            return x, y

        return None

    def plot_by_sweep(
        self,
        lengths: Union[float, List[float], np.ndarray],
        vsb: Union[float, Tuple[float, float, float]],
        vgs: Union[float, Tuple[float, float, float]],
        vds: Union[float, Tuple[float, float, float]],
        primary: str,
        x_expression: Expression,
        y_expression: Expression,
        title: str = "",
        x_limit: Optional[Tuple[float, float]] = None,
        y_limit: Optional[Tuple[float, float]] = None,
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        save_fig: str = "",
        return_result: bool = False
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Plots data by sweeping the parameters.

        Args:
            lengths: Length(s) to filter the lookup table.
            vsb: Source-to-body voltage (or sweep range).
            vgs: Gate-to-source voltage (or sweep range).
            vds: Drain-to-source voltage (or sweep range).
            primary: The primary sweep variable.
            x_expression: Expression for x-axis computation.
            y_expression: Expression for y-axis computation.
            title: Title of the plot.
            x_limit: Optional x-axis limits.
            y_limit: Optional y-axis limits.
            x_scale: Optional x-axis scale type.
            y_scale: Optional y-axis scale type.
            x_eng_format: If True, format x-axis labels in engineering units.
            y_eng_format: If True, format y-axis labels in engineering units.
            save_fig: Optional filename to save the figure.
            return_result: If True, return computed (x, y) arrays.

        Returns:
            A tuple (x, y) if return_result is True; otherwise, None.
        """
        sec_idx, filtered_vars, extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            lengths=lengths,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )
        x, x_label = evaluate_expression(x_expression, extracted_table)
        y, y_label = evaluate_expression(y_expression, extracted_table)
        _, ax = self.plotter.create_figure(
            title=title,
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit,
            y_lim=y_limit,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format
        )
        legend = None
        if sec_idx is not None:
            legend = [str(sw) for sw in filtered_vars[sec_idx]]
        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)
        if return_result:
            return x, y
        return None

    def quick_plot(
        self,
        x: np.ndarray,
        y: np.ndarray,
        x_label: str = "",
        y_label: str = "",
        x_limit: Optional[Tuple[float, float]] = None,
        y_limit: Optional[Tuple[float, float]] = None,
        x_scale: str = "",
        y_scale: str = "",
        x_eng_format: bool = False,
        y_eng_format: bool = False,
        legend: Optional[List[str]] = None,
        title: Optional[str] = None,
        save_fig: str = ""
    ) -> None:
        """
        Quickly plots the provided x and y data.

        Args:
            x: x-axis data as a numpy array.
            y: y-axis data as a numpy array.
            x_label: Label for the x-axis.
            y_label: Label for the y-axis.
            x_limit: Optional x-axis limits.
            y_limit: Optional y-axis limits.
            x_scale: Optional x-axis scale type.
            y_scale: Optional y-axis scale type.
            x_eng_format: If True, format x-axis labels in engineering units.
            y_eng_format: If True, format y-axis labels in engineering units.
            legend: Optional list of legend entries.
            title: Optional title for the plot.
            save_fig: Optional filename to save the figure.
        """
        _, ax = self.plotter.create_figure(
            title=title or "",
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit,
            y_lim=y_limit,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format
        )
        self.plotter.plot_data(ax, x, y, legend=legend, save_fig=save_fig)

    def interpolate(
        self,
        x_expression: Expression,
        x_value: Union[float, Tuple[float, float], np.ndarray],
        y_expression: Expression,
        y_value: Union[float, Tuple[float, float], np.ndarray],
        z_expression: Union[Expression, List[Expression]],
        fast: bool = False
    ) -> Union[np.ndarray, List[np.ndarray]]:
        """
        Interpolates z values given x and y expressions.

        Args:
            x_expression: Expression for the x-axis.
            x_value: Value(s) for the x-axis interpolation domain.
            y_expression: Expression for the y-axis.
            y_value: Value(s) for the y-axis interpolation domain.
            z_expression: Expression or list of expressions for the z values.
            fast: If True, use KDTree-based interpolation; otherwise, use griddata.

        Returns:
            A single interpolated numpy array if z_expression is a single Expression,
            or a list of interpolated numpy arrays if a list is provided.
        """
        if not fast:
            interpolator = GridInterpolator(self.extracted_table)
        else:
            interpolator = KDTreeInterpolator(self.extracted_table)
        return interpolator.interpolate(x_expression, x_value, y_expression, y_value, z_expression)

    def lookup_expression_from_table(
        self,
        lengths: Union[float, List[float], np.ndarray],
        vsb: Union[float, Tuple[float, float, float]],
        vgs: Union[float, Tuple[float, float, float]],
        vds: Union[float, Tuple[float, float, float]],
        primary: str,
        expression: Union[Expression, List[Expression]]
    ) -> Union[np.ndarray, List[np.ndarray]]:
        """
        Evaluates one or more expressions on the lookup table without interpolation.

        Args:
            lengths: The MOSFET lengths to filter the table.
            vsb: The source-body voltage parameter.
            vgs: The gate-source voltage parameter.
            vds: The drain-source voltage parameter.
            primary: The primary sweep variable.
            expression: A single Expression or a list of Expressions to evaluate.

        Returns:
            A numpy array if a single Expression is provided, or a list of numpy arrays if multiple
            Expressions are provided.
        """
        _, _, extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            lengths=lengths,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )

        def process_expr(expr: Expression) -> np.ndarray:
            result, _ = evaluate_expression(expr, extracted_table)
            return result

        if isinstance(expression, list):
            return [process_expr(expr) for expr in expression]
        else:
            return process_expr(expression)

