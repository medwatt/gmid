# imports <<<
from typing import List, Optional, Tuple, Union

import numpy as np

from .expressions import Expression
from .util import evaluate_expression, extract_2d_table
from .interpolation import GridInterpolator, KDTreeInterpolator
from .plot import Plotter
# >>>


class Mosfet:
    """
    Initialize a MOSFET model object.

    When creating a MOSFET object, exactly two of the parameters `length`, `vsb`, `vgs`, and `vds`
    must be provided as fixed values. The remaining parameters can be specified as a range (or set to
    `None` to use the full range from the lookup table). Any parameter not fixed is treated as the
    secondary sweep variable.

    If two parameters are given as ranges, you must indicate which one is the primary sweep variable
    using the `primary` argument.

    Args:
        lookup_table: dictionary containing MOSFET parameter data.
        mos: model name of the MOSFET.
        length: mosfet length.
        vsb: source-body voltage.
        vgs: gate-source voltage.
        vds: drain-source voltage.
        primary: name of the primary sweep.

    Examples:
        Fixed vds; sweep vgs (primary) and length (secondary)
        nmos = Mosfet(lookup_table=lookup_table, mos="nch_lvt", vsb=0.0, vds=0.4, vgs=(0.01, 1.19))
        pmos = Mosfet(lookup_table=lookup_table, mos="pch_lvt", vsb=0.0, vgs=-0.4, vds=(-1.19, -0.01))
    """
    # init <<<
    def __init__(
        self,
        *,
        lookup_table: dict,
        mos: str,
        length: Optional[Union[float, List[float], np.ndarray]] = None,
        vsb: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
        vgs: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
        vds: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
        primary: Optional[str] = None
    ) -> None:
        self.mos = mos
        self.lookup_table = lookup_table[mos]
        self.width = self.lookup_table["width"]
        self.length_all = self.lookup_table["length"]
        self.parameters = self.lookup_table["parameter_names"]

        self.secondary_var, self.filtered_variables, self.extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            length=length,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )
        self.length = self.filtered_variables["length"]
        self.vsb = self.filtered_variables["vsb"]
        self.vgs = self.filtered_variables["vgs"]
        self.vds = self.filtered_variables["vds"]

        self.plotter = Plotter()

        # Initialize basic expressions.
        self.length_expression  = Expression(variables=["length"],  label="$\\mathrm{Length}\\ (m)$")
        self.vsb_expression     = Expression(variables=["vsb"],     label="$V_{\\mathrm{SB}}\\ (V)$")
        self.vgs_expression     = Expression(variables=["vgs"],     label="$V_{\\mathrm{GS}}\\ (V)$")
        self.vds_expression     = Expression(variables=["vds"],     label="$V_{\\mathrm{DS}}\\ (V)$")
        self.id_expression      = Expression(variables=["id"],      label="$I_{D}\\ (A)$")
        self.vth_expression     = Expression(variables=["vth"],     label="$V_{\\mathrm{TH}}\\ (V)$")
        self.vdsat_expression   = Expression(variables=["vdsat"],   label="$V_{\\mathrm{DS_{\\mathrm{SAT}}}}\\ (V)$")
        self.gm_expression      = Expression(variables=["gm"],      label="$g_{m}\\ (S)$")
        self.gmbs_expression    = Expression(variables=["gmbs"],    label="$g_{\\mathrm{mbs}}\\ (S)$")
        self.gds_expression     = Expression(variables=["gds"],     label="$g_{\\mathrm{ds}}\\ (S)$")
        self.cgg_expression     = Expression(variables=["cgg"],     label="$c_{\\mathrm{gg}}\\ (F)$")
        self.cgs_expression     = Expression(variables=["cgs"],     label="$c_{\\mathrm{gs}}\\ (F)$")
        self.cbg_expression     = Expression(variables=["cbg"],     label="$c_{\\mathrm{bg}}\\ (F)$")
        self.cgd_expression     = Expression(variables=["cgd"],     label="$c_{\\mathrm{gd}}\\ (F)$")
        self.cdd_expression     = Expression(variables=["cdd"],     label="$c_{\\mathrm{dd}}\\ (F)$")

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
    # >>>

    # calculate from expression <<<
    def calculate_from_expression(
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
    # >>>

    # plot by expression <<<
    def plot_by_expression(
        self,
        *,
        x_expression: Expression,
        y_expression: Expression,
        filtered_values: Optional[Union[float, List[float], np.ndarray]] = None,
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
        Plots computed x and y expressions, filtering by the secondary sweep variable.
        """
        # Choose filtering variable: secondary if available, else 'length'
        filter_var = self.secondary_var if self.secondary_var is not None else "length"
        sec_values = self.filtered_variables[filter_var]

        if filtered_values is not None:
            if isinstance(filtered_values, (list, np.ndarray)):
                indices = np.nonzero(np.isin(sec_values, np.array(filtered_values)))[0]
            else:
                indices = np.array([np.abs(sec_values - filtered_values).argmin()])
        else:
            indices = np.arange(len(sec_values))

        legend_values = np.array(sec_values)[indices]
        legend_title_mapping = {
            "length": "Length",
            "vsb": "$V_{\\mathrm{SB}}$",
            "vgs": "$V_{\\mathrm{GS}}$",
            "vds": "$V_{\\mathrm{DS}}$"
        }
        legend_title = legend_title_mapping.get(filter_var, filter_var)

        x, x_label = self.calculate_from_expression(x_expression, indices)
        y, y_label = self.calculate_from_expression(y_expression, indices)

        fig, ax = self.plotter.create_figure(
            title="",
            x_label=x_label,
            y_label=y_label,
            x_lim=x_limit,
            y_lim=y_limit,
            x_scale=x_scale,
            y_scale=y_scale,
            x_eng_format=x_eng_format,
            y_eng_format=y_eng_format
        )

        self.plotter.plot_data(
            ax, x, y,
            legend=legend_values,
            legend_title=legend_title,
            save_fig=save_fig,
        )

        if return_result:
            return x, y
        return None
    # >>>

    # plot by sweep <<<
    def plot_by_sweep(
        self,
        *,
        x_expression: Expression,
        y_expression: Expression,
        primary: str,
        length: Optional[Union[float, List[float], np.ndarray]] = None,
        vsb: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
        vgs: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
        vds: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
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
        Exactly two of the parameters `length`, `vsb`, `vgs`, and `vds` must be provided as fixed values.
        The remaining parameters can be specified as a range (or set to `None` to use the full range from the lookup table).
        Any parameter not fixed is treated as the secondary sweep variable.

        If two parameters are given as ranges, you must indicate which one is the primary sweep variable
        using the `primary` argument.

        Args:
            length: Length(s) to filter the lookup table.
            vsb: Source-to-body voltage (or sweep range).
            vgs: Gate-to-source voltage (or sweep range).
            vds: Drain-to-source voltage (or sweep range).
            primary: The primary sweep variable.
            x_expression: Expression for x-axis computation.
            y_expression: Expression for y-axis computation.
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
        secondary_var, filtered_vars, extracted_table = extract_2d_table(
            lookup_table=self.lookup_table,
            width=self.width,
            length=length,
            vsb=vsb,
            vgs=vgs,
            vds=vds,
            primary=primary,
        )

        x, x_label = evaluate_expression(x_expression, extracted_table)
        y, y_label = evaluate_expression(y_expression, extracted_table)
        fig, ax = self.plotter.create_figure(
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
        legend_title = None
        if secondary_var is not None:
            legend = [sw for sw in filtered_vars[secondary_var]]
            legend_title_mapping = {
                "length": "Length",
                "vsb": "$V_{\\mathrm{SB}}$",
                "vgs": "$V_{\\mathrm{GS}}$",
                "vds": "$V_{\\mathrm{DS}}$"
            }
            legend_title = legend_title_mapping.get(secondary_var, secondary_var)

        self.plotter.plot_data(ax, x, y, legend=legend, legend_title=legend_title, save_fig=save_fig)
        if return_result:
            return x, y
        return None
    # >>>

    # quick_plot <<<
    def quick_plot(
        self,
        *,
        x: Union[np.ndarray, List[np.ndarray]],
        y: Union[np.ndarray, List[np.ndarray]],
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

        x and y can be either a single numpy array or a list of numpy arrays.
        This method creates a figure using the plotter and plots the data.

        Args:
            x: x-axis data as a numpy array or a list of numpy arrays.
            y: y-axis data as a numpy array or a list of numpy arrays.
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
    # >>>

    # interpolate <<<
    def interpolate(
        self,
        *,
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
    # >>>

    # lookup expression from table <<<
    def lookup_expression_from_table(
        self,
        *,
        length: Union[float, List[float], np.ndarray],
        vsb: Union[float, Tuple[float, float, float]],
        vgs: Union[float, Tuple[float, float, float]],
        vds: Union[float, Tuple[float, float, float]],
        primary: str,
        expression: Union[Expression, List[Expression]]
    ) -> Union[np.ndarray, List[np.ndarray]]:
        """
        Evaluates one or more expressions on the lookup table without interpolation.

        Args:
            length: The MOSFET length to filter the table.
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
            length=length,
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
    # >>>
