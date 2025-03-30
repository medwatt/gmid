class DesignReport:
    def __init__(self, circuit, optimizer):
        self.circuit = circuit
        self.optimizer = optimizer
        self.opt_params = optimizer.get_opt_params()
        self.full_specs = circuit.evaluate_specs(**self.opt_params)
        self.objective = optimizer.result.fun

    @staticmethod
    def format_eng(value):
        if value is None:
            return "N/A"
        if value < 1e-6:
            return f"{int(round(value * 1e9))} nm"
        return f"{value * 1e6:.2f} um"

    @staticmethod
    def format_area(value):
        if value is None:
            return "N/A"
        # Assume area is in m²; convert to µm² or nm².
        area_um2 = value * 1e12
        if area_um2 < 1:
            area_nm2 = value * 1e18
            return f"{area_nm2:.3f} nm²"
        return f"{area_um2:.3f} um²"

    def report(self):
        lines = []
        dims = getattr(self.circuit, "device_dimensions", None)
        lines.append("\nTransistor Details:")
        if dims:
            for device, item in dims.items():
                lines.append(f"    {device}:")
                length = item.get("Length")
                width = item.get("Width")
                area = item.get("Area")
                current = item.get("Current")
                gmid = item.get("GMID")
                lines.append(f"        Length: {self.format_eng(length)}")
                lines.append(f"        Width: {self.format_eng(width)}")
                lines.append(f"        Area: {self.format_area(area)}")
                lines.append(f"        Current: {self.format_eng(current)}")
                if gmid is not None:
                    lines.append(f"        GmID: {gmid:.2f}")
        else:
            lines.append("    Not available.")

        lines.append("\nFinal Circuit Specs:")
        if self.full_specs:
            for key, val in self.full_specs.items():
                if key == "Area":
                    lines.append(f"    {key}: {self.format_area(val)}")
                else:
                    try:
                        lines.append(f"    {key}: {val:.4g}")
                    except (ValueError, TypeError):
                        lines.append(f"    {key}: {val}")
        else:
            lines.append("    Not available.")
        return "\n".join(lines)
