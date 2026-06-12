from __future__ import annotations

from pathlib import Path


def export_solve_result_csv(result, output_path: str) -> None:
    import csv

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["case", result.case_name])
        writer.writerow(["converged", result.converged])
        writer.writerow([])
        writer.writerow(["Node", "Pressure (Pa)", "Pressure (kPa)"])
        for node_id in sorted(result.node_pressures_pa):
            pressure_pa = result.node_pressures_pa[node_id]
            writer.writerow([node_id, round(pressure_pa, 2), round(pressure_pa / 1000.0, 4)])
        writer.writerow([])
        writer.writerow(["Component", "Mass flow (kg/s)", "Vol. flow (m^3/h)"])
        for component in result.component_flows:
            writer.writerow(
                [
                    component.label,
                    round(component.mass_flow_kg_per_s, 4),
                    round(component.volumetric_flow_m3_per_h, 4),
                ]
            )


def export_solve_result_workbook(result, output_path: str) -> None:
    from openpyxl import Workbook

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    workbook = Workbook()

    pressures_sheet = workbook.active
    pressures_sheet.title = "Pressures"
    pressures_sheet.append(["Node", "Pressure (Pa)", "Pressure (kPa)"])
    for node_id in sorted(result.node_pressures_pa):
        pressure_pa = result.node_pressures_pa[node_id]
        pressures_sheet.append([node_id, round(pressure_pa, 2), round(pressure_pa / 1000.0, 4)])

    flows_sheet = workbook.create_sheet("Flows")
    flows_sheet.append(["Component", "Mass flow (kg/s)", "Vol. flow (m^3/h)"])
    for component in result.component_flows:
        flows_sheet.append(
            [
                component.label,
                round(component.mass_flow_kg_per_s, 4),
                round(component.volumetric_flow_m3_per_h, 4),
            ]
        )

    workbook.save(output)


def print_solve_result(result, detailed: bool = False) -> None:
    print(f"Case: {result.case_name}")
    print(f"Converged: {result.converged}")
    print("")
    _print_iteration_summary(result)
    print("")
    _print_node_table(result.node_pressures_pa)
    print("")
    _print_component_table(result.component_flows)
    if result.turbulent_metrics:
        print("")
        _print_final_turbulent_metrics(result.turbulent_metrics[-1])
    if detailed:
        print("")
        _print_history_section("Laminar Initialisation History", result.laminar_history)
        if result.laminar_metrics:
            print("")
            _print_iteration_metrics_history("Laminar Iteration Metrics", result.laminar_metrics)
        print("")
        _print_history_section("Turbulent Correction History", result.turbulent_history)
        if result.turbulent_metrics:
            print("")
            _print_iteration_metrics_history("Turbulent Iteration Metrics", result.turbulent_metrics)


def print_detailed_solve_result(result) -> None:
    print_solve_result(result, detailed=True)


def save_convergence_plot(result, output_path: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=False)
    _plot_history(axes[0], "Laminar correction history", result.laminar_history)
    _plot_history(axes[1], "Turbulent correction history", result.turbulent_history)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def _print_node_table(node_pressures_pa) -> None:
    print("Node Pressures")
    print(f"{'Node':>6}  {'Pressure (Pa)':>14}  {'Pressure (kPa)':>14}")
    for node_id in sorted(node_pressures_pa):
        pressure_pa = node_pressures_pa[node_id]
        print(f"{node_id:>6}  {pressure_pa:>14.2f}  {pressure_pa / 1000.0:>14.4f}")


def _print_component_table(component_flows) -> None:
    print("Component Flows")
    print(f"{'Component':<20}  {'Mass flow (kg/s)':>16}  {'Vol. flow (m^3/h)':>18}")
    for component in component_flows:
        print(
            f"{component.label:<20}  "
            f"{component.mass_flow_kg_per_s:>16.4f}  "
            f"{component.volumetric_flow_m3_per_h:>18.4f}"
        )


def _print_iteration_summary(result) -> None:
    print("Iteration Summary")
    print(f"  Laminar initialisation iterations: {len(result.laminar_history)}")
    print(f"  Turbulent iterations:              {len(result.turbulent_history)}")


def _print_history_section(title: str, history) -> None:
    print(title)
    print(_ascii_history_plot(history))
    print(f"{'Iter':>6}  {'Correction':>14}")
    for iteration, value in enumerate(history, start=1):
        print(f"{iteration:>6}  {value:>14.9f}")


def _print_final_turbulent_metrics(metrics) -> None:
    print("Final Turbulent Metrics")
    print(f"  Max abs pressure correction (Pa): {metrics.pressure_correction_abs_pa:.9f}")
    print(f"  Max rel pressure correction:      {metrics.pressure_correction_rel:.9e}")
    print(f"  Max nodal mass imbalance (−):     {metrics.max_nodal_mass_imbalance_rel:.9e}")


def _print_iteration_metrics_history(title: str, metrics_history) -> None:
    print(title)
    print(
        f"{'Iter':>6}  "
        f"{'Abs Corr. (Pa)':>16}  "
        f"{'Rel Corr.':>12}  "
        f"{'Max Mass Imb. (−)':>20}"
    )
    for iteration, metrics in enumerate(metrics_history, start=1):
        print(
            f"{iteration:>6}  "
            f"{metrics.pressure_correction_abs_pa:>16.9f}  "
            f"{metrics.pressure_correction_rel:>12.5e}  "
            f"{metrics.max_nodal_mass_imbalance_rel:>20.9e}"
        )


def _ascii_history_plot(history, width: int = 48) -> str:
    if not history:
        return "(no data)"

    blocks = " .:-=+*#%@"
    minimum = min(history)
    maximum = max(history)
    span = maximum - minimum
    if span == 0:
        line = blocks[-1] * len(history)
    else:
        indices = []
        for value in history:
            normalised = (value - minimum) / span
            indices.append(int(round(normalised * (len(blocks) - 1))))
        line = "".join(blocks[index] for index in indices)

    if len(line) > width:
        step = len(line) / float(width)
        compressed = []
        for idx in range(width):
            compressed.append(line[int(idx * step)])
        line = "".join(compressed)

    return f"[{line}]  min={minimum:.6g} max={maximum:.6g}"


def _plot_history(axis, title: str, history) -> None:
    iterations = list(range(1, len(history) + 1))
    axis.plot(iterations, history, marker="o", linewidth=1.5)
    axis.set_title(title)
    axis.set_xlabel("Iteration")
    axis.set_ylabel("Correction")
    axis.grid(True, alpha=0.3)
