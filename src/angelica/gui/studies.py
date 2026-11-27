"""
GUI for sensitivity studies and nodal analysis (``angelica.studies``).

``study_targets`` and ``study_outputs`` list what can be varied and what can
be plotted for a case built from the scene; the dialogs solve and draw a
line plot on a Tk canvas (axes from ``plot_axes``).
"""
from __future__ import annotations

import dataclasses
import math
import tkinter as tk
from tkinter import messagebox, ttk

from .plot_axes import nice_ticks, tick_label, tick_values

_SKIP_FIELDS = {"n_thermal_segments", "tubing_segments", "n_steps", "num_segments"}


def study_targets(case) -> list[tuple[str, str]]:
    """(label, target) pairs of the numeric parameters of a case."""
    out: list[tuple[str, str]] = []
    for b in tuple(case.pressure_inlets) + tuple(case.pressure_outlets):
        out.append((f"Pressure of node {b.node_id} (Pa)", f"pressure:{b.node_id}"))
    for b in tuple(case.flow_inlets) + tuple(case.flow_outlets):
        out.append((f"Mass flow at node {b.node_id} (kg/s)", f"flow:{b.node_id}"))
    for t in case.thermal_inlets:
        if t.bc_type == "fixed_temperature":
            out.append((f"Temperature of node {t.node_id} (°C)", f"temperature:{t.node_id}"))
    for c in case.components:
        for f in dataclasses.fields(c):
            if f.name in _SKIP_FIELDS or f.name in ("start_node", "end_node"):
                continue
            v = getattr(c, f.name)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(float(v)):
                continue
            out.append((f"{c.component_id}: {f.name}", f"component:{c.component_id}.{f.name}"))
    return out


def study_outputs(case) -> list[tuple[str, tuple[str, object]]]:
    """(label, (kind, key)) pairs of what a study can plot."""
    out: list[tuple[str, tuple[str, object]]] = [
        (f"Mass flow through {c.component_id} (kg/s)", ("flow", c.component_id)) for c in case.components]
    nodes = sorted({n for c in case.components for n in (c.start_node, c.end_node)})
    out += [(f"Pressure at node {n} (Pa)", ("pressure", n)) for n in nodes]
    return out


def output_value(point, output: tuple[str, object]) -> float | None:
    """The plotted quantity of a study point (None when it did not converge)."""
    if not point.converged:
        return None                                  # not converged, or no result (see point.error)
    kind, key = output
    return point.flow(str(key)) if kind == "flow" else point.pressure(int(key))   # type: ignore[call-overload]


def draw_plot(canvas: tk.Canvas, series, x_label: str, y_label: str, width: int = 640, height: int = 400,
              markers=()) -> None:
    """Line plot of ``series`` = [(name, colour, [(x, y), …]), …]; ``markers``
    = [(x, y, colour, text), …]."""
    canvas.delete("all")
    series = [(n, c, [(x, y) for x, y in s if y is not None and math.isfinite(x) and math.isfinite(y)])
              for n, c, s in series]
    pts = [(x, y) for _n, _c, s in series for x, y in s] + [(m[0], m[1]) for m in markers]
    if not pts:
        canvas.create_text(width / 2, height / 2, text="No converged point to plot")
        return
    ML, MR, MT, MB = 80, 20, 20, 50
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, x1, xs_step = nice_ticks(min(xs), max(xs) if max(xs) > min(xs) else min(xs) + 1.0)
    y0, y1, ys_step = nice_ticks(min(ys), max(ys) if max(ys) > min(ys) else min(ys) + 1.0)
    pw, ph = width - ML - MR, height - MT - MB

    def cx(x):
        return ML + (x - x0) / (x1 - x0) * pw

    def cy(y):
        return MT + (1.0 - (y - y0) / (y1 - y0)) * ph

    canvas.create_rectangle(ML, MT, ML + pw, MT + ph)
    for v in tick_values(x0, x1, xs_step):
        canvas.create_line(cx(v), MT + ph, cx(v), MT + ph + 4)
        canvas.create_text(cx(v), MT + ph + 14, text=tick_label(v, xs_step), font=("TkDefaultFont", 8))
    for v in tick_values(y0, y1, ys_step):
        canvas.create_line(ML - 4, cy(v), ML, cy(v))
        canvas.create_text(ML - 8, cy(v), text=tick_label(v, ys_step), anchor="e", font=("TkDefaultFont", 8))
    canvas.create_text(ML + pw / 2, height - 12, text=x_label)
    canvas.create_text(14, MT + ph / 2, text=y_label, angle=90)
    ly = MT + 10
    for name, colour, s in series:
        seg = [(cx(x), cy(y)) for x, y in s if y is not None]
        for (a, b), (c, d) in zip(seg, seg[1:]):
            canvas.create_line(a, b, c, d, fill=colour, width=2)
        for a, b in seg:
            canvas.create_oval(a - 2, b - 2, a + 2, b + 2, outline=colour)
        canvas.create_text(ML + pw - 8, ly, text=name, fill=colour, anchor="e")
        ly += 14
    for x, y, colour, text in markers:
        canvas.create_oval(cx(x) - 5, cy(y) - 5, cx(x) + 5, cy(y) + 5, outline=colour, width=2)
        canvas.create_text(cx(x) + 8, cy(y) - 8, text=text, fill=colour, anchor="w")


def _linspace(a: float, b: float, n: int) -> list[float]:
    return [a + (b - a) * i / (n - 1) for i in range(n)] if n > 1 else [a]


def open_sensitivity_dialog(app) -> None:
    """Studies → Sensitivity: vary one parameter, plot one output."""
    from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene
    from angelica.studies import sensitivity
    try:
        case = build_network_case_from_scene(app.scene)
    except Exception as exc:
        messagebox.showerror("Sensitivity", str(exc))
        return
    targets, outputs = study_targets(case), study_outputs(case)
    win = tk.Toplevel(app.root)
    win.title("Sensitivity study")
    frm = ttk.Frame(win, padding=10)
    frm.pack(fill="both", expand=True)
    t_var = tk.StringVar(value=targets[0][0] if targets else "")
    o_var = tk.StringVar(value=outputs[0][0] if outputs else "")
    a_var, b_var, n_var = tk.StringVar(), tk.StringVar(), tk.StringVar(value="8")
    ttk.Label(frm, text="Vary").grid(row=0, column=0, sticky="w")
    ttk.Combobox(frm, textvariable=t_var, values=[t[0] for t in targets], state="readonly", width=48).grid(
        row=0, column=1, columnspan=3, sticky="ew")
    for i, (lbl, var) in enumerate((("From", a_var), ("To", b_var), ("Points", n_var))):
        ttk.Label(frm, text=lbl).grid(row=1, column=2 * i if i < 2 else 4, sticky="w")
        ttk.Entry(frm, textvariable=var, width=12).grid(row=1, column=2 * i + 1 if i < 2 else 5)
    ttk.Label(frm, text="Plot").grid(row=2, column=0, sticky="w")
    ttk.Combobox(frm, textvariable=o_var, values=[o[0] for o in outputs], state="readonly", width=48).grid(
        row=2, column=1, columnspan=3, sticky="ew")
    canvas = tk.Canvas(frm, width=640, height=400, background="white")
    canvas.grid(row=4, column=0, columnspan=6, pady=(8, 0))
    status = tk.StringVar(value="SI units (the solver's).")
    ttk.Label(frm, textvariable=status).grid(row=5, column=0, columnspan=6, sticky="w")

    def run() -> None:
        try:
            target = dict(targets)[t_var.get()]
            output = dict(outputs)[o_var.get()]
            values = _linspace(float(a_var.get()), float(b_var.get()), max(int(n_var.get()), 2))
        except (KeyError, ValueError):
            messagebox.showerror("Sensitivity", "Choose a parameter and an output, and give From, To and Points.",
                                 parent=win)
            return
        status.set(f"Solving {len(values)} cases…")
        win.update_idletasks()
        try:
            points = sensitivity(case, lambda: build_solver_from_scene(app.scene), target, values)
        except ValueError as exc:
            messagebox.showerror("Sensitivity", str(exc), parent=win)
            status.set("")
            return
        data = [(p.values[0], output_value(p, output)) for p in points]
        bad = [f"{x:g}" for x, y in data if y is None]
        draw_plot(canvas, [(o_var.get(), "#1f77b4", data)], t_var.get(), o_var.get())
        status.set("All converged." if not bad else f"Not converged at: {', '.join(bad)} (gaps in the line).")

    ttk.Button(frm, text="Run", command=run).grid(row=3, column=5, sticky="e", pady=4)


def open_nodal_dialog(app) -> None:
    """Studies → Nodal analysis: inflow and outflow curves of a well."""
    from angelica.core.components import Well
    from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene
    from angelica.studies import nodal_analysis
    try:
        case = build_network_case_from_scene(app.scene)
    except Exception as exc:
        messagebox.showerror("Nodal analysis", str(exc))
        return
    wells = [c.component_id for c in case.components if isinstance(c, Well)]
    if not wells:
        messagebox.showinfo("Nodal analysis", "The network has no Well component.")
        return
    win = tk.Toplevel(app.root)
    win.title("Nodal analysis")
    frm = ttk.Frame(win, padding=10)
    frm.pack(fill="both", expand=True)
    w_var = tk.StringVar(value=wells[0])
    ttk.Label(frm, text="Well").grid(row=0, column=0, sticky="w")
    ttk.Combobox(frm, textvariable=w_var, values=wells, state="readonly", width=40).grid(row=0, column=1, sticky="w")
    canvas = tk.Canvas(frm, width=640, height=400, background="white")
    canvas.grid(row=1, column=0, columnspan=3, pady=(8, 0))
    status = tk.StringVar(value="Wellhead pressure (Pa) against the well's mass rate (kg/s).")
    ttk.Label(frm, textvariable=status, wraplength=620).grid(row=2, column=0, columnspan=3, sticky="w")

    def run() -> None:
        status.set("Solving the outflow curve…")
        win.update_idletasks()
        try:
            na = nodal_analysis(case, lambda: build_solver_from_scene(app.scene), w_var.get())
        except ValueError as exc:
            messagebox.showerror("Nodal analysis", str(exc), parent=win)
            return
        markers = []
        if na.operating_point is not None:
            m, p = na.operating_point
            markers.append((m, p, "#d62728", f"{m:.3g} kg/s, {p / 1e5:.4g} bar"))
        draw_plot(canvas, [("Inflow (IPR + tubing)", "#2ca02c", list(na.inflow)),
                           ("Outflow (network)", "#1f77b4", list(na.outflow))],
                  "Mass rate (kg/s)", "Wellhead pressure (Pa)", markers=markers)
        status.set(na.note or "Operating point marked in red (where the network solve puts the well).")

    ttk.Button(frm, text="Run", command=run).grid(row=0, column=2, sticky="e")
