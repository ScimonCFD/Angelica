from __future__ import annotations

import copy
import math
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import sv_ttk

from angelica.core.case import NetworkCase
from angelica.core.components import FITTING_PRESET_LIBRARY
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.io import export_solve_result_csv, export_solve_result_workbook

from .io import (
    build_network_case_from_scene,
    build_solver_from_scene,
    load_scene_and_results_from_file,
    load_scene_from_file,
    save_scene_to_file,
    scene_from_dict,
    scene_to_dict,
)
from .model import (
    CanvasLink,
    CanvasLinkComponent,
    CanvasNode,
    CanvasScene,
    DEFAULT_LIBRARY_MATERIAL,
    DEFAULT_PRESSURE_DROP_MODEL,
)


class ToolTip:
    def __init__(self, widget: tk.Widget, text: str) -> None:
        self.widget = widget
        self.text = text
        self.tip_window: tk.Toplevel | None = None

        self.widget.bind("<Enter>", self._show)
        self.widget.bind("<Leave>", self._hide)

    def _show(self, _event: tk.Event) -> None:
        if self.tip_window is not None:
            return

        x = self.widget.winfo_rootx() + self.widget.winfo_width() + 8
        y = self.widget.winfo_rooty() + 4

        self.tip_window = tk.Toplevel(self.widget)
        self.tip_window.wm_overrideredirect(True)
        self.tip_window.wm_geometry(f"+{x}+{y}")

        label = tk.Label(
            self.tip_window,
            text=self.text,
            background="#fff8dc",
            relief="solid",
            borderwidth=1,
            padx=6,
            pady=3,
        )
        label.pack()

    def _hide(self, _event: tk.Event) -> None:
        if self.tip_window is not None:
            self.tip_window.destroy()
            self.tip_window = None


class NetSimGui:
    _THEMES: dict[str, dict[str, str]] = {
        "light": {
            "canvas_bg":      "#fbfaf4",
            "canvas_hl":      "#b8b2a7",
            "node_source":    "#8ecae6",
            "node_sink":      "#f28482",
            "node_junction":  "#d9d9d9",
            "node_outline":   "#2e2a24",
            "node_text":      "#000000",
            "node_summary":   "#3d3a35",
            "link":           "#4f5d75",
            "drag_line":      "#6c757d",
            "plot_bg":        "#ffffff",
            "plot_axis":      "#333333",
            "plot_grid":      "#e6e6e6",
            "plot_text":      "#333333",
            "plot_muted":     "#555555",
            "plot_faint":     "#999999",
            "plot_faint2":    "#666666",
            "plot_laminar":     "#1d3557",
            "plot_turbulent":   "#c1121f",
            "plot_temperature": "#1a7a3c",
        },
        "dark": {
            "canvas_bg":      "#070c17",
            "canvas_hl":      "#2a3a4a",
            "node_source":    "#3a9fd4",
            "node_sink":      "#e8633a",
            "node_junction":  "#4a6a7a",
            "node_outline":   "#c8d4dc",
            "node_text":      "#e8e8e8",
            "node_summary":   "#9ab0bc",
            "link":           "#8aaab8",
            "drag_line":      "#8aaab8",
            "plot_bg":        "#0d1a24",
            "plot_axis":      "#c8d4dc",
            "plot_grid":      "#1e3040",
            "plot_text":      "#c8d4dc",
            "plot_muted":     "#7a9aaa",
            "plot_faint":     "#4a6a7a",
            "plot_faint2":    "#5a7a8a",
            "plot_laminar":     "#3a9fd4",
            "plot_turbulent":   "#e8633a",
            "plot_temperature": "#4dd476",
        },
    }

    FITTING_MODE_LIBRARY = {
        "manual": {"name": "Manual K"},
        "preset": {"name": "Preset from table"},
    }
    FITTING_PRESET_LIBRARY = FITTING_PRESET_LIBRARY
    MATERIAL_LIBRARY = {
        "water_liquid": {
            "definition_mode": "library",
            "name": "Water",
            "density_kg_per_m3": "998.25",
            "viscosity_pa_s": "0.001",
            "specific_heat_j_per_kg_k": "4182.0",
            "thermal_conductivity_w_per_m_k": "0.598",
        }
    }
    PRESSURE_DROP_MODEL_LIBRARY = {
        "colebrook_white": {
            "name": "Colebrook-White",
        },
        "hazen_williams": {
            "name": "Hazen-Williams",
        }
    }
    VELOCITY_LOOP_METHOD_LIBRARY = {
        "fixed_point": {
            "name": "Fixed-point",
        },
        "secant": {
            "name": "Secant",
        },
    }
    FRICTION_FACTOR_METHOD_LIBRARY = {
        "newton": {
            "name": "Newton",
        },
    }
    COLEBROOK_FRICTION_STRATEGY_LIBRARY = {
        "transformed": {
            "name": "Protected (log f)",
        },
        "direct": {
            "name": "Direct f",
        },
    }
    METRIC_OPTIONS = (
        ("Max abs ΔP correction (Pa)", "pressure_correction_abs_pa"),
        ("Max rel mass imbalance (−)", "max_nodal_mass_imbalance_rel"),
    )
    UNIT_SYSTEMS: dict[str, dict] = {
        "si": {
            "name": "SI  (m · Pa · kg/s)",
            "quantities": {
                "pressure":      ("Pa",   1.0),
                "flow":          ("kg/s", 1.0),
                "length":        ("m",    1.0),
                "diameter":      ("m",    1.0),
                "roughness":     ("m",    1.0),
                "height_change": ("m",    1.0),
            },
        },
        "us": {
            "name": "US Customary  (ft · psi · lb/s)",
            "quantities": {
                "pressure":      ("psi",  1.0 / 6894.757),
                "flow":          ("lb/s", 1.0 / 0.453592),
                "length":        ("ft",   1.0 / 0.3048),
                "diameter":      ("in",   1.0 / 0.0254),
                "roughness":     ("in",   1.0 / 0.0254),
                "height_change": ("ft",   1.0 / 0.3048),
            },
        },
    }
    _FIELD_QUANTITY: dict[str, str] = {
        "length_m":        "length",
        "diameter_m":      "diameter",
        "roughness_m":     "roughness",
        "height_change_m": "height_change",
    }

    def __init__(self) -> None:
        self.scene = CanvasScene()
        self.drag_source_node_id: int | None = None
        self.drag_line_id: int | None = None
        self.moving_node_id: int | None = None
        self.selected_node_id: int | None = None
        self._last_dir: str = self._default_open_dir()
        self.middle_pan_anchor: tuple[float, float] | None = None
        self.view_offset_x = 0.0
        self.view_offset_y = 0.0
        self.view_scale = 1.0
        self.latest_result = None
        self.latest_boundary_results: dict[int, dict[str, float]] = {}
        self.current_file_path: str | None = None
        self._undo_stack: list[dict] = []
        self._redo_stack: list[dict] = []
        self._move_pre_snapshot: dict | None = None
        self.convergence_window: tk.Toplevel | None = None
        self.convergence_canvas: tk.Canvas | None = None
        self.temperature_canvas: tk.Canvas | None = None
        self.temperature_history: list[float] = []
        self.density_history: list[float] = []
        self.outer_turbulent_final_metrics: list = []
        self._dark = False
        self._unit_system_key = "si"
        self.root = tk.Tk()
        sv_ttk.set_theme("light")
        self.root.title("Angelica")
        self.root.geometry("1100x700")
        self.root.minsize(900, 600)
        self._set_window_icon()

        self.metric_label_to_name = {label: name for label, name in self.METRIC_OPTIONS}
        self.metric_name_to_label = {name: label for label, name in self.METRIC_OPTIONS}
        self.convergence_metric_var = tk.StringVar(
            master=self.root,
            value=self.metric_name_to_label["pressure_correction_abs_pa"],
        )
        self.show_hydraulic_detail_var = tk.BooleanVar(master=self.root, value=False)
        self.convergence_history = {"laminar": [], "turbulent": []}
        self.status_var = tk.StringVar(value="Select a node type from the palette.")
        self.tool_var = tk.StringVar(value="No tool selected")
        self.material_summary_var = tk.StringVar(value=self._material_summary_text())
        self.pressure_drop_summary_var = tk.StringVar(value=self._pressure_drop_summary_text())
        self.numerics_summary_var = tk.StringVar(value=self._numerics_summary_text())

        self._build_menu()
        self._build_layout()

    @property
    def _t(self) -> dict[str, str]:
        return self._THEMES["dark" if self._dark else "light"]

    def _toggle_theme(self) -> None:
        self._dark = not self._dark
        sv_ttk.set_theme("dark" if self._dark else "light")
        self.canvas.configure(
            background=self._t["canvas_bg"],
            highlightbackground=self._t["canvas_hl"],
        )
        if self.convergence_canvas is not None:
            self.convergence_canvas.configure(
                background=self._t["plot_bg"],
                highlightbackground=self._t["canvas_hl"],
            )
        self._redraw_scene()
        self._redraw_convergence_plot()

    def _build_menu(self) -> None:
        menu_bar = tk.Menu(self.root)
        file_menu = tk.Menu(menu_bar, tearoff=False)
        file_menu.add_command(label="New", command=self._new_scene)
        file_menu.add_command(label="Open…", command=self._open_scene)
        file_menu.add_separator()
        file_menu.add_command(label="Save", accelerator="Ctrl+S", command=self._save_scene)
        file_menu.add_command(label="Save As…", command=self._save_scene_as)
        file_menu.add_separator()
        file_menu.add_command(label="Export Results Report (Excel)", command=self._export_results_report)
        file_menu.add_command(label="Export Results Report (CSV)", command=self._export_results_report_csv)
        file_menu.add_separator()
        file_menu.add_command(label="Close", command=self.root.destroy)
        menu_bar.add_cascade(label="File", menu=file_menu)
        self.root.bind("<Control-s>", lambda _event: self._save_scene())

        edit_menu = tk.Menu(menu_bar, tearoff=False)
        edit_menu.add_command(label="Undo", accelerator="Ctrl+Z", command=self._undo)
        edit_menu.add_command(label="Redo", accelerator="Ctrl+Y", command=self._redo)
        menu_bar.add_cascade(label="Edit", menu=edit_menu)
        self.root.bind("<Control-z>", lambda _event: self._undo())
        self.root.bind("<Control-y>", lambda _event: self._redo())

        material_menu = tk.Menu(menu_bar, tearoff=False)
        material_menu.add_command(label="Define Material", command=self._open_material_dialog)
        menu_bar.add_cascade(label="Material", menu=material_menu)

        physics_menu = tk.Menu(menu_bar, tearoff=False)
        physics_menu.add_command(
            label="Case Type…",
            command=self._open_case_type_dialog,
        )
        physics_menu.add_separator()
        physics_menu.add_command(
            label="Define Pressure-Drop Model",
            command=self._open_pressure_drop_model_dialog,
        )
        menu_bar.add_cascade(label="Physics", menu=physics_menu)

        numerics_menu = tk.Menu(menu_bar, tearoff=False)
        numerics_menu.add_command(
            label="Define Numerics",
            command=self._open_numerics_dialog,
        )
        menu_bar.add_cascade(label="Numerics", menu=numerics_menu)

        view_menu = tk.Menu(menu_bar, tearoff=False)
        view_menu.add_command(label="Toggle Dark / Light Theme", command=self._toggle_theme)
        view_menu.add_command(label="Convergence Window", command=self._open_convergence_window)
        menu_bar.add_cascade(label="View", menu=view_menu)

        settings_menu = tk.Menu(menu_bar, tearoff=False)
        settings_menu.add_command(label="Unit System…", command=self._open_unit_system_dialog)
        menu_bar.add_cascade(label="Settings", menu=settings_menu)

        self.root.config(menu=menu_bar)

    def _build_layout(self) -> None:
        container = ttk.Frame(self.root, padding=8)
        container.pack(fill="both", expand=True)

        palette = ttk.Frame(container, padding=(8, 8, 12, 8))
        palette.pack(side="left", fill="y")

        palette_title = ttk.Label(palette, text="Node Palette")
        palette_title.pack(anchor="w", pady=(0, 8))

        ttk.Button(
            palette,
            text="▶ Run",
            command=self._run_simulation,
            width=12,
        ).pack(anchor="w", pady=(0, 10))

        source_button = ttk.Button(
            palette,
            text="▲",
            command=lambda: self._select_tool("source"),
            width=4,
        )
        source_button.pack(anchor="w", pady=4)
        ToolTip(source_button, "Source")

        sink_button = ttk.Button(
            palette,
            text="▼",
            command=lambda: self._select_tool("sink"),
            width=4,
        )
        sink_button.pack(anchor="w", pady=4)
        ToolTip(sink_button, "Sink")

        junction_button = ttk.Button(
            palette,
            text="○",
            command=lambda: self._select_tool("junction"),
            width=4,
        )
        junction_button.pack(anchor="w", pady=4)
        ToolTip(junction_button, "Junction")

        ttk.Separator(palette, orient="horizontal").pack(fill="x", pady=10)

        ttk.Label(palette, text="Active Tool").pack(anchor="w")
        ttk.Label(
            palette,
            textvariable=self.tool_var,
            width=20,
            relief="groove",
            padding=6,
        ).pack(anchor="w", pady=(4, 0))

        ttk.Separator(palette, orient="horizontal").pack(fill="x", pady=10)

        ttk.Label(palette, text="Material").pack(anchor="w")
        ttk.Label(
            palette,
            textvariable=self.material_summary_var,
            width=26,
            relief="groove",
            padding=6,
            justify="left",
        ).pack(anchor="w", pady=(4, 0), fill="x")

        ttk.Label(palette, text="Pipe Model").pack(anchor="w", pady=(10, 0))
        ttk.Label(
            palette,
            textvariable=self.pressure_drop_summary_var,
            width=26,
            relief="groove",
            padding=6,
            justify="left",
        ).pack(anchor="w", pady=(4, 0), fill="x")

        ttk.Label(palette, text="Numerics").pack(anchor="w", pady=(10, 0))
        ttk.Label(
            palette,
            textvariable=self.numerics_summary_var,
            width=26,
            relief="groove",
            padding=6,
            justify="left",
        ).pack(anchor="w", pady=(4, 0), fill="x")

        canvas_frame = ttk.Frame(container)
        canvas_frame.pack(side="left", fill="both", expand=True)

        self.canvas = tk.Canvas(
            canvas_frame,
            background=self._t["canvas_bg"],
            highlightthickness=1,
            highlightbackground=self._t["canvas_hl"],
        )
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<ButtonPress-1>", self._on_canvas_press)
        self.canvas.bind("<B1-Motion>", self._on_canvas_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_canvas_release)
        self.canvas.bind("<Shift-ButtonPress-1>", self._on_canvas_shift_press)
        self.canvas.bind("<Shift-B1-Motion>", self._on_canvas_shift_drag)
        self.canvas.bind("<Shift-ButtonRelease-1>", self._on_canvas_shift_release)
        self.canvas.bind("<ButtonPress-3>", self._on_canvas_right_press)
        self.canvas.bind("<B3-Motion>", self._on_canvas_right_drag)
        self.canvas.bind("<ButtonRelease-3>", self._on_canvas_right_release)
        self.canvas.bind("<ButtonPress-2>", self._on_canvas_middle_press)
        self.canvas.bind("<B2-Motion>", self._on_canvas_middle_drag)
        self.canvas.bind("<ButtonRelease-2>", self._on_canvas_middle_release)
        self.canvas.bind("<Double-Button-1>", self._on_canvas_double_click)
        self.canvas.bind("<MouseWheel>", self._on_canvas_scroll)
        self.canvas.bind("<Button-4>", self._on_canvas_scroll)
        self.canvas.bind("<Button-5>", self._on_canvas_scroll)
        self.root.bind("<Delete>", self._on_delete_key)
        self.root.bind("<BackSpace>", self._on_delete_key)

        status = ttk.Label(
            self.root,
            textvariable=self.status_var,
            relief="sunken",
            anchor="w",
            padding=(8, 4),
        )
        status.pack(fill="x", side="bottom")

    def _select_tool(self, tool: str) -> None:
        self.scene.set_active_tool(tool)
        self.tool_var.set(tool.capitalize())
        self.status_var.set(f"{tool.capitalize()} selected. Click on the canvas to place it.")

    def _new_scene(self) -> None:
        self.scene.clear()
        self.canvas.delete("all")
        self.drag_source_node_id = None
        self.drag_line_id = None
        self.moving_node_id = None
        self.selected_node_id = None
        self.middle_pan_anchor = None
        self.view_offset_x = 0.0
        self.view_offset_y = 0.0
        self.view_scale = 1.0
        self.latest_result = None
        self.latest_boundary_results = {}
        self.temperature_history = []
        self.density_history = []
        self.outer_turbulent_final_metrics = []
        self.current_file_path = None
        self._undo_stack.clear()
        self._redo_stack.clear()
        self.tool_var.set("No tool selected")
        self._update_title()
        self._refresh_global_summaries()
        self.status_var.set("New scene created. Select a node type from the palette.")

    def _open_scene(self) -> None:
        file_path = filedialog.askopenfilename(
            title="Open Angelica GUI case",
            initialdir=self._last_dir,
            filetypes=(
                ("Angelica GUI case", "*.gui.json"),
                ("JSON files", "*.json"),
                ("All files", "*.*"),
            ),
        )
        if not file_path:
            return

        self._last_dir = os.path.dirname(file_path)

        try:
            scene, cached_results = load_scene_and_results_from_file(file_path)
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Open failed", f"Could not open case:\n{exc}")
            return

        self.scene = scene
        self.canvas.delete("all")
        self.drag_source_node_id = None
        self.drag_line_id = None
        self.moving_node_id = None
        self.selected_node_id = None
        self.middle_pan_anchor = None
        self.view_offset_x = 0.0
        self.view_offset_y = 0.0
        self.view_scale = 1.0
        self.latest_result = None
        self.current_file_path = file_path
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._update_title()

        if cached_results is not None:
            self.latest_boundary_results = cached_results["boundary_results"]
            self.convergence_history = cached_results["convergence_history"]
            self.temperature_history = []
            self.density_history = []
            self.outer_turbulent_final_metrics = []
            if cached_results["component_flows"]:
                node_pressures = {
                    node_id: entry["pressure_pa"]
                    for node_id, entry in cached_results["boundary_results"].items()
                }
                self.latest_result = SolveResult(
                    case_name="GUI scene",
                    converged=cached_results["converged"],
                    node_pressures_pa=node_pressures,
                    component_flows=cached_results["component_flows"],
                    laminar_history=[],
                    laminar_metrics=cached_results["convergence_history"].get("laminar", []),
                    turbulent_history=[],
                    turbulent_metrics=cached_results["convergence_history"].get("turbulent", []),
                )
            status_suffix = " — results restored" if cached_results["converged"] else " — unconverged results restored"
        else:
            self.latest_boundary_results = {}
            self.convergence_history = {"laminar": [], "turbulent": []}
            self.temperature_history = []
            self.density_history = []
            self.outer_turbulent_final_metrics = []
            status_suffix = ""

        self.tool_var.set("No tool selected")
        self._refresh_global_summaries()
        self._redraw_scene()
        self.status_var.set(f"Opened: {os.path.basename(file_path)}{status_suffix}")

    @staticmethod
    def _is_tutorial_file(file_path: str) -> bool:
        """Return True if file_path lives inside any directory named 'tutorials'."""
        parts = os.path.normcase(os.path.abspath(file_path)).split(os.sep)
        return "tutorials" in parts

    def _update_title(self) -> None:
        if self.current_file_path:
            name = os.path.basename(self.current_file_path)
            suffix = " [tutorial]" if self._is_tutorial_file(self.current_file_path) else ""
            self.root.title(f"Angelica — {name}{suffix}")
        else:
            self.root.title("Angelica")

    def _save_scene(self) -> None:
        if self.current_file_path and not self._is_tutorial_file(self.current_file_path):
            self._do_save(self.current_file_path)
        else:
            if self.current_file_path:
                messagebox.showinfo(
                    "Tutorial — read-only",
                    "Tutorial files cannot be overwritten.\n\n"
                    "Choose a different location to save your work.",
                )
            self._save_scene_as()

    def _save_scene_as(self) -> None:
        while True:
            file_path = filedialog.asksaveasfilename(
                title="Save Angelica GUI case",
                initialdir=self._last_dir,
                defaultextension=".gui.json",
                filetypes=(
                    ("Angelica GUI case", "*.gui.json"),
                    ("JSON files", "*.json"),
                    ("All files", "*.*"),
                ),
            )
            if not file_path:
                return
            if not self._is_tutorial_file(file_path):
                break
            messagebox.showwarning(
                "Tutorial folder — read-only",
                "The tutorials folder is protected and cannot be used as a save location.\n\n"
                "Please choose a different folder.",
            )
        self._last_dir = os.path.dirname(file_path)
        self.current_file_path = file_path
        self._update_title()
        self._do_save(file_path)

    def _do_save(self, file_path: str) -> None:
        converged = self.latest_result.converged if self.latest_result is not None else False
        component_flows = (
            list(self.latest_result.component_flows) if self.latest_result is not None else None
        )
        try:
            save_scene_to_file(
                self.scene,
                file_path,
                boundary_results=self.latest_boundary_results or None,
                convergence_history=self.convergence_history if self.latest_boundary_results else None,
                converged=converged,
                component_flows=component_flows,
            )
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Save failed", f"Could not save case:\n{exc}")
            return
        self.status_var.set(f"Saved: {os.path.basename(file_path)}")

    def _make_snapshot(self) -> dict:
        return copy.deepcopy(scene_to_dict(self.scene))

    def _push_undo(self) -> None:
        self._undo_stack.append(self._make_snapshot())
        self._redo_stack.clear()

    def _commit_move_if_changed(self) -> None:
        if self._move_pre_snapshot is None:
            return
        snap = self._move_pre_snapshot
        self._move_pre_snapshot = None
        current_nodes = scene_to_dict(self.scene)["nodes"]
        if current_nodes != snap["nodes"]:
            self._undo_stack.append(snap)
            self._redo_stack.clear()

    def _restore_snapshot(self, snap: dict) -> None:
        self.scene = scene_from_dict(snap)
        self.selected_node_id = None
        self.moving_node_id = None
        self.latest_boundary_results = {}
        self.latest_result = None
        self._redraw_scene()

    def _undo(self) -> None:
        if not self._undo_stack:
            self.status_var.set("Nothing to undo.")
            return
        self._redo_stack.append(self._make_snapshot())
        self._restore_snapshot(self._undo_stack.pop())
        self.status_var.set("Undo.")

    def _redo(self) -> None:
        if not self._redo_stack:
            self.status_var.set("Nothing to redo.")
            return
        self._undo_stack.append(self._make_snapshot())
        self._restore_snapshot(self._redo_stack.pop())
        self.status_var.set("Redo.")

    def _open_material_dialog(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title("Define Material")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=12)
        frame.pack(fill="both", expand=True)

        material = dict(self.scene.material)
        definition_mode = material.get("definition_mode", "library" if material.get("library_key") else "custom")

        library_var = tk.StringVar(
            master=dialog,
            value=material.get("library_key", DEFAULT_LIBRARY_MATERIAL["library_key"]),
        )
        mode_var = tk.StringVar(master=dialog, value=definition_mode)
        name_var = tk.StringVar(master=dialog, value=material.get("name", ""))
        density_var = tk.StringVar(master=dialog, value=material.get("density_kg_per_m3", ""))
        viscosity_var = tk.StringVar(master=dialog, value=material.get("viscosity_pa_s", ""))
        api_var = tk.StringVar(master=dialog, value=material.get("api_gravity", "30.0"))
        temperature_var = tk.StringVar(master=dialog, value=material.get("temperature_c", "60.0"))
        cp_var = tk.StringVar(master=dialog, value=material.get("specific_heat_j_per_kg_k", ""))
        k_var = tk.StringVar(master=dialog, value=material.get("thermal_conductivity_w_per_m_k", ""))
        mw_var = tk.StringVar(master=dialog, value=material.get("molecular_weight_kg_per_mol", ""))
        tc_var = tk.StringVar(master=dialog, value=material.get("critical_temperature_k", ""))
        pc_var = tk.StringVar(master=dialog, value=material.get("critical_pressure_pa", ""))
        omega_var = tk.StringVar(master=dialog, value=material.get("acentric_factor", ""))

        ttk.Label(frame, text="Definition").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        mode_row = ttk.Frame(frame)
        mode_row.grid(row=0, column=1, sticky="w", pady=4)
        ttk.Radiobutton(mode_row, text="Library", variable=mode_var, value="library").pack(side="left")
        ttk.Radiobutton(mode_row, text="Custom", variable=mode_var, value="custom").pack(side="left", padx=(8, 0))
        ttk.Radiobutton(mode_row, text="Crude oil", variable=mode_var, value="crude_oil").pack(side="left", padx=(8, 0))
        ttk.Radiobutton(mode_row, text="Gas (ideal)", variable=mode_var, value="gas").pack(side="left", padx=(8, 0))
        ttk.Radiobutton(mode_row, text="Gas (PR)", variable=mode_var, value="gas_pr").pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Material Library").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        library_box = ttk.Combobox(
            frame,
            textvariable=library_var,
            state="readonly",
            values=("water_liquid",),
            width=24,
        )
        library_box.grid(row=1, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Name").grid(row=2, column=0, sticky="w", padx=(0, 8), pady=4)
        name_entry = ttk.Entry(frame, textvariable=name_var, width=26)
        name_entry.grid(row=2, column=1, sticky="ew", pady=4)

        density_label = ttk.Label(frame, text="Density (kg/m³)")
        density_label.grid(row=3, column=0, sticky="w", padx=(0, 8), pady=4)
        density_entry = ttk.Entry(frame, textvariable=density_var, width=26)
        density_entry.grid(row=3, column=1, sticky="ew", pady=4)

        mw_label = ttk.Label(frame, text="Molecular Weight M (kg/mol)")
        mw_label.grid(row=3, column=0, sticky="w", padx=(0, 8), pady=4)
        mw_entry = ttk.Entry(frame, textvariable=mw_var, width=26)
        mw_entry.grid(row=3, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Viscosity (Pa·s)").grid(row=4, column=0, sticky="w", padx=(0, 8), pady=4)
        viscosity_entry = ttk.Entry(frame, textvariable=viscosity_var, width=26)
        viscosity_entry.grid(row=4, column=1, sticky="ew", pady=4)

        api_label = ttk.Label(frame, text="API Gravity (°)")
        api_label.grid(row=5, column=0, sticky="w", padx=(0, 8), pady=4)
        api_entry = ttk.Entry(frame, textvariable=api_var, width=26)
        api_entry.grid(row=5, column=1, sticky="ew", pady=4)

        tc_label = ttk.Label(frame, text="Critical Temperature Tc (K)")
        tc_label.grid(row=5, column=0, sticky="w", padx=(0, 8), pady=4)
        tc_entry = ttk.Entry(frame, textvariable=tc_var, width=26)
        tc_entry.grid(row=5, column=1, sticky="ew", pady=4)

        temperature_label = ttk.Label(frame, text="Temperature (°C)")
        temperature_label.grid(row=6, column=0, sticky="w", padx=(0, 8), pady=4)
        temperature_entry = ttk.Entry(frame, textvariable=temperature_var, width=26)
        temperature_entry.grid(row=6, column=1, sticky="ew", pady=4)

        pc_label = ttk.Label(frame, text="Critical Pressure Pc (Pa)")
        pc_label.grid(row=6, column=0, sticky="w", padx=(0, 8), pady=4)
        pc_entry = ttk.Entry(frame, textvariable=pc_var, width=26)
        pc_entry.grid(row=6, column=1, sticky="ew", pady=4)

        omega_label = ttk.Label(frame, text="Acentric Factor ω")
        omega_label.grid(row=7, column=0, sticky="w", padx=(0, 8), pady=4)
        omega_entry = ttk.Entry(frame, textvariable=omega_var, width=26)
        omega_entry.grid(row=7, column=1, sticky="ew", pady=4)

        cp_label = ttk.Label(frame, text="Specific Heat cp (J/kg·K)")
        cp_label.grid(row=8, column=0, sticky="w", padx=(0, 8), pady=4)
        cp_entry = ttk.Entry(frame, textvariable=cp_var, width=26)
        cp_entry.grid(row=8, column=1, sticky="ew", pady=4)

        k_label = ttk.Label(frame, text="Thermal Conductivity k (W/m·K)")
        k_label.grid(row=9, column=0, sticky="w", padx=(0, 8), pady=4)
        k_entry = ttk.Entry(frame, textvariable=k_var, width=26)
        k_entry.grid(row=9, column=1, sticky="ew", pady=4)

        frame.columnconfigure(1, weight=1)

        def apply_library_selection(_event: tk.Event | None = None) -> None:
            preset = self.MATERIAL_LIBRARY[library_var.get()]
            name_var.set(preset["name"])
            density_var.set(preset["density_kg_per_m3"])
            viscosity_var.set(preset["viscosity_pa_s"])
            cp_var.set(preset.get("specific_heat_j_per_kg_k", ""))
            k_var.set(preset.get("thermal_conductivity_w_per_m_k", ""))

        is_non_isothermal = self.scene.physics_mode == "non_isothermal"
        is_compressible = self.scene.physics_mode == "compressible"

        def sync_mode_state(*_args: object) -> None:
            mode = mode_var.get()
            is_library = mode == "library"
            is_crude_oil = mode == "crude_oil"
            is_gas = mode == "gas"
            is_gas_pr = mode == "gas_pr"
            is_any_gas = is_gas or is_gas_pr
            library_box.configure(state="readonly" if is_library else "disabled")
            std_state = "normal" if mode == "custom" else "disabled"
            name_entry.configure(state=std_state if not is_any_gas else "normal")
            viscosity_entry.configure(state=std_state if not is_any_gas else "normal")
            # density vs molecular weight
            if is_any_gas:
                density_label.grid_remove()
                density_entry.grid_remove()
                mw_label.grid()
                mw_entry.grid()
                mw_entry.configure(state="normal")
            else:
                mw_label.grid_remove()
                mw_entry.grid_remove()
                density_label.grid()
                density_entry.grid()
                density_entry.configure(state=std_state)
            if is_crude_oil:
                api_label.grid()
                api_entry.grid()
                temperature_label.grid()
                temperature_entry.grid()
            else:
                api_label.grid_remove()
                api_entry.grid_remove()
                temperature_label.grid_remove()
                temperature_entry.grid_remove()
            if is_gas_pr:
                tc_label.grid()
                tc_entry.grid()
                tc_entry.configure(state="normal")
                pc_label.grid()
                pc_entry.grid()
                pc_entry.configure(state="normal")
                omega_label.grid()
                omega_entry.grid()
                omega_entry.configure(state="normal")
            else:
                tc_label.grid_remove()
                tc_entry.grid_remove()
                pc_label.grid_remove()
                pc_entry.grid_remove()
                omega_label.grid_remove()
                omega_entry.grid_remove()
            if is_non_isothermal or is_any_gas:
                cp_label.grid()
                cp_entry.grid()
                k_label.grid()
                k_entry.grid()
                cp_entry.configure(state="normal" if is_any_gas else std_state)
                k_entry.configure(state="normal" if is_any_gas else std_state)
            else:
                cp_label.grid_remove()
                cp_entry.grid_remove()
                k_label.grid_remove()
                k_entry.grid_remove()
            if is_library and library_var.get():
                apply_library_selection()

        library_box.bind("<<ComboboxSelected>>", apply_library_selection)
        mode_var.trace_add("write", sync_mode_state)
        sync_mode_state()

        button_row = ttk.Frame(frame)
        button_row.grid(row=10, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right")
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_material_definition(
                dialog,
                mode_var,
                library_var,
                name_var,
                density_var,
                viscosity_var,
                api_var,
                temperature_var,
                cp_var,
                k_var,
                mw_var,
                tc_var,
                pc_var,
                omega_var,
            ),
        ).pack(side="right", padx=(0, 8))

        dialog.update_idletasks()
        dialog.wait_visibility()
        dialog.grab_set()
        dialog.focus_set()
        name_entry.focus_set()

    def _open_pressure_drop_model_dialog(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title("Define Pressure-Drop Model")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=12)
        frame.pack(fill="both", expand=True)

        model = dict(DEFAULT_PRESSURE_DROP_MODEL)
        model.update(self.scene.pressure_drop_model)

        library_var = tk.StringVar(
            master=dialog,
            value=model.get("library_key", DEFAULT_PRESSURE_DROP_MODEL["library_key"]),
        )

        ttk.Label(frame, text="Pipe Model Library").grid(
            row=0, column=0, sticky="w", padx=(0, 8), pady=4
        )
        library_box = ttk.Combobox(
            frame,
            textvariable=library_var,
            state="readonly",
            values=tuple(self.PRESSURE_DROP_MODEL_LIBRARY.keys()),
            width=24,
        )
        library_box.grid(row=0, column=1, sticky="ew", pady=4)

        selected_name_var = tk.StringVar(
            master=dialog,
            value=self.PRESSURE_DROP_MODEL_LIBRARY[library_var.get()]["name"],
        )
        ttk.Label(frame, text="Selected Model").grid(
            row=1, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Label(
            frame,
            textvariable=selected_name_var,
            relief="groove",
            padding=6,
            width=24,
        ).grid(row=1, column=1, sticky="ew", pady=4)

        def apply_model_selection(_event: tk.Event | None = None) -> None:
            selected_name_var.set(
                self.PRESSURE_DROP_MODEL_LIBRARY[library_var.get()]["name"]
            )

        library_box.bind("<<ComboboxSelected>>", apply_model_selection)

        button_row = ttk.Frame(frame)
        button_row.grid(row=2, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right")
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_pressure_drop_model_definition(dialog, library_var),
        ).pack(side="right", padx=(0, 8))

        frame.columnconfigure(1, weight=1)
        dialog.update_idletasks()
        dialog.wait_visibility()
        dialog.grab_set()
        dialog.focus_set()
        library_box.focus_set()

    def _open_numerics_dialog(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title("Define Numerics")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=12)
        frame.pack(fill="both", expand=True)

        _li = self.scene.solver_settings.get("laminar_iterations")
        current_laminar_iterations = "auto" if _li is None else str(_li)
        current_turbulent_iterations = str(
            self.scene.solver_settings.get("turbulent_iterations", "60")
        )
        current_alpha = str(self.scene.solver_settings.get("pressure_relaxation", "1.0"))
        current_friction_max_iterations = str(
            self.scene.solver_settings.get("friction_factor_max_iterations", "50")
        )
        current_colebrook_strategy = str(
            self.scene.solver_settings.get("colebrook_friction_strategy", "transformed")
        )
        current_velocity_method = str(
            self.scene.solver_settings.get("velocity_loop_method", "fixed_point")
        )
        current_velocity_max_iterations = str(
            self.scene.solver_settings.get("velocity_loop_max_iterations", "50")
        )
        current_friction_method = str(
            self.scene.solver_settings.get("friction_factor_method", "newton")
        )
        current_colebrook_tol = self._fmt_sci(self.scene.solver_settings.get("colebrook_residual_tolerance", 1e-4))
        current_velocity_loop_tol = self._fmt_sci(self.scene.solver_settings.get("velocity_loop_tolerance", 1e-4))
        current_dp_tol = self._fmt_sci(self.scene.solver_settings.get("pressure_correction_abs_tolerance_pa", 1e-3))
        current_continuity_tol = self._fmt_sci(self.scene.solver_settings.get("nodal_mass_imbalance_rel_tolerance", 1e-3))

        # Non-isothermal energy settings
        current_convection_scheme = str(self.scene.solver_settings.get("convection_scheme", "upwind"))
        current_max_temp_iter = str(self.scene.solver_settings.get("max_temperature_iterations", "50"))
        current_temp_tol = self._fmt_sci(self.scene.solver_settings.get("temperature_tolerance_k", 0.01))
        current_temp_relax = str(self.scene.solver_settings.get("temperature_relaxation", "1.0"))

        convection_scheme_var = tk.StringVar(master=dialog, value=current_convection_scheme)
        max_temp_iter_var = tk.StringVar(master=dialog, value=current_max_temp_iter)
        temp_tol_var = tk.StringVar(master=dialog, value=current_temp_tol)
        temp_relax_var = tk.StringVar(master=dialog, value=current_temp_relax)

        alpha_var = tk.StringVar(master=dialog, value=current_alpha)
        friction_max_iterations_var = tk.StringVar(
            master=dialog,
            value=current_friction_max_iterations,
        )
        colebrook_strategy_var = tk.StringVar(
            master=dialog,
            value=current_colebrook_strategy,
        )
        velocity_method_var = tk.StringVar(master=dialog, value=current_velocity_method)
        velocity_max_iterations_var = tk.StringVar(
            master=dialog,
            value=current_velocity_max_iterations,
        )
        friction_method_var = tk.StringVar(master=dialog, value=current_friction_method)
        colebrook_tol_var = tk.StringVar(master=dialog, value=current_colebrook_tol)
        velocity_loop_tol_var = tk.StringVar(master=dialog, value=current_velocity_loop_tol)
        dp_tol_var = tk.StringVar(master=dialog, value=current_dp_tol)
        continuity_tol_var = tk.StringVar(master=dialog, value=current_continuity_tol)
        laminar_iterations_var = tk.StringVar(
            master=dialog,
            value=current_laminar_iterations,
        )
        turbulent_iterations_var = tk.StringVar(
            master=dialog,
            value=current_turbulent_iterations,
        )

        ttk.Label(frame, text="Laminar Iterations").grid(
            row=0, column=0, sticky="w", padx=(0, 8), pady=4
        )
        laminar_iterations_entry = ttk.Entry(
            frame,
            textvariable=laminar_iterations_var,
            width=26,
        )
        laminar_iterations_entry.grid(row=0, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Turbulent Iterations").grid(
            row=1, column=0, sticky="w", padx=(0, 8), pady=4
        )
        turbulent_iterations_entry = ttk.Entry(
            frame,
            textvariable=turbulent_iterations_var,
            width=26,
        )
        turbulent_iterations_entry.grid(row=1, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Pressure Relaxation").grid(
            row=2, column=0, sticky="w", padx=(0, 8), pady=4
        )
        alpha_entry = ttk.Entry(frame, textvariable=alpha_var, width=26)
        alpha_entry.grid(row=2, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Velocity Loop").grid(
            row=3, column=0, sticky="w", padx=(0, 8), pady=4
        )
        velocity_method_box = ttk.Combobox(
            frame,
            textvariable=velocity_method_var,
            state="readonly",
            values=tuple(self.VELOCITY_LOOP_METHOD_LIBRARY.keys()),
            width=24,
        )
        velocity_method_box.grid(row=3, column=1, sticky="ew", pady=4)

        velocity_name_var = tk.StringVar(
            master=dialog,
            value=self.VELOCITY_LOOP_METHOD_LIBRARY[velocity_method_var.get()]["name"],
        )
        ttk.Label(frame, text="Selected Velocity Loop").grid(
            row=4, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Label(
            frame,
            textvariable=velocity_name_var,
            relief="groove",
            padding=6,
            width=24,
        ).grid(row=4, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Velocity Max Iterations").grid(
            row=5, column=0, sticky="w", padx=(0, 8), pady=4
        )
        velocity_max_iterations_entry = ttk.Entry(
            frame,
            textvariable=velocity_max_iterations_var,
            width=26,
        )
        velocity_max_iterations_entry.grid(row=5, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Colebrook Strategy").grid(
            row=6, column=0, sticky="w", padx=(0, 8), pady=4
        )
        colebrook_strategy_box = ttk.Combobox(
            frame,
            textvariable=colebrook_strategy_var,
            state="readonly",
            values=tuple(self.COLEBROOK_FRICTION_STRATEGY_LIBRARY.keys()),
            width=24,
        )
        colebrook_strategy_box.grid(row=6, column=1, sticky="ew", pady=4)

        colebrook_strategy_name_var = tk.StringVar(
            master=dialog,
            value=self.COLEBROOK_FRICTION_STRATEGY_LIBRARY[colebrook_strategy_var.get()]["name"],
        )
        ttk.Label(frame, text="Selected Colebrook Strategy").grid(
            row=7, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Label(
            frame,
            textvariable=colebrook_strategy_name_var,
            relief="groove",
            padding=6,
            width=24,
        ).grid(row=7, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Friction Factor").grid(
            row=8, column=0, sticky="w", padx=(0, 8), pady=4
        )
        friction_method_box = ttk.Combobox(
            frame,
            textvariable=friction_method_var,
            state="readonly",
            values=tuple(self.FRICTION_FACTOR_METHOD_LIBRARY.keys()),
            width=24,
        )
        friction_method_box.grid(row=8, column=1, sticky="ew", pady=4)

        friction_name_var = tk.StringVar(
            master=dialog,
            value=self.FRICTION_FACTOR_METHOD_LIBRARY[friction_method_var.get()]["name"],
        )
        ttk.Label(frame, text="Selected Friction Method").grid(
            row=9, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Label(
            frame,
            textvariable=friction_name_var,
            relief="groove",
            padding=6,
            width=24,
        ).grid(row=9, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Friction Max Iterations").grid(
            row=10, column=0, sticky="w", padx=(0, 8), pady=4
        )
        friction_max_iterations_entry = ttk.Entry(
            frame,
            textvariable=friction_max_iterations_var,
            width=26,
        )
        friction_max_iterations_entry.grid(row=10, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="f Loop Tolerance (−)").grid(
            row=11, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=colebrook_tol_var, width=26).grid(
            row=11, column=1, sticky="ew", pady=4
        )

        ttk.Label(frame, text="V* Loop Tolerance (m/s)").grid(
            row=12, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=velocity_loop_tol_var, width=26).grid(
            row=12, column=1, sticky="ew", pady=4
        )

        ttk.Separator(frame, orient="horizontal").grid(
            row=13, column=0, columnspan=2, sticky="ew", pady=(6, 2)
        )
        ttk.Label(frame, text="— Convergence criteria —", foreground="gray").grid(
            row=14, column=0, columnspan=2, pady=(0, 4)
        )

        ttk.Label(frame, text="ΔP Correction Tol (Pa)").grid(
            row=15, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=dp_tol_var, width=26).grid(
            row=15, column=1, sticky="ew", pady=4
        )

        ttk.Label(frame, text="Continuity Tol (−)").grid(
            row=16, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=continuity_tol_var, width=26).grid(
            row=16, column=1, sticky="ew", pady=4
        )

        ttk.Separator(frame, orient="horizontal").grid(
            row=17, column=0, columnspan=2, sticky="ew", pady=(6, 2)
        )
        ttk.Label(frame, text="— Non-isothermal energy —", foreground="gray").grid(
            row=18, column=0, columnspan=2, pady=(0, 4)
        )

        ttk.Label(frame, text="Convection Scheme").grid(
            row=19, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Combobox(
            frame,
            textvariable=convection_scheme_var,
            state="readonly",
            values=("upwind", "hybrid", "power_law"),
            width=24,
        ).grid(row=19, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Max Temperature Iterations").grid(
            row=20, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=max_temp_iter_var, width=26).grid(
            row=20, column=1, sticky="ew", pady=4
        )

        ttk.Label(frame, text="Temperature Tol (K)").grid(
            row=21, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=temp_tol_var, width=26).grid(
            row=21, column=1, sticky="ew", pady=4
        )

        ttk.Label(frame, text="Temperature Relaxation").grid(
            row=22, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(frame, textvariable=temp_relax_var, width=26).grid(
            row=22, column=1, sticky="ew", pady=4
        )

        def apply_velocity_method_selection(_event: tk.Event | None = None) -> None:
            velocity_name_var.set(
                self.VELOCITY_LOOP_METHOD_LIBRARY[velocity_method_var.get()]["name"]
            )

        def apply_colebrook_strategy_selection(_event: tk.Event | None = None) -> None:
            colebrook_strategy_name_var.set(
                self.COLEBROOK_FRICTION_STRATEGY_LIBRARY[colebrook_strategy_var.get()]["name"]
            )

        def apply_friction_method_selection(_event: tk.Event | None = None) -> None:
            friction_name_var.set(
                self.FRICTION_FACTOR_METHOD_LIBRARY[friction_method_var.get()]["name"]
            )

        velocity_method_box.bind("<<ComboboxSelected>>", apply_velocity_method_selection)
        colebrook_strategy_box.bind("<<ComboboxSelected>>", apply_colebrook_strategy_selection)
        friction_method_box.bind("<<ComboboxSelected>>", apply_friction_method_selection)

        button_row = ttk.Frame(frame)
        button_row.grid(row=23, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right")
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_numerics_definition(
                dialog,
                laminar_iterations_var,
                turbulent_iterations_var,
                alpha_var,
                colebrook_strategy_var,
                friction_method_var,
                friction_max_iterations_var,
                velocity_method_var,
                velocity_max_iterations_var,
                colebrook_tol_var,
                velocity_loop_tol_var,
                dp_tol_var,
                continuity_tol_var,
                convection_scheme_var,
                max_temp_iter_var,
                temp_tol_var,
                temp_relax_var,
            ),
        ).pack(side="right", padx=(0, 8))

        frame.columnconfigure(1, weight=1)
        dialog.update_idletasks()
        dialog.wait_visibility()
        dialog.grab_set()
        dialog.focus_set()
        laminar_iterations_entry.focus_set()

    def _save_material_definition(
        self,
        dialog: tk.Toplevel,
        mode_var: tk.StringVar,
        library_var: tk.StringVar,
        name_var: tk.StringVar,
        density_var: tk.StringVar,
        viscosity_var: tk.StringVar,
        api_var: tk.StringVar,
        temperature_var: tk.StringVar,
        cp_var: tk.StringVar,
        k_var: tk.StringVar,
        mw_var: tk.StringVar | None = None,
        tc_var: tk.StringVar | None = None,
        pc_var: tk.StringVar | None = None,
        omega_var: tk.StringVar | None = None,
    ) -> None:
        mode = mode_var.get()
        is_non_isothermal = self.scene.physics_mode == "non_isothermal"

        if mode == "gas":
            mw_text = (mw_var.get() if mw_var else "").strip()
            visc_text = viscosity_var.get().strip()
            try:
                mw = float(mw_text)
                visc = float(visc_text)
            except ValueError:
                messagebox.showerror(
                    "Invalid material",
                    "Molecular weight and viscosity must be valid numbers.",
                    parent=dialog,
                )
                return
            if mw <= 0:
                messagebox.showerror(
                    "Invalid material",
                    "Molecular weight must be positive.",
                    parent=dialog,
                )
                return
            material: dict = {
                "definition_mode": "gas",
                "library_key": "",
                "name": name_var.get().strip() or "Gas",
                "molecular_weight_kg_per_mol": str(mw),
                "viscosity_pa_s": str(visc),
            }
            cp_text = cp_var.get().strip()
            k_text = k_var.get().strip()
            if cp_text:
                material["specific_heat_j_per_kg_k"] = cp_text
            if k_text:
                material["thermal_conductivity_w_per_m_k"] = k_text
            self.scene.update_material(material)
            self._refresh_global_summaries()
            self.status_var.set(f"Material set to {material['name']}.")
            dialog.destroy()
            return

        if mode == "gas_pr":
            mw_text = (mw_var.get() if mw_var else "").strip()
            visc_text = viscosity_var.get().strip()
            tc_text = (tc_var.get() if tc_var else "").strip()
            pc_text = (pc_var.get() if pc_var else "").strip()
            omega_text = (omega_var.get() if omega_var else "").strip()
            try:
                mw = float(mw_text)
                visc = float(visc_text)
                tc = float(tc_text)
                pc = float(pc_text)
                omega = float(omega_text)
            except ValueError:
                messagebox.showerror(
                    "Invalid material",
                    "Molecular weight, viscosity, Tc, Pc and ω must be valid numbers.",
                    parent=dialog,
                )
                return
            if mw <= 0 or tc <= 0 or pc <= 0:
                messagebox.showerror(
                    "Invalid material",
                    "Molecular weight, Tc and Pc must be positive.",
                    parent=dialog,
                )
                return
            material = {
                "definition_mode": "gas_pr",
                "library_key": "",
                "name": name_var.get().strip() or "Gas (PR)",
                "molecular_weight_kg_per_mol": str(mw),
                "viscosity_pa_s": str(visc),
                "critical_temperature_k": str(tc),
                "critical_pressure_pa": str(pc),
                "acentric_factor": str(omega),
            }
            cp_text2 = cp_var.get().strip()
            k_text2 = k_var.get().strip()
            if cp_text2:
                material["specific_heat_j_per_kg_k"] = cp_text2
            if k_text2:
                material["thermal_conductivity_w_per_m_k"] = k_text2
            self.scene.update_material(material)
            self._refresh_global_summaries()
            self.status_var.set(f"Material set to {material['name']}.")
            dialog.destroy()
            return

        if mode == "crude_oil":
            try:
                api = float(api_var.get().strip())
                temp = float(temperature_var.get().strip())
            except ValueError:
                messagebox.showerror(
                    "Invalid material",
                    "API gravity and temperature must be valid numbers.",
                    parent=dialog,
                )
                return
            from angelica.properties.dead_oil import dead_oil_density_kg_per_m3, dead_oil_viscosity_pa_s
            try:
                density = dead_oil_density_kg_per_m3(api)
                viscosity = dead_oil_viscosity_pa_s(api, temp)
            except ValueError as exc:
                messagebox.showerror("Invalid material", str(exc), parent=dialog)
                return
            material = {
                "definition_mode": "crude_oil",
                "library_key": "",
                "name": f"Crude oil ({api:.1f}°API, {temp:.1f}°C)",
                "api_gravity": str(api),
                "temperature_c": str(temp),
                "density_kg_per_m3": f"{density:.4f}",
                "viscosity_pa_s": f"{viscosity:.6f}",
            }
            if is_non_isothermal:
                if not self._validate_and_attach_thermal_props(
                    material, cp_var, k_var, dialog
                ):
                    return
            self.scene.update_material(material)
            self._refresh_global_summaries()
            self.status_var.set(f"Material set to {material['name']}.")
            dialog.destroy()
            return

        try:
            float(density_var.get().strip())
            float(viscosity_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid material",
                "Density and viscosity must be valid numbers.",
                parent=dialog,
            )
            return

        material = {
            "definition_mode": mode,
            "library_key": library_var.get().strip() if mode == "library" else "",
            "name": name_var.get().strip(),
            "density_kg_per_m3": density_var.get().strip(),
            "viscosity_pa_s": viscosity_var.get().strip(),
        }
        if mode == "library" and not material["library_key"]:
            messagebox.showerror(
                "Invalid material",
                "Select a material from the library.",
                parent=dialog,
            )
            return
        if not material["name"]:
            messagebox.showerror(
                "Invalid material",
                "Material name cannot be empty.",
                parent=dialog,
            )
            return
        if is_non_isothermal:
            if not self._validate_and_attach_thermal_props(
                material, cp_var, k_var, dialog
            ):
                return
        self.scene.update_material(material)
        self._refresh_global_summaries()
        self.status_var.set(f"Material set to {material['name']}.")
        dialog.destroy()

    def _validate_and_attach_thermal_props(
        self,
        material: dict,
        cp_var: tk.StringVar,
        k_var: tk.StringVar,
        dialog: tk.Toplevel,
    ) -> bool:
        cp_text = cp_var.get().strip()
        k_text = k_var.get().strip()
        try:
            float(cp_text)
            float(k_text)
        except ValueError:
            messagebox.showerror(
                "Invalid material",
                "Specific heat (cp) and thermal conductivity (k) must be valid numbers.",
                parent=dialog,
            )
            return False
        material["specific_heat_j_per_kg_k"] = cp_text
        material["thermal_conductivity_w_per_m_k"] = k_text
        return True

    def _save_pressure_drop_model_definition(
        self,
        dialog: tk.Toplevel,
        library_var: tk.StringVar,
    ) -> None:
        model_key = library_var.get().strip()
        if model_key not in self.PRESSURE_DROP_MODEL_LIBRARY:
            messagebox.showerror(
                "Invalid model",
                "Select a valid pressure-drop model from the library.",
                parent=dialog,
            )
            return

        model_definition = {
            "library_key": model_key,
            "name": self.PRESSURE_DROP_MODEL_LIBRARY[model_key]["name"],
        }
        self.scene.update_pressure_drop_model(model_definition)
        self._refresh_global_summaries()
        self.status_var.set(f"Pipe pressure-drop model set to {model_definition['name']}.")
        dialog.destroy()

    def _save_numerics_definition(
        self,
        dialog: tk.Toplevel,
        laminar_iterations_var: tk.StringVar,
        turbulent_iterations_var: tk.StringVar,
        alpha_var: tk.StringVar,
        colebrook_strategy_var: tk.StringVar,
        friction_method_var: tk.StringVar,
        friction_max_iterations_var: tk.StringVar,
        velocity_method_var: tk.StringVar,
        velocity_max_iterations_var: tk.StringVar,
        colebrook_tol_var: tk.StringVar,
        velocity_loop_tol_var: tk.StringVar,
        dp_tol_var: tk.StringVar,
        continuity_tol_var: tk.StringVar,
        convection_scheme_var: tk.StringVar,
        max_temp_iter_var: tk.StringVar,
        temp_tol_var: tk.StringVar,
        temp_relax_var: tk.StringVar,
    ) -> None:
        laminar_iterations_text = laminar_iterations_var.get().strip()
        if laminar_iterations_text.lower() in ("", "auto"):
            laminar_iterations = None
        else:
            try:
                laminar_iterations = int(laminar_iterations_text)
            except ValueError:
                messagebox.showerror(
                    "Invalid numerics",
                    "Laminar iterations must be an integer or 'auto'.",
                    parent=dialog,
                )
                return
            if laminar_iterations <= 0:
                messagebox.showerror(
                    "Invalid numerics",
                    "Laminar iterations must be greater than zero.",
                    parent=dialog,
                )
                return

        try:
            turbulent_iterations = int(turbulent_iterations_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Turbulent iterations must be an integer.",
                parent=dialog,
            )
            return

        if turbulent_iterations <= 0:
            messagebox.showerror(
                "Invalid numerics",
                "Turbulent iterations must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            alpha = float(alpha_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Pressure relaxation must be a valid number.",
                parent=dialog,
            )
            return

        if alpha <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "Pressure relaxation must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            friction_max_iterations = int(friction_max_iterations_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Friction max iterations must be an integer.",
                parent=dialog,
            )
            return

        if friction_max_iterations <= 0:
            messagebox.showerror(
                "Invalid numerics",
                "Friction max iterations must be greater than zero.",
                parent=dialog,
            )
            return

        velocity_method = velocity_method_var.get().strip()
        if velocity_method not in self.VELOCITY_LOOP_METHOD_LIBRARY:
            messagebox.showerror(
                "Invalid numerics",
                "Select a valid velocity loop method.",
                parent=dialog,
            )
            return

        try:
            velocity_max_iterations = int(velocity_max_iterations_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Velocity max iterations must be an integer.",
                parent=dialog,
            )
            return

        if velocity_max_iterations <= 0:
            messagebox.showerror(
                "Invalid numerics",
                "Velocity max iterations must be greater than zero.",
                parent=dialog,
            )
            return

        colebrook_strategy = colebrook_strategy_var.get().strip()
        if colebrook_strategy not in self.COLEBROOK_FRICTION_STRATEGY_LIBRARY:
            messagebox.showerror(
                "Invalid numerics",
                "Select a valid Colebrook friction strategy.",
                parent=dialog,
            )
            return

        friction_method = friction_method_var.get().strip()
        if friction_method not in self.FRICTION_FACTOR_METHOD_LIBRARY:
            messagebox.showerror(
                "Invalid numerics",
                "Select a valid friction-factor method.",
                parent=dialog,
            )
            return

        try:
            colebrook_tol = float(colebrook_tol_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "f loop tolerance must be a valid number.",
                parent=dialog,
            )
            return

        if colebrook_tol <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "f loop tolerance must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            velocity_loop_tol = float(velocity_loop_tol_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "V* loop tolerance must be a valid number.",
                parent=dialog,
            )
            return

        if velocity_loop_tol <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "V* loop tolerance must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            dp_tol = float(dp_tol_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "ΔP correction tolerance must be a valid number.",
                parent=dialog,
            )
            return

        if dp_tol <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "ΔP correction tolerance must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            continuity_tol = float(continuity_tol_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Continuity tolerance must be a valid number.",
                parent=dialog,
            )
            return

        if continuity_tol <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "Continuity tolerance must be greater than zero.",
                parent=dialog,
            )
            return

        convection_scheme = convection_scheme_var.get().strip()
        if convection_scheme not in ("upwind", "hybrid", "power_law"):
            messagebox.showerror(
                "Invalid numerics",
                "Select a valid convection scheme (upwind, hybrid, power_law).",
                parent=dialog,
            )
            return

        try:
            max_temp_iter = int(max_temp_iter_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Max temperature iterations must be an integer.",
                parent=dialog,
            )
            return
        if max_temp_iter <= 0:
            messagebox.showerror(
                "Invalid numerics",
                "Max temperature iterations must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            temp_tol = float(temp_tol_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Temperature tolerance must be a valid number.",
                parent=dialog,
            )
            return
        if temp_tol <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "Temperature tolerance must be greater than zero.",
                parent=dialog,
            )
            return

        try:
            temp_relax = float(temp_relax_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid numerics",
                "Temperature relaxation must be a valid number.",
                parent=dialog,
            )
            return
        if temp_relax <= 0.0:
            messagebox.showerror(
                "Invalid numerics",
                "Temperature relaxation must be greater than zero.",
                parent=dialog,
            )
            return

        self.scene.update_solver_settings(
            {
                "laminar_iterations": laminar_iterations,
                "turbulent_iterations": turbulent_iterations,
                "pressure_relaxation": alpha,
                "colebrook_friction_strategy": colebrook_strategy,
                "friction_factor_method": friction_method,
                "friction_factor_max_iterations": friction_max_iterations,
                "velocity_loop_method": velocity_method,
                "velocity_loop_max_iterations": velocity_max_iterations,
                "colebrook_residual_tolerance": colebrook_tol,
                "velocity_loop_tolerance": velocity_loop_tol,
                "pressure_correction_abs_tolerance_pa": dp_tol,
                "nodal_mass_imbalance_rel_tolerance": continuity_tol,
                "convection_scheme": convection_scheme,
                "max_temperature_iterations": max_temp_iter,
                "temperature_tolerance_k": temp_tol,
                "temperature_relaxation": temp_relax,
            }
        )
        self._refresh_global_summaries()
        self.status_var.set(
            "Numerics updated: "
            f"laminar={'auto' if laminar_iterations is None else laminar_iterations}, "
            f"turbulent={turbulent_iterations}, "
            f"alpha={alpha:g}, "
            f"f-tol={self._fmt_sci(colebrook_tol)}, V*-tol={self._fmt_sci(velocity_loop_tol)}, "
            f"ΔP-tol={self._fmt_sci(dp_tol)} Pa, continuity-tol={self._fmt_sci(continuity_tol)}."
        )
        dialog.destroy()

    def _open_unit_system_dialog(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title("Unit System")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=16)
        frame.pack(fill="both", expand=True)

        unit_var = tk.StringVar(value=self._unit_system_key)
        for row, (key, info) in enumerate(self.UNIT_SYSTEMS.items()):
            ttk.Radiobutton(
                frame,
                text=info["name"],
                value=key,
                variable=unit_var,
            ).grid(row=row, column=0, sticky="w", pady=4)

        button_row = ttk.Frame(frame)
        button_row.grid(row=len(self.UNIT_SYSTEMS), column=0, sticky="e", pady=(12, 0))
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right")
        ttk.Button(
            button_row,
            text="Apply",
            command=lambda: self._apply_unit_system(unit_var.get(), dialog),
        ).pack(side="right", padx=(0, 8))

        dialog.update_idletasks()
        dialog.wait_visibility()
        dialog.grab_set()
        dialog.focus_set()

    def _apply_unit_system(self, key: str, dialog: tk.Toplevel) -> None:
        self._unit_system_key = key
        self._redraw_scene()
        self.status_var.set(f"Unit system: {self.UNIT_SYSTEMS[key]['name']}.")
        dialog.destroy()

    def _open_case_type_dialog(self) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title("Case Type")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=16)
        frame.pack(fill="both", expand=True)

        current_mode = self.scene.physics_mode
        is_currently_compressible = current_mode == "compressible"

        compressibility_var = tk.StringVar(
            value="compressible" if is_currently_compressible else "incompressible"
        )
        energy_var = tk.StringVar(
            value="non_isothermal" if current_mode == "non_isothermal" else "isothermal"
        )

        ttk.Label(frame, text="Compressibility").grid(row=0, column=0, sticky="w", pady=(0, 4))
        comp_frame = ttk.Frame(frame)
        comp_frame.grid(row=0, column=1, sticky="w", padx=(8, 0), pady=(0, 4))
        ttk.Radiobutton(
            comp_frame, text="Incompressible", variable=compressibility_var, value="incompressible"
        ).pack(side="left", padx=(0, 12))
        ttk.Radiobutton(
            comp_frame, text="Compressible", variable=compressibility_var, value="compressible"
        ).pack(side="left")

        ttk.Separator(frame, orient="horizontal").grid(
            row=1, column=0, columnspan=2, sticky="ew", pady=8
        )

        ttk.Label(frame, text="Energy").grid(row=2, column=0, sticky="w", pady=4)
        energy_frame = ttk.Frame(frame)
        energy_frame.grid(row=2, column=1, sticky="w", padx=(8, 0))
        rb_isothermal = ttk.Radiobutton(
            energy_frame, text="Isothermal", variable=energy_var, value="isothermal"
        )
        rb_isothermal.pack(side="left", padx=(0, 12))
        rb_non_isothermal = ttk.Radiobutton(
            energy_frame, text="Non-isothermal", variable=energy_var, value="non_isothermal"
        )
        rb_non_isothermal.pack(side="left")

        note_ni = ttk.Label(
            frame,
            text=(
                "Non-isothermal requires cp and k in the material,\n"
                "inlet temperature on source nodes, and U / T_amb on pipes.\n"
                "Compressible requires molecular weight (M) in the material."
            ),
            foreground="gray",
            justify="left",
        )
        note_ni.grid(row=3, column=0, columnspan=2, sticky="w", pady=(4, 8))

        def _sync_energy_state(*_args: object) -> None:
            state = "disabled" if compressibility_var.get() == "compressible" else "normal"
            rb_isothermal.configure(state=state)
            rb_non_isothermal.configure(state=state)

        compressibility_var.trace_add("write", _sync_energy_state)
        _sync_energy_state()

        button_row = ttk.Frame(frame)
        button_row.grid(row=4, column=0, columnspan=2, sticky="e", pady=(4, 0))
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right")
        ttk.Button(
            button_row,
            text="Apply",
            command=lambda: self._apply_case_type(
                compressibility_var.get(), energy_var.get(), dialog
            ),
        ).pack(side="right", padx=(0, 8))

        frame.columnconfigure(1, weight=1)
        dialog.update_idletasks()
        dialog.wait_visibility()
        dialog.grab_set()
        dialog.focus_set()

    def _apply_case_type(
        self, compressibility: str, energy: str, dialog: tk.Toplevel
    ) -> None:
        if compressibility == "compressible":
            mode = "compressible"
        elif energy == "non_isothermal":
            mode = "non_isothermal"
        else:
            mode = "isothermal"

        self.scene.physics_mode = mode
        if mode == "non_isothermal" and self.scene.material:
            library_key = self.scene.material.get("library_key", "")
            if library_key in self.MATERIAL_LIBRARY:
                lib = self.MATERIAL_LIBRARY[library_key]
                updated = dict(self.scene.material)
                changed = False
                for field in ("specific_heat_j_per_kg_k", "thermal_conductivity_w_per_m_k"):
                    if not updated.get(field, "").strip() and lib.get(field, ""):
                        updated[field] = lib[field]
                        changed = True
                if changed:
                    self.scene.update_material(updated)
                    self._refresh_global_summaries()

        _LABELS = {
            "isothermal": "Isothermal",
            "non_isothermal": "Non-isothermal",
            "compressible": "Compressible",
        }
        self.status_var.set(f"Case type set to: {_LABELS.get(mode, mode)}.")
        dialog.destroy()

    def _unit_quantities(self) -> dict:
        return self.UNIT_SYSTEMS[self._unit_system_key]["quantities"]

    def _unit_label(self, quantity: str) -> str:
        return self._unit_quantities()[quantity][0]

    def _default_open_dir(self) -> str:
        """Start the Open dialog in the tutorials folder when one can be found."""
        if getattr(sys, "frozen", False):
            exe_dir = os.path.dirname(sys.executable)
            tutorials = os.path.join(exe_dir, "tutorials")
            if os.path.isdir(tutorials):
                return tutorials
        else:
            # Dev / editable install: tutorials/ is 3 dirs above src/angelica/gui/
            here = os.path.dirname(os.path.abspath(__file__))
            candidate = os.path.normpath(os.path.join(here, "..", "..", "..", "tutorials"))
            if os.path.isdir(candidate):
                return candidate
        return os.path.expanduser("~")

    def _set_window_icon(self) -> None:
        icon_path = self._resource_path("angelica_32.png")
        if not os.path.exists(icon_path):
            return
        try:
            photo = tk.PhotoImage(file=icon_path)
            self.root.iconphoto(True, photo)
            self._icon_photo = photo  # keep a reference so GC doesn't collect it
        except Exception:
            pass

    @staticmethod
    def _resource_path(filename: str) -> str:
        """Resolve a bundled data file: PyInstaller _MEIPASS at runtime, repo root in dev."""
        if getattr(sys, "frozen", False):
            return os.path.join(sys._MEIPASS, filename)
        here = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(here, "..", "..", "..", "..", "installer", filename)

    @staticmethod
    def _fmt(value: float) -> str:
        """European format: period thousands separator, comma decimal, 2 decimal places."""
        s = f"{value:,.2f}"
        return s.replace(",", "\x00").replace(".", ",").replace("\x00", ".")

    @staticmethod
    def _fmt_sci(value: float) -> str:
        """Scientific notation without decimals: 1e-4, 5e-3, 1e+3."""
        s = f"{float(value):.0e}"
        mantissa, exp_part = s.split("e")
        sign = exp_part[0]
        exp_val = str(int(exp_part[1:]))
        return f"{mantissa}e{sign}{exp_val}"

    def _si_to_display(self, si_str: str, quantity: str) -> str:
        if not si_str.strip():
            return si_str
        try:
            si_val = float(si_str)
        except ValueError:
            return si_str
        _unit, factor = self._unit_quantities()[quantity]
        return f"{si_val * factor:.6g}"

    def _display_to_si(self, display_str: str, quantity: str) -> str:
        if not display_str.strip():
            return display_str
        try:
            display_val = float(display_str)
        except ValueError:
            return display_str
        _unit, factor = self._unit_quantities()[quantity]
        return repr(display_val / factor)

    def _field_from_si(self, field_name: str, si_str: str) -> str:
        quantity = self._FIELD_QUANTITY.get(field_name)
        return self._si_to_display(si_str, quantity) if quantity else si_str

    def _field_to_si(self, field_name: str, display_str: str) -> str:
        quantity = self._FIELD_QUANTITY.get(field_name)
        return self._display_to_si(display_str, quantity) if quantity else display_str

    def _refresh_global_summaries(self) -> None:
        self.material_summary_var.set(self._material_summary_text())
        self.pressure_drop_summary_var.set(self._pressure_drop_summary_text())
        self.numerics_summary_var.set(self._numerics_summary_text())

    def _material_summary_text(self) -> str:
        if not self.scene.material:
            return "Not defined"

        name = self.scene.material.get("name", "Unnamed")
        lines = [name]
        if self.scene.physics_mode == "compressible":
            mw = self.scene.material.get("molecular_weight_kg_per_mol", "").strip()
            viscosity = self.scene.material.get("viscosity_pa_s", "").strip()
            if mw:
                lines.append(f"M={mw} kg/mol")
            if viscosity:
                lines.append(f"mu={viscosity} Pa·s")
        else:
            density = self.scene.material.get("density_kg_per_m3", "").strip()
            viscosity = self.scene.material.get("viscosity_pa_s", "").strip()
            if density:
                lines.append(f"rho={density} kg/m^3")
            if viscosity:
                lines.append(f"mu={viscosity} Pa·s")
        return "\n".join(lines)

    def _pressure_drop_summary_text(self) -> str:
        model_name = self.scene.pressure_drop_model.get("name", "").strip()
        if model_name:
            return model_name
        return "Not defined"

    def _numerics_summary_text(self) -> str:
        laminar_iterations = self.scene.solver_settings.get("laminar_iterations")
        turbulent_iterations = self.scene.solver_settings.get("turbulent_iterations", 60)
        alpha = self.scene.solver_settings.get("pressure_relaxation", 1.0)
        velocity_method = str(
            self.scene.solver_settings.get("velocity_loop_method", "fixed_point")
        )
        colebrook_strategy = str(
            self.scene.solver_settings.get("colebrook_friction_strategy", "transformed")
        )
        friction_method = str(
            self.scene.solver_settings.get("friction_factor_method", "newton")
        )
        friction_max_iterations = self.scene.solver_settings.get(
            "friction_factor_max_iterations",
            50,
        )
        friction_method_name = self.FRICTION_FACTOR_METHOD_LIBRARY.get(
            friction_method,
            {},
        ).get("name", friction_method)
        colebrook_strategy_name = self.COLEBROOK_FRICTION_STRATEGY_LIBRARY.get(
            colebrook_strategy,
            {},
        ).get("name", colebrook_strategy)
        velocity_method_name = self.VELOCITY_LOOP_METHOD_LIBRARY.get(
            velocity_method,
            {},
        ).get("name", velocity_method)
        velocity_max_iterations = self.scene.solver_settings.get(
            "velocity_loop_max_iterations",
            50,
        )
        colebrook_tol = self.scene.solver_settings.get("colebrook_residual_tolerance", 1e-4)
        velocity_loop_tol = self.scene.solver_settings.get("velocity_loop_tolerance", 1e-4)
        dp_tol = self.scene.solver_settings.get("pressure_correction_abs_tolerance_pa", 1e-3)
        continuity_tol = self.scene.solver_settings.get("nodal_mass_imbalance_rel_tolerance", 1e-3)
        return (
            f"laminar={'auto' if laminar_iterations is None else laminar_iterations}\n"
            f"turbulent={turbulent_iterations}\n"
            f"alpha={alpha}\n"
            f"colebrook={colebrook_strategy_name}\n"
            f"friction={friction_method_name} ({friction_max_iterations})\n"
            f"velocity={velocity_method_name} ({velocity_max_iterations})\n"
            f"f-tol={self._fmt_sci(colebrook_tol)}\n"
            f"V*-tol={self._fmt_sci(velocity_loop_tol)} m/s\n"
            f"ΔP-tol={self._fmt_sci(dp_tol)} Pa\n"
            f"cont-tol={self._fmt_sci(continuity_tol)}"
        )

    def _validate_scene(self) -> list[str]:
        errors: list[str] = []
        scene = self.scene

        if not scene.nodes:
            return ["The network has no nodes."]
        if not scene.links:
            return ["The network has no connections between nodes."]

        # adjacency map from canvas node IDs only (no internal solver nodes)
        adjacency: dict[int, set[int]] = {node.node_id: set() for node in scene.nodes}
        linked_ids: set[int] = set()
        for link in scene.links:
            linked_ids.add(link.start_node_id)
            linked_ids.add(link.end_node_id)
            adjacency[link.start_node_id].add(link.end_node_id)
            adjacency[link.end_node_id].add(link.start_node_id)

        # isolated nodes (not in any link)
        for node in scene.nodes:
            if node.node_id not in linked_ids:
                errors.append(
                    f"{node.node_type.capitalize()} #{node.node_id} is not connected to any other node."
                )

        # fully connected graph — BFS from the first node
        if not errors:
            all_ids = {node.node_id for node in scene.nodes}
            start = next(iter(all_ids))
            visited: set[int] = {start}
            queue = [start]
            while queue:
                current = queue.pop()
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            unreachable = all_ids - visited
            if unreachable:
                errors.append(
                    "The network has disconnected sub-networks. "
                    "All nodes must form a single connected graph."
                )

        # at least one pressure boundary condition
        has_pressure_bc = any(
            node.node_type in {"source", "sink"}
            and node.properties.get("condition_type", "pressure") == "pressure"
            for node in scene.nodes
        )
        if not has_pressure_bc:
            errors.append(
                "No pressure boundary condition found. "
                "At least one source or sink must have a fixed pressure."
            )

        return errors

    def _run_simulation(self) -> None:
        errors = self._validate_scene()
        if errors:
            messagebox.showerror(
                "Cannot run simulation",
                "\n\n".join(f"• {e}" for e in errors),
            )
            return

        try:
            case = self._build_network_case_from_scene()
        except ValueError as exc:
            messagebox.showerror("Run failed", str(exc))
            return

        solver = build_solver_from_scene(self.scene)
        self.convergence_history = {"laminar": [], "turbulent": []}
        self._prepare_convergence_window()
        result = solver.solve(case, progress_callback=self._on_solver_progress)
        self.latest_result = result
        self.latest_boundary_results = self._build_boundary_results(case, result)
        self.convergence_history = {
            "laminar": list(result.laminar_metrics),
            "turbulent": list(result.turbulent_metrics),
        }
        self.temperature_history = list(result.temperature_history)
        self.density_history = list(result.density_history)
        self.outer_turbulent_final_metrics = list(result.outer_turbulent_final_metrics)
        self._redraw_scene()
        self._redraw_convergence_plot()

        if result.converged:
            self.status_var.set(f"Simulation converged for case '{case.name}'.")
        else:
            self.status_var.set(f"Simulation did not converge for case '{case.name}'.")

    def _export_results_report(self) -> None:
        if self.latest_result is None:
            messagebox.showerror(
                "Export failed",
                "Run a simulation before exporting a results report.",
            )
            return
        if not self.latest_result.converged:
            messagebox.showerror(
                "Export failed",
                "Only converged simulations can be exported to a results report.",
            )
            return

        file_path = filedialog.asksaveasfilename(
            title="Export Results Report",
            defaultextension=".xlsx",
            filetypes=(
                ("Excel workbook", "*.xlsx"),
                ("All files", "*.*"),
            ),
        )
        if not file_path:
            return

        try:
            export_solve_result_workbook(self.latest_result, file_path)
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Export failed", f"Could not export report:\n{exc}")
            return

        self.status_var.set(f"Results report exported to {file_path}.")

    def _export_results_report_csv(self) -> None:
        if self.latest_result is None:
            messagebox.showerror(
                "Export failed",
                "Run a simulation before exporting a results report.",
            )
            return
        if not self.latest_result.converged:
            messagebox.showerror(
                "Export failed",
                "Only converged simulations can be exported to a results report.",
            )
            return

        file_path = filedialog.asksaveasfilename(
            title="Export Results Report (CSV)",
            defaultextension=".csv",
            filetypes=(
                ("CSV file", "*.csv"),
                ("All files", "*.*"),
            ),
        )
        if not file_path:
            return

        try:
            export_solve_result_csv(self.latest_result, file_path)
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Export failed", f"Could not export report:\n{exc}")
            return

        self.status_var.set(f"Results report exported to {file_path}.")

    def _on_canvas_press(self, event: tk.Event) -> None:
        node_id = self._node_id_at(event.x, event.y)

        if self.scene.active_tool is not None:
            return

        if node_id is None:
            if self.selected_node_id is not None:
                self.selected_node_id = None
                self._redraw_scene()
            return

        node = self.scene.get_node(node_id)
        if node is None:
            return

        prev_selected = self.selected_node_id
        self.selected_node_id = node_id
        self.moving_node_id = node_id
        self._move_pre_snapshot = self._make_snapshot()
        self.canvas.focus_set()
        if self.selected_node_id != prev_selected:
            self._redraw_scene()
        self.status_var.set(
            f"Selected {node.node_type} #{node.node_id}. "
            "Drag to move · Delete/Backspace to remove."
        )

    def _on_canvas_right_press(self, event: tk.Event) -> None:
        if self.scene.active_tool is not None:
            return

        node_id = self._node_id_at(event.x, event.y)
        if node_id is not None:
            node = self.scene.get_node(node_id)
            if node is None:
                return
            self.drag_source_node_id = node_id
            source_x, source_y = self._scene_to_canvas(node.x, node.y)
            self.drag_line_id = self.canvas.create_line(
                source_x,
                source_y,
                event.x,
                event.y,
                fill=self._t["drag_line"],
                width=2,
                dash=(5, 3),
            )
            self.status_var.set(
                f"Connecting from {node.node_type} #{node.node_id}. Drag with right mouse button."
            )
            return

        link_id = self._link_id_at(event.x, event.y)
        if link_id is None:
            return
        link = self.scene.get_link(link_id)
        if link is None:
            return
        self._show_link_context_menu(event, link)

    def _on_canvas_shift_press(self, event: tk.Event) -> None:
        self._on_canvas_right_press(event)

    def _on_canvas_drag(self, event: tk.Event) -> None:
        if self.moving_node_id is None:
            return

        scene_x, scene_y = self._canvas_to_scene(event.x, event.y)
        updated_node = self.scene.move_node(self.moving_node_id, scene_x, scene_y)
        self._redraw_scene()
        self.status_var.set(
            f"Moving {updated_node.node_type} #{updated_node.node_id} to "
            f"({int(updated_node.x)}, {int(updated_node.y)})."
        )

    def _on_canvas_right_drag(self, event: tk.Event) -> None:
        if self.drag_line_id is None or self.drag_source_node_id is None:
            return

        source = self.scene.get_node(self.drag_source_node_id)
        if source is None:
            return

        source_x, source_y = self._scene_to_canvas(source.x, source.y)
        self.canvas.coords(self.drag_line_id, source_x, source_y, event.x, event.y)

    def _on_canvas_shift_drag(self, event: tk.Event) -> None:
        self._on_canvas_right_drag(event)

    def _on_canvas_middle_press(self, event: tk.Event) -> None:
        self.middle_pan_anchor = (float(event.x), float(event.y))
        self.status_var.set("Panning view. Drag with middle mouse button.")

    def _on_canvas_middle_drag(self, event: tk.Event) -> None:
        if self.middle_pan_anchor is None:
            return

        anchor_x, anchor_y = self.middle_pan_anchor
        self.view_offset_x += float(event.x) - anchor_x
        self.view_offset_y += float(event.y) - anchor_y
        self.middle_pan_anchor = (float(event.x), float(event.y))
        self._redraw_scene()

    def _on_canvas_middle_release(self, _event: tk.Event) -> None:
        self.middle_pan_anchor = None
        self.status_var.set("View panned.")

    def _on_canvas_scroll(self, event: tk.Event) -> None:
        if event.num == 4 or (hasattr(event, "delta") and event.delta > 0):
            factor = 1.1
        elif event.num == 5 or (hasattr(event, "delta") and event.delta < 0):
            factor = 1.0 / 1.1
        else:
            return
        cx, cy = float(event.x), float(event.y)
        self.view_offset_x = cx - (cx - self.view_offset_x) * factor
        self.view_offset_y = cy - (cy - self.view_offset_y) * factor
        self.view_scale *= factor
        self._redraw_scene()

    def _on_delete_key(self, _event: tk.Event) -> None:
        # Ignore if a dialog or entry widget has focus
        focused = self.root.focus_get()
        if focused is not None and focused is not self.canvas and focused is not self.root:
            return
        self._delete_selected_node()

    def _delete_selected_node(self) -> None:
        if self.selected_node_id is None:
            return
        node = self.scene.get_node(self.selected_node_id)
        if node is None:
            self.selected_node_id = None
            return
        node_id = self.selected_node_id
        node_type = node.node_type
        self._push_undo()
        self.scene.remove_node(node_id)
        self.selected_node_id = None
        self.latest_boundary_results = {}
        self.latest_result = None
        self._redraw_scene()
        self.status_var.set(f"Deleted {node_type} #{node_id} and its connections.")

    def _show_link_context_menu(self, event: tk.Event, link: CanvasLink) -> None:
        menu = tk.Menu(self.root, tearoff=False)
        menu.add_command(
            label=f"Edit Components — Connection #{link.link_id}",
            command=lambda: self._open_link_properties_dialog(link),
        )
        if self.latest_boundary_results:
            menu.add_command(
                label="Pressure Profile",
                command=lambda: self._show_link_pressure_profile(link),
            )
        if self._has_temperature_results():
            menu.add_command(
                label="Temperature Profile",
                command=lambda: self._show_link_temperature_profile(link),
            )
        menu.add_separator()
        menu.add_command(
            label=f"Delete Connection #{link.link_id}",
            command=lambda: self._delete_link(link.link_id),
        )
        menu.tk_popup(event.x_root, event.y_root)

    def _delete_link(self, link_id: int) -> None:
        self._push_undo()
        self.scene.remove_link(link_id)
        self.latest_boundary_results = {}
        self.latest_result = None
        self._redraw_scene()
        self.status_var.set(f"Deleted connection #{link_id}.")

    def _open_convergence_window(self) -> None:
        self._prepare_convergence_window()

    def _on_canvas_release(self, event: tk.Event) -> None:
        if self.moving_node_id is not None:
            moved_node = self.scene.get_node(self.moving_node_id)
            self.moving_node_id = None
            self._commit_move_if_changed()
            if moved_node is not None:
                self.status_var.set(
                    f"Moved {moved_node.node_type} #{moved_node.node_id} to "
                    f"({int(moved_node.x)}, {int(moved_node.y)})."
                )
            return

        if self.scene.active_tool is None:
            self.status_var.set("Select Source, Sink, or Junction before placing a node.")
            return

        if self._node_id_at(event.x, event.y) is not None:
            self.status_var.set("Release on empty canvas space to place a new node.")
            return

        self._push_undo()
        scene_x, scene_y = self._canvas_to_scene(event.x, event.y)
        node = self.scene.add_node(scene_x, scene_y)
        self._draw_node(node)
        placed_tool = node.node_type
        self.scene.set_active_tool(None)
        self.tool_var.set("No tool selected")
        self.status_var.set(
            f"Placed {placed_tool} #{node.node_id} at ({int(node.x)}, {int(node.y)}). "
            "Select a node type to place another."
        )

    def _on_canvas_right_release(self, event: tk.Event) -> None:
        if self.drag_source_node_id is not None:
            self._finish_connection(event)

    def _on_canvas_shift_release(self, event: tk.Event) -> None:
        self._on_canvas_right_release(event)

    def _on_canvas_double_click(self, event: tk.Event) -> None:
        node_id = self._node_id_at(event.x, event.y)
        if node_id is not None:
            node = self.scene.get_node(node_id)
            if node is None:
                return

            self._open_node_properties_dialog(node)
            return

        link_id = self._link_id_at(event.x, event.y)
        if link_id is None:
            return

        link = self.scene.get_link(link_id)
        if link is None:
            return

        if self.latest_boundary_results:
            self._show_link_pressure_profile(link)
        else:
            self._open_link_properties_dialog(link)

    def _finish_connection(self, event: tk.Event) -> None:
        source_node_id = self.drag_source_node_id
        target_node_id = self._node_id_at(event.x, event.y)
        self._clear_drag_line()

        if source_node_id is None:
            return
        if target_node_id is None:
            self.status_var.set("Connection cancelled. Release over another node to connect.")
            return
        if target_node_id == source_node_id:
            self.status_var.set("Connection cancelled. Choose a different target node.")
            return

        self._push_undo()
        try:
            link = self.scene.add_link(source_node_id, target_node_id)
        except ValueError as exc:
            self._undo_stack.pop()
            self.status_var.set(str(exc))
            return

        self._draw_link(link)
        source = self.scene.get_node(source_node_id)
        target = self.scene.get_node(target_node_id)
        if source is not None and target is not None:
            self.status_var.set(
                f"Connected {source.node_type} #{source.node_id} to "
                f"{target.node_type} #{target.node_id}."
            )

    def _draw_node(self, node: CanvasNode) -> None:
        s = self.view_scale
        radius = max(8, round(24 * s))
        canvas_x, canvas_y = self._scene_to_canvas(node.x, node.y)
        x0 = canvas_x - radius
        y0 = canvas_y - radius
        x1 = canvas_x + radius
        y1 = canvas_y + radius

        fill_color = self._node_fill(node.node_type)
        label = self._node_label(node)
        is_selected = node.node_id == self.selected_node_id

        node_tag = f"node_{node.node_id}"

        self.canvas.create_oval(
            x0,
            y0,
            x1,
            y1,
            fill=fill_color,
            outline="#ffcc00" if is_selected else self._t["node_outline"],
            width=max(2, round(3 * s)) if is_selected else max(1, round(2 * s)),
            tags=(node_tag, "node"),
        )
        self.canvas.create_text(
            canvas_x,
            canvas_y - round(2 * s),
            text=label,
            fill=self._t["node_text"],
            font=("TkDefaultFont", max(6, round(9 * s)), "bold"),
            tags=(node_tag, "node"),
        )
        self.canvas.create_text(
            canvas_x,
            canvas_y + round(12 * s),
            text=str(node.node_id),
            fill=self._t["node_text"],
            font=("TkDefaultFont", max(5, round(8 * s))),
            tags=(node_tag, "node"),
        )

        summary = self._node_summary_text(node)
        if summary:
            self.canvas.create_text(
                canvas_x,
                canvas_y + round(40 * s),
                text=summary,
                font=("TkDefaultFont", max(5, round(8 * s))),
                fill=self._t["node_summary"],
                justify="center",
                tags=(node_tag, "node"),
            )

    def _draw_link(self, link: CanvasLink) -> None:
        start = self.scene.get_node(link.start_node_id)
        end = self.scene.get_node(link.end_node_id)
        if start is None or end is None:
            return
        start_x, start_y = self._scene_to_canvas(start.x, start.y)
        end_x, end_y = self._scene_to_canvas(end.x, end.y)

        types = {c.component_type for c in link.components}
        if "heat_source" in types:
            link_color = "#c0397d"
        elif "pump" in types:
            link_color = self._t["plot_turbulent"]
        else:
            link_color = self._t["link"]

        self.canvas.create_line(
            start_x,
            start_y,
            end_x,
            end_y,
            fill=link_color,
            width=3,
            tags=("link", f"link_{link.link_id}"),
        )
        self.canvas.tag_lower("link")

    def _redraw_scene(self) -> None:
        self.canvas.delete("all")
        for link in self.scene.links:
            self._draw_link(link)
        for node in self.scene.nodes:
            self._draw_node(node)

    def _node_id_at(self, x: float, y: float) -> int | None:
        overlapping = self.canvas.find_overlapping(x - 1, y - 1, x + 1, y + 1)
        for item_id in reversed(overlapping):
            for tag in self.canvas.gettags(item_id):
                if tag.startswith("node_"):
                    return int(tag.split("_", 1)[1])
        return None

    def _link_id_at(self, x: float, y: float) -> int | None:
        overlapping = self.canvas.find_overlapping(x - 3, y - 3, x + 3, y + 3)
        for item_id in reversed(overlapping):
            for tag in self.canvas.gettags(item_id):
                if tag.startswith("link_"):
                    return int(tag.split("_", 1)[1])
        return None

    def _clear_drag_line(self) -> None:
        if self.drag_line_id is not None:
            self.canvas.delete(self.drag_line_id)
        self.drag_line_id = None
        self.drag_source_node_id = None

    def _scene_to_canvas(self, x: float, y: float) -> tuple[float, float]:
        return x * self.view_scale + self.view_offset_x, y * self.view_scale + self.view_offset_y

    def _canvas_to_scene(self, x: float, y: float) -> tuple[float, float]:
        return (x - self.view_offset_x) / self.view_scale, (y - self.view_offset_y) / self.view_scale

    def _open_node_properties_dialog(self, node: CanvasNode) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title(f"{node.node_type.capitalize()} #{node.node_id}")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        container = ttk.Frame(dialog, padding=12)
        container.pack(fill="both", expand=True)
        container.columnconfigure(1, weight=1)

        entries: dict[str, tk.StringVar] = {}

        if node.node_type in {"source", "sink"}:
            ttk.Label(container, text="Boundary Type").grid(row=0, column=0, sticky="w", pady=4)
            condition_var = tk.StringVar(value=node.properties.get("condition_type", "pressure"))
            entries["condition_type"] = condition_var

            condition_frame = ttk.Frame(container)
            condition_frame.grid(row=0, column=1, sticky="w", pady=4)

            ttk.Radiobutton(
                condition_frame,
                text="Pressure",
                value="pressure",
                variable=condition_var,
            ).pack(side="left", padx=(0, 8))
            ttk.Radiobutton(
                condition_frame,
                text="Flow",
                value="flow",
                variable=condition_var,
            ).pack(side="left")

            ttk.Label(container, text=f"Pressure ({self._unit_label('pressure')})").grid(row=1, column=0, sticky="w", pady=4)
            pressure_var = tk.StringVar(value=self._si_to_display(node.properties.get("pressure", ""), "pressure"))
            pressure_entry = ttk.Entry(container, textvariable=pressure_var, width=20)
            pressure_entry.grid(row=1, column=1, sticky="ew", pady=4)
            entries["pressure"] = pressure_var

            ttk.Label(container, text=f"Flow ({self._unit_label('flow')})").grid(row=2, column=0, sticky="w", pady=4)
            flow_var = tk.StringVar(value=self._si_to_display(node.properties.get("flow", ""), "flow"))
            flow_entry = ttk.Entry(container, textvariable=flow_var, width=20)
            flow_entry.grid(row=2, column=1, sticky="ew", pady=4)
            entries["flow"] = flow_var

            self._sync_boundary_entries(
                condition_var,
                pressure_entry,
                flow_entry,
            )
            condition_var.trace_add(
                "write",
                lambda *_args: self._sync_boundary_entries(
                    condition_var,
                    pressure_entry,
                    flow_entry,
                ),
            )

            if self.scene.physics_mode == "non_isothermal":
                ttk.Separator(container, orient="horizontal").grid(
                    row=3, column=0, columnspan=2, sticky="ew", pady=(6, 2)
                )
                ttk.Label(
                    container, text="— Thermal boundary —", foreground="gray"
                ).grid(row=4, column=0, columnspan=2, pady=(0, 4))

                # Infer default bc_type for old files that lack the property
                default_bc = node.properties.get("thermal_bc_type", "")
                if not default_bc:
                    has_t = node.properties.get("inlet_temperature_c", "").strip()
                    default_bc = "fixed_temperature" if (node.node_type == "source" and has_t) else "zero_gradient"

                thermal_bc_var = tk.StringVar(value=default_bc)
                entries["thermal_bc_type"] = thermal_bc_var

                bc_frame = ttk.Frame(container)
                bc_frame.grid(row=5, column=0, columnspan=2, sticky="w", pady=(0, 4))
                ttk.Radiobutton(bc_frame, text="Zero gradient  (∂T/∂x = 0)",
                                value="zero_gradient", variable=thermal_bc_var).pack(anchor="w")
                ttk.Radiobutton(bc_frame, text="Fixed temperature",
                                value="fixed_temperature", variable=thermal_bc_var).pack(anchor="w")
                ttk.Radiobutton(bc_frame, text="Fixed gradient  (∂T/∂x = g)",
                                value="fixed_gradient", variable=thermal_bc_var).pack(anchor="w")

                # Temperature entry (row 6)
                t_label = ttk.Label(container, text="Temperature (°C)")
                t_label.grid(row=6, column=0, sticky="w", pady=4)
                t_in_var = tk.StringVar(value=node.properties.get("inlet_temperature_c", ""))
                t_entry = ttk.Entry(container, textvariable=t_in_var, width=20)
                t_entry.grid(row=6, column=1, sticky="ew", pady=4)
                entries["inlet_temperature_c"] = t_in_var

                # Gradient entry (row 7)
                g_label = ttk.Label(container, text="Gradient g (°C/m)")
                g_label.grid(row=7, column=0, sticky="w", pady=4)
                g_var = tk.StringVar(value=node.properties.get("thermal_gradient_dc_per_m", "0.0"))
                g_entry = ttk.Entry(container, textvariable=g_var, width=20)
                g_entry.grid(row=7, column=1, sticky="ew", pady=4)
                entries["thermal_gradient_dc_per_m"] = g_var

                def _sync_thermal_bc(*_args: object) -> None:
                    bc = thermal_bc_var.get()
                    if bc == "fixed_temperature":
                        t_label.grid()
                        t_entry.grid()
                        g_label.grid_remove()
                        g_entry.grid_remove()
                    elif bc == "fixed_gradient":
                        t_label.grid_remove()
                        t_entry.grid_remove()
                        g_label.grid()
                        g_entry.grid()
                    else:  # zero_gradient
                        t_label.grid_remove()
                        t_entry.grid_remove()
                        g_label.grid_remove()
                        g_entry.grid_remove()

                thermal_bc_var.trace_add("write", _sync_thermal_bc)
                _sync_thermal_bc()
        else:
            ttk.Label(container, text="Label").grid(row=0, column=0, sticky="w", pady=4)
            label_var = tk.StringVar(value=node.properties.get("label", ""))
            ttk.Entry(container, textvariable=label_var, width=20).grid(
                row=0, column=1, sticky="ew", pady=4
            )
            entries["label"] = label_var

        button_row = ttk.Frame(container)
        button_row.grid(row=10, column=0, columnspan=2, sticky="e", pady=(10, 0))

        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="right", padx=(8, 0))
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_node_properties(node.node_id, entries, dialog),
        ).pack(side="right")

        dialog.update_idletasks()
        dialog.grab_set()
        dialog.focus_set()

    def _open_link_properties_dialog(self, link: CanvasLink) -> None:
        dialog = tk.Toplevel(self.root)
        dialog.title(f"Connection #{link.link_id}")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        container = ttk.Frame(dialog, padding=12)
        container.pack(fill="both", expand=True)
        container.columnconfigure(1, weight=1)
        container.columnconfigure(2, weight=1)

        ttk.Label(container, text="Component Palette").grid(
            row=0, column=0, sticky="w", padx=(0, 12), pady=(0, 8)
        )
        ttk.Label(container, text="Added Components").grid(
            row=0, column=1, sticky="w", pady=(0, 8)
        )
        ttk.Label(container, text="Component Properties").grid(
            row=0, column=2, sticky="w", padx=(12, 0), pady=(0, 8)
        )

        palette = ttk.Frame(container)
        palette.grid(row=1, column=0, sticky="ns")

        components_list = tk.Listbox(container, width=28, height=8)
        components_list.grid(row=1, column=1, sticky="nsew")
        properties_frame = ttk.LabelFrame(container, text="Selected Component", padding=10)
        properties_frame.grid(row=1, column=2, sticky="nsew", padx=(12, 0))

        for component_index, component in enumerate(link.components, start=1):
            components_list.insert("end", self._component_list_label(component, component_index))

        ttk.Button(
            palette,
            text="Pipe",
            command=lambda: self._add_component_to_link(
                link.link_id,
                "pipe",
                components_list,
                properties_frame,
            ),
            width=12,
        ).pack(anchor="w", pady=4)
        ttk.Button(
            palette,
            text="Fitting",
            command=lambda: self._add_component_to_link(
                link.link_id,
                "fitting",
                components_list,
                properties_frame,
            ),
            width=12,
        ).pack(anchor="w", pady=4)
        ttk.Button(
            palette,
            text="Pump",
            command=lambda: self._add_component_to_link(
                link.link_id,
                "pump",
                components_list,
                properties_frame,
            ),
            width=12,
        ).pack(anchor="w", pady=4)
        if self.scene.physics_mode == "non_isothermal":
            ttk.Button(
                palette,
                text="Heat Source",
                command=lambda: self._add_component_to_link(
                    link.link_id,
                    "heat_source",
                    components_list,
                    properties_frame,
                ),
                width=12,
            ).pack(anchor="w", pady=4)

        components_list.bind(
            "<<ListboxSelect>>",
            lambda _event: self._render_link_component_properties(
                link.link_id,
                components_list,
                properties_frame,
            ),
        )

        list_actions = ttk.Frame(container)
        list_actions.grid(row=2, column=1, sticky="w", pady=(4, 0))

        def _refresh_list() -> None:
            updated = self.scene.get_link(link.link_id)
            if updated is None:
                return
            sel = components_list.curselection()
            components_list.delete(0, "end")
            for ci, comp in enumerate(updated.components, start=1):
                components_list.insert("end", self._component_list_label(comp, ci))
            if sel:
                new_idx = min(sel[0], components_list.size() - 1)
                if new_idx >= 0:
                    components_list.selection_set(new_idx)
                    self._render_link_component_properties(link.link_id, components_list, properties_frame)

        def _move_selected(direction: int) -> None:
            sel = components_list.curselection()
            if not sel:
                return
            current_link = self.scene.get_link(link.link_id)
            if current_link is None or sel[0] >= len(current_link.components):
                return
            comp = current_link.components[sel[0]]
            self._push_undo()
            self.scene.move_link_component(link.link_id, comp.component_id, direction)
            _refresh_list()
            new_idx = sel[0] + direction
            if 0 <= new_idx < components_list.size():
                components_list.selection_clear(0, "end")
                components_list.selection_set(new_idx)
                self._render_link_component_properties(link.link_id, components_list, properties_frame)

        def _delete_selected() -> None:
            sel = components_list.curselection()
            if not sel:
                return
            current_link = self.scene.get_link(link.link_id)
            if current_link is None or sel[0] >= len(current_link.components):
                return
            comp = current_link.components[sel[0]]
            self._push_undo()
            self.scene.remove_link_component(link.link_id, comp.component_id)
            self.latest_boundary_results = {}
            self.latest_result = None
            _refresh_list()

        ttk.Button(list_actions, text="↑", width=3, command=lambda: _move_selected(-1)).pack(side="left", padx=(0, 2))
        ttk.Button(list_actions, text="↓", width=3, command=lambda: _move_selected(1)).pack(side="left", padx=(0, 6))
        ttk.Button(list_actions, text="Delete", command=_delete_selected).pack(side="left")

        button_row = ttk.Frame(container)
        button_row.grid(row=3, column=0, columnspan=3, sticky="e", pady=(12, 0))
        ttk.Button(button_row, text="Close", command=dialog.destroy).pack(side="right", padx=(6, 0))
        ttk.Button(
            button_row,
            text="Temperature Profile",
            command=lambda: self._show_link_temperature_profile(link),
            state="normal" if self._has_temperature_results() else "disabled",
        ).pack(side="right", padx=(0, 4))
        ttk.Button(
            button_row,
            text="Pressure Profile",
            command=lambda: self._show_link_pressure_profile(link),
            state="normal" if self.latest_boundary_results else "disabled",
        ).pack(side="right")

        dialog.update_idletasks()
        dialog.grab_set()
        dialog.focus_set()

    def _link_node_sequence(self, link) -> list[int]:
        """Return [start, internal_1, ..., end] node IDs matching solver result indices."""
        next_id = max(node.node_id for node in self.scene.nodes) + 1
        for scene_link in self.scene.links:
            node_ids: list[int] = [scene_link.start_node_id]
            for comp_idx, _comp in enumerate(scene_link.components):
                is_last = comp_idx == len(scene_link.components) - 1
                end = scene_link.end_node_id if is_last else next_id
                if not is_last:
                    next_id += 1
                node_ids.append(end)
            if scene_link.link_id == link.link_id:
                return node_ids
        return [link.start_node_id, link.end_node_id]

    def _show_link_pressure_profile(self, link) -> None:
        win = tk.Toplevel(self.root)
        win.title(f"Pressure Profile — Connection #{link.link_id}")
        win.geometry("760x480")

        frame = ttk.Frame(win, padding=10)
        frame.pack(fill="both", expand=True)

        plot_canvas = tk.Canvas(
            frame,
            background=self._t["plot_bg"],
            highlightthickness=1,
            highlightbackground=self._t["canvas_hl"],
        )
        plot_canvas.pack(fill="both", expand=True)
        plot_canvas.bind(
            "<Configure>",
            lambda _event: self._draw_pressure_profile_plot(plot_canvas, link),
        )

        btn_row = ttk.Frame(frame)
        btn_row.pack(fill="x", pady=(8, 0))
        ttk.Button(btn_row, text="Edit Components…",
                   command=lambda: self._open_link_properties_dialog(link)).pack(side="left")
        ttk.Button(btn_row, text="Close", command=win.destroy).pack(side="right")

        win.transient(self.root)
        win.focus_set()
        win.update_idletasks()
        self._draw_pressure_profile_plot(plot_canvas, link)

    def _draw_pressure_profile_plot(self, canvas: tk.Canvas, link) -> None:
        canvas.delete("all")
        W = int(canvas.winfo_width() or 720)
        H = int(canvas.winfo_height() or 380)

        ML, MR, MT, MB = 75, 20, 20, 50
        plot_w = W - ML - MR
        plot_h = H - MT - MB

        if plot_w < 100 or plot_h < 80:
            return

        node_ids = self._link_node_sequence(link)
        results = self.latest_boundary_results

        if not results or not all(nid in results for nid in node_ids):
            canvas.create_text(
                W // 2, H // 2,
                text="Run the solver first to see the pressure profile.",
                fill=self._t["plot_muted"],
                font=("TkDefaultFont", 11),
                anchor="center",
            )
            return

        p_unit = self._unit_label("pressure")
        _, p_factor = self._unit_quantities()["pressure"]

        pressures_disp = [results[nid]["pressure_pa"] * p_factor for nid in node_ids]

        # cumulative distance: pipes contribute length_m; fittings and pumps are point elements
        cum_dist: list[float] = [0.0]
        for comp in link.components:
            if comp.component_type == "pipe":
                try:
                    length = float(comp.properties.get("length_m", "0") or "0")
                except ValueError:
                    length = 0.0
                cum_dist.append(cum_dist[-1] + length)
            else:
                cum_dist.append(cum_dist[-1])

        x_max = max(cum_dist) if max(cum_dist) > 0 else 1.0
        x_min = 0.0

        p_min = min(pressures_disp)
        p_max = max(pressures_disp)
        p_span = p_max - p_min
        if p_span < 1e-10:
            p_span = max(abs(p_max) * 0.1, 1.0)
        pad = p_span * 0.08
        p_lo = p_min - pad
        p_hi = p_max + pad
        p_span = p_hi - p_lo

        def _cx(x): return ML + (x - x_min) / (x_max - x_min) * plot_w
        def _cy(p): return MT + plot_h - (p - p_lo) / p_span * plot_h

        # plot border
        canvas.create_rectangle(ML, MT, ML + plot_w, MT + plot_h,
                                 outline=self._t["plot_axis"], width=1)

        # horizontal grid + y-tick labels
        for i in range(6):
            p_val = p_lo + i / 5 * p_span
            y = _cy(p_val)
            canvas.create_line(ML, y, ML + plot_w, y, fill=self._t["plot_grid"])
            canvas.create_text(ML - 6, y, text=f"{p_val:.4g}", anchor="e",
                                font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        # vertical grid + x-tick labels
        for i in range(6):
            x_val = x_min + i / 5 * (x_max - x_min)
            x = _cx(x_val)
            canvas.create_line(x, MT, x, MT + plot_h, fill=self._t["plot_grid"])
            canvas.create_text(x, MT + plot_h + 4, text=f"{x_val:.4g}", anchor="n",
                                font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        # axes
        canvas.create_line(ML, MT, ML, MT + plot_h, fill=self._t["plot_axis"], width=1.5)
        canvas.create_line(ML, MT + plot_h, ML + plot_w, MT + plot_h,
                            fill=self._t["plot_axis"], width=1.5)

        # axis labels
        canvas.create_text(ML + plot_w // 2, MT + plot_h + 36,
                            text="Cumulative distance (m)", anchor="s",
                            font=("TkDefaultFont", 9), fill=self._t["plot_text"])
        canvas.create_text(12, MT + plot_h // 2,
                            text=f"Pressure ({p_unit})", anchor="center",
                            font=("TkDefaultFont", 9), fill=self._t["plot_text"],
                            angle=90)

        pipe_color = self._t["plot_laminar"]
        fitting_color = "#e69c00"
        pump_color = self._t["plot_turbulent"]
        heater_color = "#c0397d"

        legend_seen: dict[str, tuple[str, tuple]] = {}

        for i, comp in enumerate(link.components):
            x0c, y0c = _cx(cum_dist[i]), _cy(pressures_disp[i])
            x1c, y1c = _cx(cum_dist[i + 1]), _cy(pressures_disp[i + 1])

            if comp.component_type == "pipe":
                canvas.create_line(x0c, y0c, x1c, y1c, fill=pipe_color, width=2)
                legend_seen.setdefault("Pipe", (pipe_color, ()))
            elif comp.component_type == "fitting":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=fitting_color, width=2, dash=(6, 3))
                legend_seen.setdefault("Fitting", (fitting_color, (6, 3)))
            elif comp.component_type == "pump":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=pump_color, width=2, dash=(6, 3))
                legend_seen.setdefault("Pump", (pump_color, (6, 3)))
            elif comp.component_type == "heat_source":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=heater_color, width=2, dash=(4, 4))
                legend_seen.setdefault("Heat Source", (heater_color, (4, 4)))

        # markers: filled circle at start/end, × at intermediate device boundaries
        n_nodes = len(cum_dist)
        for idx, (x_m, p_d) in enumerate(zip(cum_dist, pressures_disp)):
            xc, yc = _cx(x_m), _cy(p_d)
            if idx == 0 or idx == n_nodes - 1:
                r = 4
                canvas.create_oval(xc - r, yc - r, xc + r, yc + r,
                                   fill=self._t["plot_text"],
                                   outline=self._t["plot_axis"], width=1)
            else:
                r = 5
                canvas.create_line(xc - r, yc - r, xc + r, yc + r,
                                   fill=self._t["plot_axis"], width=2)
                canvas.create_line(xc + r, yc - r, xc - r, yc + r,
                                   fill=self._t["plot_axis"], width=2)

        # segment type labels (midpoint of each component's x-range)
        for i, comp in enumerate(link.components):
            if comp.component_type == "pipe":
                mid_x = _cx(0.5 * (cum_dist[i] + cum_dist[i + 1]))
                mid_y = _cy(0.5 * (pressures_disp[i] + pressures_disp[i + 1]))
                canvas.create_text(mid_x, mid_y - 9, text="pipe",
                                   anchor="s", font=("TkDefaultFont", 7),
                                   fill=pipe_color)
            else:
                lbl_map = {"pump": "pump", "fitting": "fitting", "heat_source": "heater"}
                clr_map = {"pump": pump_color, "fitting": fitting_color, "heat_source": heater_color}
                lbl = lbl_map.get(comp.component_type, comp.component_type)
                clr = clr_map.get(comp.component_type, fitting_color)
                mid_y = _cy(0.5 * (pressures_disp[i] + pressures_disp[i + 1]))
                lbl_x = _cx(cum_dist[i]) + 7
                canvas.create_text(lbl_x, mid_y, text=lbl,
                                   anchor="w", font=("TkDefaultFont", 7),
                                   fill=clr)

        # legend (top-left of plot area)
        lx, ly = ML + 10, MT + 10
        for label, (color, dash) in legend_seen.items():
            canvas.create_line(lx, ly + 5, lx + 22, ly + 5, fill=color, width=2, dash=dash)
            canvas.create_text(lx + 27, ly + 5, text=label, anchor="w",
                                font=("TkDefaultFont", 9), fill=self._t["plot_text"])
            ly += 20

    def _has_temperature_results(self) -> bool:
        return any(
            "temperature_c" in v for v in self.latest_boundary_results.values()
        )

    def _show_link_temperature_profile(self, link) -> None:
        win = tk.Toplevel(self.root)
        win.title(f"Temperature Profile — Connection #{link.link_id}")
        win.geometry("760x480")

        frame = ttk.Frame(win, padding=10)
        frame.pack(fill="both", expand=True)

        plot_canvas = tk.Canvas(
            frame,
            background=self._t["plot_bg"],
            highlightthickness=1,
            highlightbackground=self._t["canvas_hl"],
        )
        plot_canvas.pack(fill="both", expand=True)
        plot_canvas.bind(
            "<Configure>",
            lambda _event: self._draw_temperature_profile_plot(plot_canvas, link),
        )

        btn_row = ttk.Frame(frame)
        btn_row.pack(fill="x", pady=(8, 0))
        ttk.Button(btn_row, text="Edit Components…",
                   command=lambda: self._open_link_properties_dialog(link)).pack(side="left")
        ttk.Button(btn_row, text="Close", command=win.destroy).pack(side="right")

        win.transient(self.root)
        win.focus_set()
        win.update_idletasks()
        self._draw_temperature_profile_plot(plot_canvas, link)

    def _draw_temperature_profile_plot(self, canvas: tk.Canvas, link) -> None:
        canvas.delete("all")
        W = int(canvas.winfo_width() or 720)
        H = int(canvas.winfo_height() or 380)

        ML, MR, MT, MB = 75, 20, 20, 50
        plot_w = W - ML - MR
        plot_h = H - MT - MB

        if plot_w < 100 or plot_h < 80:
            return

        node_ids = self._link_node_sequence(link)
        results = self.latest_boundary_results

        if not results or not all(
            nid in results and "temperature_c" in results[nid] for nid in node_ids
        ):
            canvas.create_text(
                W // 2, H // 2,
                text="Run a non-isothermal simulation first to see the temperature profile.",
                fill=self._t["plot_muted"],
                font=("TkDefaultFont", 11),
                anchor="center",
            )
            return

        temps_disp = [results[nid]["temperature_c"] for nid in node_ids]

        cum_dist: list[float] = [0.0]
        for comp in link.components:
            if comp.component_type == "pipe":
                try:
                    length = float(comp.properties.get("length_m", "0") or "0")
                except ValueError:
                    length = 0.0
                cum_dist.append(cum_dist[-1] + length)
            else:
                cum_dist.append(cum_dist[-1])

        x_max = max(cum_dist) if max(cum_dist) > 0 else 1.0
        x_min = 0.0

        t_min = min(temps_disp)
        t_max = max(temps_disp)
        t_span = t_max - t_min
        if t_span < 1e-10:
            t_span = max(abs(t_max) * 0.1, 1.0)
        pad = t_span * 0.08
        t_lo = t_min - pad
        t_hi = t_max + pad
        t_span = t_hi - t_lo

        def _cx(x): return ML + (x - x_min) / (x_max - x_min) * plot_w
        def _cy(t): return MT + plot_h - (t - t_lo) / t_span * plot_h

        canvas.create_rectangle(ML, MT, ML + plot_w, MT + plot_h,
                                 outline=self._t["plot_axis"], width=1)

        for i in range(6):
            t_val = t_lo + i / 5 * t_span
            y = _cy(t_val)
            canvas.create_line(ML, y, ML + plot_w, y, fill=self._t["plot_grid"])
            canvas.create_text(ML - 6, y, text=f"{t_val:.4g}", anchor="e",
                                font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        for i in range(6):
            x_val = x_min + i / 5 * (x_max - x_min)
            x = _cx(x_val)
            canvas.create_line(x, MT, x, MT + plot_h, fill=self._t["plot_grid"])
            canvas.create_text(x, MT + plot_h + 4, text=f"{x_val:.4g}", anchor="n",
                                font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        canvas.create_line(ML, MT, ML, MT + plot_h, fill=self._t["plot_axis"], width=1.5)
        canvas.create_line(ML, MT + plot_h, ML + plot_w, MT + plot_h,
                            fill=self._t["plot_axis"], width=1.5)

        canvas.create_text(ML + plot_w // 2, MT + plot_h + 36,
                            text="Cumulative distance (m)", anchor="s",
                            font=("TkDefaultFont", 9), fill=self._t["plot_text"])
        canvas.create_text(12, MT + plot_h // 2,
                            text="Temperature (°C)", anchor="center",
                            font=("TkDefaultFont", 9), fill=self._t["plot_text"],
                            angle=90)

        temp_color = self._t["plot_temperature"]
        fitting_color = "#e69c00"
        pump_color = self._t["plot_turbulent"]
        heater_color = "#c0397d"

        legend_seen: dict[str, tuple[str, tuple]] = {}

        for i, comp in enumerate(link.components):
            x0c, y0c = _cx(cum_dist[i]), _cy(temps_disp[i])
            x1c, y1c = _cx(cum_dist[i + 1]), _cy(temps_disp[i + 1])

            if comp.component_type == "pipe":
                canvas.create_line(x0c, y0c, x1c, y1c, fill=temp_color, width=2)
                legend_seen.setdefault("Pipe", (temp_color, ()))
            elif comp.component_type == "fitting":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=fitting_color, width=2, dash=(6, 3))
                legend_seen.setdefault("Fitting", (fitting_color, (6, 3)))
            elif comp.component_type == "pump":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=pump_color, width=2, dash=(6, 3))
                legend_seen.setdefault("Pump", (pump_color, (6, 3)))
            elif comp.component_type == "heat_source":
                canvas.create_line(x0c, y0c, x0c, y1c, fill=heater_color, width=2, dash=(4, 4))
                legend_seen.setdefault("Heat Source", (heater_color, (4, 4)))

        n_nodes = len(cum_dist)
        for idx, (x_m, t_d) in enumerate(zip(cum_dist, temps_disp)):
            xc, yc = _cx(x_m), _cy(t_d)
            if idx == 0 or idx == n_nodes - 1:
                r = 4
                canvas.create_oval(xc - r, yc - r, xc + r, yc + r,
                                   fill=self._t["plot_text"],
                                   outline=self._t["plot_axis"], width=1)
            else:
                r = 5
                canvas.create_line(xc - r, yc - r, xc + r, yc + r,
                                   fill=self._t["plot_axis"], width=2)
                canvas.create_line(xc + r, yc - r, xc - r, yc + r,
                                   fill=self._t["plot_axis"], width=2)

        for i, comp in enumerate(link.components):
            if comp.component_type == "pipe":
                mid_x = _cx(0.5 * (cum_dist[i] + cum_dist[i + 1]))
                mid_y = _cy(0.5 * (temps_disp[i] + temps_disp[i + 1]))
                canvas.create_text(mid_x, mid_y - 9, text="pipe",
                                   anchor="s", font=("TkDefaultFont", 7),
                                   fill=temp_color)
            else:
                lbl_map = {"pump": "pump", "fitting": "fitting", "heat_source": "heater"}
                clr_map = {"pump": pump_color, "fitting": fitting_color, "heat_source": heater_color}
                lbl = lbl_map.get(comp.component_type, comp.component_type)
                clr = clr_map.get(comp.component_type, fitting_color)
                mid_y = _cy(0.5 * (temps_disp[i] + temps_disp[i + 1]))
                lbl_x = _cx(cum_dist[i]) + 7
                canvas.create_text(lbl_x, mid_y, text=lbl,
                                   anchor="w", font=("TkDefaultFont", 7),
                                   fill=clr)

        lx, ly = ML + 10, MT + 10
        for label, (color, dash) in legend_seen.items():
            canvas.create_line(lx, ly + 5, lx + 22, ly + 5, fill=color, width=2, dash=dash)
            canvas.create_text(lx + 27, ly + 5, text=label, anchor="w",
                                font=("TkDefaultFont", 9), fill=self._t["plot_text"])
            ly += 20

    def _add_component_to_link(
        self,
        link_id: int,
        component_type: str,
        components_list: tk.Listbox,
        properties_frame: ttk.LabelFrame,
    ) -> None:
        updated_link = self.scene.add_link_component(link_id, component_type)
        components_list.delete(0, "end")
        for component_index, component in enumerate(updated_link.components, start=1):
            components_list.insert("end", self._component_list_label(component, component_index))
        components_list.selection_clear(0, "end")
        components_list.selection_set("end")
        self._render_link_component_properties(link_id, components_list, properties_frame)
        self.status_var.set(
            f"Added {component_type} to connection #{updated_link.link_id}."
        )

    def _render_link_component_properties(
        self,
        link_id: int,
        components_list: tk.Listbox,
        properties_frame: ttk.LabelFrame,
    ) -> None:
        for child in properties_frame.winfo_children():
            child.destroy()

        selected = components_list.curselection()
        if not selected:
            ttk.Label(properties_frame, text="Select a component to edit its data.").pack(
                anchor="w"
            )
            return

        link = self.scene.get_link(link_id)
        if link is None:
            return

        component = link.components[selected[0]]
        component_index = selected[0] + 1
        entries: dict[str, tk.StringVar] = {}

        ttk.Label(
            properties_frame,
            text=self._component_list_label(component, component_index),
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))

        if component.component_type == "fitting":
            self._render_fitting_component_properties(
                component,
                entries,
                properties_frame,
                link_id,
                components_list,
            )
            return
        if component.component_type == "pump":
            self._render_pump_component_properties(
                component,
                properties_frame,
                link_id,
                components_list,
            )
            return
        if component.component_type == "heat_source":
            self._render_heat_source_component_properties(
                component,
                properties_frame,
                link_id,
                components_list,
            )
            return

        row = 1
        for key, value in component.properties.items():
            ttk.Label(properties_frame, text=self._pretty_field_name(key)).grid(
                row=row, column=0, sticky="w", pady=4
            )
            var = tk.StringVar(value=self._field_from_si(key, value))
            ttk.Entry(properties_frame, textvariable=var, width=18).grid(
                row=row, column=1, sticky="ew", pady=4
            )
            entries[key] = var
            row += 1

        if self.scene.physics_mode == "non_isothermal":
            ttk.Separator(properties_frame, orient="horizontal").grid(
                row=row, column=0, columnspan=2, sticky="ew", pady=(4, 2)
            )
            row += 1
            ttk.Label(
                properties_frame, text="— Heat transfer —", foreground="gray"
            ).grid(row=row, column=0, columnspan=2, pady=(0, 2))
            row += 1
            thermal_fields = [
                ("heat_transfer_coefficient_w_per_m2k", "U — Heat transfer coeff. (W/m²K)"),
                ("ambient_temperature_c", "T_amb — Ambient temperature (°C)"),
                ("n_thermal_segments", "Thermal segments"),
            ]
            defaults = {
                "heat_transfer_coefficient_w_per_m2k": "0.0",
                "ambient_temperature_c": "20.0",
                "n_thermal_segments": "10",
            }
            for key, label in thermal_fields:
                ttk.Label(properties_frame, text=label).grid(
                    row=row, column=0, sticky="w", pady=4
                )
                var = tk.StringVar(
                    value=component.properties.get(key, defaults[key])
                )
                ttk.Entry(properties_frame, textvariable=var, width=18).grid(
                    row=row, column=1, sticky="ew", pady=4
                )
                entries[key] = var
                row += 1

        button_row = ttk.Frame(properties_frame)
        button_row.grid(row=row, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_link_component_properties(
                link_id,
                component.component_id,
                entries,
                components_list,
                properties_frame,
            ),
        ).pack(side="right")

    def _render_fitting_component_properties(
        self,
        component: CanvasLinkComponent,
        entries: dict[str, tk.StringVar],
        properties_frame: ttk.LabelFrame,
        link_id: int,
        components_list: tk.Listbox,
    ) -> None:
        diameter_var = tk.StringVar(value=self._field_from_si("diameter_m", component.properties.get("diameter_m", "")))
        mode_var = tk.StringVar(value=component.properties.get("fitting_mode", "manual"))
        preset_var = tk.StringVar(
            value=component.properties.get("fitting_preset", "regular_90_flanged")
        )
        loss_var = tk.StringVar(value=component.properties.get("loss_coefficient", "1.5"))

        entries["diameter_m"] = diameter_var
        entries["fitting_mode"] = mode_var
        entries["fitting_preset"] = preset_var
        entries["loss_coefficient"] = loss_var

        ttk.Label(properties_frame, text=f"Diameter ({self._unit_label('diameter')})").grid(
            row=1, column=0, sticky="w", pady=4
        )
        ttk.Entry(properties_frame, textvariable=diameter_var, width=18).grid(
            row=1, column=1, sticky="ew", pady=4
        )

        ttk.Label(properties_frame, text="Fitting Mode").grid(
            row=2, column=0, sticky="w", pady=4
        )
        mode_box = ttk.Combobox(
            properties_frame,
            textvariable=mode_var,
            state="readonly",
            values=tuple(self.FITTING_MODE_LIBRARY.keys()),
            width=18,
        )
        mode_box.grid(row=2, column=1, sticky="ew", pady=4)

        mode_name_var = tk.StringVar(
            value=self.FITTING_MODE_LIBRARY[mode_var.get()]["name"]
        )
        ttk.Label(properties_frame, text="Selected Mode").grid(
            row=3, column=0, sticky="w", pady=4
        )
        ttk.Label(
            properties_frame,
            textvariable=mode_name_var,
            relief="groove",
            padding=6,
            width=20,
        ).grid(row=3, column=1, sticky="ew", pady=4)

        ttk.Label(properties_frame, text="Fitting Preset").grid(
            row=4, column=0, sticky="w", pady=4
        )
        preset_box = ttk.Combobox(
            properties_frame,
            textvariable=preset_var,
            state="readonly",
            values=tuple(self.FITTING_PRESET_LIBRARY.keys()),
            width=18,
        )
        preset_box.grid(row=4, column=1, sticky="ew", pady=4)

        preset_name_var = tk.StringVar(
            value=self.FITTING_PRESET_LIBRARY[preset_var.get()]["name"]
        )
        ttk.Label(properties_frame, text="Selected Preset").grid(
            row=5, column=0, sticky="w", pady=4
        )
        ttk.Label(
            properties_frame,
            textvariable=preset_name_var,
            relief="groove",
            padding=6,
            width=20,
        ).grid(row=5, column=1, sticky="ew", pady=4)

        ttk.Label(properties_frame, text="Loss Coefficient K").grid(
            row=6, column=0, sticky="w", pady=4
        )
        loss_entry = ttk.Entry(properties_frame, textvariable=loss_var, width=18)
        loss_entry.grid(row=6, column=1, sticky="ew", pady=4)

        def sync_fitting_controls(*_args: object) -> None:
            mode_name_var.set(self.FITTING_MODE_LIBRARY[mode_var.get()]["name"])
            preset_name_var.set(self.FITTING_PRESET_LIBRARY[preset_var.get()]["name"])
            if mode_var.get() == "preset":
                preset_box.configure(state="readonly")
                loss_entry.configure(state="disabled")
                preset_loss = self.FITTING_PRESET_LIBRARY[preset_var.get()]["loss_coefficient"]
                if math.isinf(preset_loss):
                    loss_var.set("inf")
                else:
                    loss_var.set(f"{preset_loss:g}")
            else:
                preset_box.configure(state="disabled")
                loss_entry.configure(state="normal")

        mode_var.trace_add("write", sync_fitting_controls)
        preset_var.trace_add("write", sync_fitting_controls)
        sync_fitting_controls()

        button_row = ttk.Frame(properties_frame)
        button_row.grid(row=7, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_link_component_properties(
                link_id,
                component.component_id,
                entries,
                components_list,
                properties_frame,
            ),
        ).pack(side="right")

    def _render_pump_component_properties(
        self,
        component: CanvasLinkComponent,
        properties_frame: ttk.LabelFrame,
        link_id: int,
        components_list: tk.Listbox,
    ) -> None:
        diameter_var = tk.StringVar(value=self._field_from_si("diameter_m", component.properties.get("diameter_m", "")))

        ttk.Label(properties_frame, text=f"Diameter ({self._unit_label('diameter')})").grid(
            row=1, column=0, sticky="w", pady=4
        )
        ttk.Entry(properties_frame, textvariable=diameter_var, width=18).grid(
            row=1, column=1, sticky="ew", pady=4
        )

        ttk.Label(properties_frame, text="Q-Head Table").grid(
            row=2, column=0, sticky="nw", pady=4
        )
        curve_text = tk.Text(properties_frame, width=24, height=7, wrap="none")
        curve_text.grid(row=2, column=1, sticky="ew", pady=4)
        existing_curve_text = component.properties.get("curve_points_q_head", "").strip()
        if existing_curve_text:
            curve_text.insert("1.0", existing_curve_text)

        ttk.Label(
            properties_frame,
            text="One pair per line: Q (m^3/h), Head (m)",
        ).grid(row=3, column=0, columnspan=2, sticky="w", pady=(0, 6))

        button_row = ttk.Frame(properties_frame)
        button_row.grid(row=4, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_pump_component_properties(
                link_id,
                component.component_id,
                diameter_var,
                curve_text,
                components_list,
                properties_frame,
            ),
        ).pack(side="right")

    def _save_link_component_properties(
        self,
        link_id: int,
        component_id: int,
        entries: dict[str, tk.StringVar],
        components_list: tk.Listbox,
        properties_frame: ttk.LabelFrame,
    ) -> None:
        properties = {key: self._field_to_si(key, value.get().strip()) for key, value in entries.items()}
        updated_link = self.scene.update_link_component_properties(
            link_id,
            component_id,
            properties,
        )
        selection = components_list.curselection()
        components_list.delete(0, "end")
        for component_index, component in enumerate(updated_link.components, start=1):
            components_list.insert("end", self._component_list_label(component, component_index))
        if selection:
            components_list.selection_set(selection[0])
        self._render_link_component_properties(link_id, components_list, properties_frame)
        self.status_var.set(f"Updated component #{component_id} in connection #{link_id}.")

    def _save_pump_component_properties(
        self,
        link_id: int,
        component_id: int,
        diameter_var: tk.StringVar,
        curve_text: tk.Text,
        components_list: tk.Listbox,
        properties_frame: ttk.LabelFrame,
    ) -> None:
        properties = {
            "diameter_m": self._field_to_si("diameter_m", diameter_var.get().strip()),
            "curve_points_q_head": curve_text.get("1.0", "end").strip(),
        }
        updated_link = self.scene.update_link_component_properties(
            link_id,
            component_id,
            properties,
        )
        selection = components_list.curselection()
        components_list.delete(0, "end")
        for component_index, component in enumerate(updated_link.components, start=1):
            components_list.insert("end", self._component_list_label(component, component_index))
        if selection:
            components_list.selection_set(selection[0])
        self._render_link_component_properties(link_id, components_list, properties_frame)
        self.status_var.set(f"Updated component #{component_id} in connection #{link_id}.")

    def _render_heat_source_component_properties(
        self,
        component: CanvasLinkComponent,
        properties_frame: ttk.LabelFrame,
        link_id: int,
        components_list: tk.Listbox,
    ) -> None:
        diameter_var = tk.StringVar(value=self._field_from_si("diameter_m", component.properties.get("diameter_m", "")))
        power_var = tk.StringVar(value=component.properties.get("power_w", "0.0"))
        mode_var = tk.StringVar(value=component.properties.get("pressure_drop_mode", "rated"))
        dp_var = tk.StringVar(value=component.properties.get("pressure_drop_pa", "0.0"))
        rated_mdot_var = tk.StringVar(value=component.properties.get("rated_mass_flow_kg_per_s", "1.0"))
        n_segs_var = tk.StringVar(value=component.properties.get("n_thermal_segments", "10"))

        row = 1
        ttk.Label(properties_frame, text=f"Diameter ({self._unit_label('diameter')})").grid(
            row=row, column=0, sticky="w", pady=4
        )
        ttk.Entry(properties_frame, textvariable=diameter_var, width=18).grid(
            row=row, column=1, sticky="ew", pady=4
        )
        row += 1

        ttk.Label(properties_frame, text="Power (W)").grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(properties_frame, textvariable=power_var, width=18).grid(
            row=row, column=1, sticky="ew", pady=4
        )
        row += 1

        ttk.Label(properties_frame, text="ΔP mode").grid(row=row, column=0, sticky="w", pady=4)
        mode_box = ttk.Combobox(
            properties_frame,
            textvariable=mode_var,
            values=["rated", "fixed"],
            state="readonly",
            width=16,
        )
        mode_box.grid(row=row, column=1, sticky="ew", pady=4)
        row += 1

        ttk.Label(properties_frame, text="Pressure drop (Pa)").grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(properties_frame, textvariable=dp_var, width=18).grid(
            row=row, column=1, sticky="ew", pady=4
        )
        row += 1

        rated_label = ttk.Label(properties_frame, text="Rated flow (kg/s)")
        rated_label.grid(row=row, column=0, sticky="w", pady=4)
        rated_entry = ttk.Entry(properties_frame, textvariable=rated_mdot_var, width=18)
        rated_entry.grid(row=row, column=1, sticky="ew", pady=4)
        row += 1

        def _sync_mode(*_args: object) -> None:
            if mode_var.get() == "rated":
                rated_label.grid()
                rated_entry.grid()
            else:
                rated_label.grid_remove()
                rated_entry.grid_remove()

        mode_var.trace_add("write", _sync_mode)
        _sync_mode()

        if self.scene.physics_mode == "non_isothermal":
            ttk.Separator(properties_frame, orient="horizontal").grid(
                row=row, column=0, columnspan=2, sticky="ew", pady=(4, 2)
            )
            row += 1
            ttk.Label(
                properties_frame, text="— Heat transfer —", foreground="gray"
            ).grid(row=row, column=0, columnspan=2, pady=(0, 2))
            row += 1
            ttk.Label(properties_frame, text="Thermal segments").grid(
                row=row, column=0, sticky="w", pady=4
            )
            ttk.Entry(properties_frame, textvariable=n_segs_var, width=18).grid(
                row=row, column=1, sticky="ew", pady=4
            )
            row += 1

        button_row = ttk.Frame(properties_frame)
        button_row.grid(row=row, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(
            button_row,
            text="Save",
            command=lambda: self._save_heat_source_component_properties(
                link_id,
                component.component_id,
                diameter_var,
                power_var,
                mode_var,
                dp_var,
                rated_mdot_var,
                n_segs_var,
                components_list,
                properties_frame,
            ),
        ).pack(side="right")

    def _save_heat_source_component_properties(
        self,
        link_id: int,
        component_id: int,
        diameter_var: tk.StringVar,
        power_var: tk.StringVar,
        mode_var: tk.StringVar,
        dp_var: tk.StringVar,
        rated_mdot_var: tk.StringVar,
        n_segs_var: tk.StringVar,
        components_list: tk.Listbox,
        properties_frame: ttk.LabelFrame,
    ) -> None:
        properties = {
            "diameter_m": self._field_to_si("diameter_m", diameter_var.get().strip()),
            "power_w": power_var.get().strip(),
            "pressure_drop_mode": mode_var.get().strip(),
            "pressure_drop_pa": dp_var.get().strip(),
            "rated_mass_flow_kg_per_s": rated_mdot_var.get().strip(),
            "n_thermal_segments": n_segs_var.get().strip(),
        }
        updated_link = self.scene.update_link_component_properties(
            link_id,
            component_id,
            properties,
        )
        selection = components_list.curselection()
        components_list.delete(0, "end")
        for component_index, component in enumerate(updated_link.components, start=1):
            components_list.insert("end", self._component_list_label(component, component_index))
        if selection:
            components_list.selection_set(selection[0])
        self._render_link_component_properties(link_id, components_list, properties_frame)
        self.status_var.set(f"Updated component #{component_id} in connection #{link_id}.")

    @staticmethod
    def _sync_boundary_entries(
        condition_var: tk.StringVar,
        pressure_entry: ttk.Entry,
        flow_entry: ttk.Entry,
    ) -> None:
        if condition_var.get() == "pressure":
            pressure_entry.state(["!disabled"])
            flow_entry.state(["disabled"])
        else:
            pressure_entry.state(["disabled"])
            flow_entry.state(["!disabled"])

    def _save_node_properties(
        self,
        node_id: int,
        entries: dict[str, tk.StringVar],
        dialog: tk.Toplevel,
    ) -> None:
        properties: dict[str, str] = {}
        for key, var in entries.items():
            val = var.get().strip()
            if key == "pressure":
                properties[key] = self._display_to_si(val, "pressure")
            elif key == "flow":
                properties[key] = self._display_to_si(val, "flow")
            else:
                properties[key] = val
        updated_node = self.scene.update_node_properties(node_id, properties)
        self._redraw_scene()
        self.status_var.set(
            f"Updated properties for {updated_node.node_type} #{updated_node.node_id}."
        )
        dialog.destroy()

    def _node_fill(self, node_type: str) -> str:
        return self._t[f"node_{node_type}"]

    @staticmethod
    def _node_label(node: CanvasNode) -> str:
        labels = {
            "source": "S",
            "sink": "K",
            "junction": "J",
        }
        return labels[node.node_type]

    @staticmethod
    def _component_list_label(component: CanvasLinkComponent, display_index: int) -> str:
        if component.component_type == "fitting":
            fitting_mode = component.properties.get("fitting_mode", "manual")
            if fitting_mode == "preset":
                preset_key = component.properties.get("fitting_preset", "")
                preset = FITTING_PRESET_LIBRARY.get(preset_key)
                if preset is not None:
                    return f"{preset['name']} #{display_index}"
            return f"Manual fitting #{display_index}"
        if component.component_type == "pump":
            return f"Pump #{display_index}"
        if component.component_type == "heat_source":
            return f"Heat Source #{display_index}"
        return f"{component.component_type.capitalize()} #{display_index}"

    def _pretty_field_name(self, field_name: str) -> str:
        quantity = self._FIELD_QUANTITY.get(field_name)
        if quantity:
            unit = self._unit_label(quantity)
            base = {
                "length_m": "Length",
                "diameter_m": "Diameter",
                "roughness_m": "Roughness",
                "height_change_m": "Height change",
            }[field_name]
            return f"{base} ({unit})"
        if field_name == "hazen_williams_c":
            return "Hazen-Williams C (−)"
        return field_name.replace("_", " ").capitalize()

    def _node_summary_text(self, node: CanvasNode) -> str:
        lines: list[str] = []
        if node.node_type in {"source", "sink"}:
            condition_type = node.properties.get("condition_type", "pressure")
            _, p_factor = self._unit_quantities()["pressure"]
            _, q_factor = self._unit_quantities()["flow"]
            p_unit = self._unit_label("pressure")
            q_unit = self._unit_label("flow")
            if condition_type == "pressure":
                pressure_si = node.properties.get("pressure", "").strip()
                try:
                    lines.append(f"P={self._fmt(float(pressure_si) * p_factor)} {p_unit}")
                except ValueError:
                    pass
            else:
                flow_si = node.properties.get("flow", "").strip()
                try:
                    lines.append(f"Q={self._fmt(float(flow_si) * q_factor)} {q_unit}")
                except ValueError:
                    pass

            result_data = self.latest_boundary_results.get(node.node_id)
            if result_data is not None:
                if condition_type == "pressure" and "flow_kg_per_s" in result_data:
                    lines.append(f"Q={self._fmt(result_data['flow_kg_per_s'] * q_factor)} {q_unit}")
                elif condition_type == "flow" and "pressure_pa" in result_data:
                    lines.append(f"P={self._fmt(result_data['pressure_pa'] * p_factor)} {p_unit}")
                if "temperature_c" in result_data:
                    lines.append(f"T={result_data['temperature_c']:.1f} °C")
            return "\n".join(lines)

        result_data = self.latest_boundary_results.get(node.node_id)
        if result_data:
            _, p_factor = self._unit_quantities()["pressure"]
            p_unit = self._unit_label("pressure")
            summary_lines = []
            if "pressure_pa" in result_data:
                summary_lines.append(f"P={self._fmt(result_data['pressure_pa'] * p_factor)} {p_unit}")
            if "temperature_c" in result_data:
                summary_lines.append(f"T={result_data['temperature_c']:.1f} °C")
            return "\n".join(summary_lines)
        return ""

    def _build_network_case_from_scene(self) -> NetworkCase:
        return build_network_case_from_scene(self.scene)

    def _build_boundary_results(self, case: NetworkCase, result) -> dict[int, dict[str, float]]:
        boundary_results: dict[int, dict[str, float]] = {}
        for node_id, pressure in result.node_pressures_pa.items():
            boundary_results[node_id] = {"pressure_pa": pressure}

        for node_id, temp_c in result.node_temperatures_c.items():
            boundary_results.setdefault(node_id, {})["temperature_c"] = temp_c

        for component, component_result in zip(case.components, result.component_flows):
            start_entry = boundary_results.setdefault(component.start_node, {"pressure_pa": result.node_pressures_pa.get(component.start_node, 0.0)})
            end_entry = boundary_results.setdefault(component.end_node, {"pressure_pa": result.node_pressures_pa.get(component.end_node, 0.0)})
            start_entry["flow_kg_per_s"] = start_entry.get("flow_kg_per_s", 0.0) - component_result.mass_flow_kg_per_s
            end_entry["flow_kg_per_s"] = end_entry.get("flow_kg_per_s", 0.0) + component_result.mass_flow_kg_per_s

        return boundary_results

    def _prepare_convergence_window(self) -> None:
        if self.convergence_window is None or not self.convergence_window.winfo_exists():
            self.convergence_window = tk.Toplevel(self.root)
            self.convergence_window.title("Convergence Metrics")
            self.convergence_window.transient(self.root)
            self.convergence_window.geometry("820x520")

            frame = ttk.Frame(self.convergence_window, padding=12)
            frame.pack(fill="both", expand=True)

            control_row = ttk.Frame(frame)
            control_row.pack(fill="x", pady=(0, 8))
            ttk.Label(control_row, text="Metric").pack(side="left", padx=(0, 8))
            metric_box = ttk.Combobox(
                control_row,
                textvariable=self.convergence_metric_var,
                state="readonly",
                values=tuple(label for label, _name in self.METRIC_OPTIONS),
                width=32,
            )
            metric_box.pack(side="left")
            metric_box.bind("<<ComboboxSelected>>", lambda _event: self._redraw_convergence_plot())

            ttk.Checkbutton(
                control_row,
                text="Show hydraulic detail",
                variable=self.show_hydraulic_detail_var,
                command=self._redraw_convergence_plot,
            ).pack(side="left", padx=(16, 0))

            self.convergence_canvas = tk.Canvas(
                frame,
                background=self._t["plot_bg"],
                highlightthickness=1,
                highlightbackground=self._t["canvas_hl"],
                width=760,
                height=300,
            )
            self.convergence_canvas.pack(fill="both", expand=True)
            self.convergence_canvas.bind(
                "<Configure>",
                lambda _event: self._redraw_convergence_plot(),
            )
            # Legend is drawn inside the canvas by _draw_history_plot
        else:
            self.convergence_window.deiconify()
            self.convergence_window.lift()

        self.convergence_window.update_idletasks()
        if self.convergence_canvas is not None:
            self.convergence_canvas.after_idle(self._redraw_convergence_plot)

    def _on_solver_progress(self, stage: str, _iteration_index: int, metrics) -> None:
        self.convergence_history[stage].append(metrics)
        if self.convergence_window is not None and self.convergence_window.winfo_exists():
            self.convergence_window.update_idletasks()
        self._redraw_convergence_plot()
        self.root.update()

    def _redraw_convergence_plot(self) -> None:
        if self.convergence_canvas is None:
            return

        show_detail = self.show_hydraulic_detail_var.get()

        # Simple view: outer-iteration convergence (non-isothermal / compressible only).
        # Not used for isothermal (1 outer iteration) — falls through to detail view,
        # which shows the non-zero laminar corrections.
        if not show_detail and len(self.outer_turbulent_final_metrics) > 1:
            metric_name = self.metric_label_to_name[self.convergence_metric_var.get()]

            # Primary series: actual outer convergence criteria (ΔT, Δρ/ρ).
            # Hydraulic final corrections are often 0 (inner solver converges immediately
            # each outer iteration), so they are secondary — shown on the right axis only
            # if non-zero.  Exact zeros in ΔT / Δρ/ρ are clamped to 1e-10 so the log
            # scale shows the convergence drop rather than filtering the point out.
            def _log_safe(vals: list[float]) -> list[float]:
                return [v if v > 0.0 else 1e-10 for v in vals]

            outer_primary: list[tuple[str, list[float], str, int]] = []
            if self.density_history:
                outer_primary.append(("Δρ/ρ (−)", _log_safe(self.density_history), self._t["plot_faint2"], 0))
            if self.temperature_history:
                outer_primary.append(("ΔT (K)", _log_safe(self.temperature_history), self._t["plot_temperature"], 0))

            if outer_primary:
                hydraulic_vals = [getattr(m, metric_name) for m in self.outer_turbulent_final_metrics]
                secondary = None
                if any(v > 0.0 for v in hydraulic_vals):
                    secondary = [("Hydraulic (final)", hydraulic_vals, self._t["plot_turbulent"])]
                self._draw_history_plot(
                    self.convergence_canvas,
                    outer_primary,
                    "outer_residual",
                    secondary_series=secondary,
                    x_label="Outer iteration",
                )
                return

            # No outer criteria available: show hydraulic final corrections directly
            hydraulic_vals = [getattr(m, metric_name) for m in self.outer_turbulent_final_metrics]
            self._draw_history_plot(
                self.convergence_canvas,
                [("Hydraulic (final)", hydraulic_vals, self._t["plot_turbulent"], 0)],
                metric_name,
                secondary_series=None,
                x_label="Outer iteration",
            )
            return

        # Detail view (or isothermal): hydraulic series + optional temperature overlay
        metric_name = self.metric_label_to_name[self.convergence_metric_var.get()]
        history_series: list[tuple[str, list[float], str, int]] = []
        laminar_values = [getattr(metric, metric_name) for metric in self.convergence_history["laminar"]]
        turbulent_values = [getattr(metric, metric_name) for metric in self.convergence_history["turbulent"]]
        if laminar_values:
            history_series.append(("Laminar", laminar_values, self._t["plot_laminar"], 0))
        if turbulent_values:
            history_series.append(
                ("Turbulent", turbulent_values, self._t["plot_turbulent"], len(laminar_values))
            )
        if not history_series:
            self.convergence_canvas.delete("all")
            self.convergence_canvas.create_text(
                20,
                20,
                anchor="nw",
                text="No convergence data yet.",
                fill=self._t["plot_muted"],
            )
            return

        secondary = None
        if self.scene.physics_mode == "compressible":
            sec = []
            if self.temperature_history:
                sec.append(("ΔT (K)", self.temperature_history, self._t["plot_temperature"]))
            if self.density_history:
                sec.append(("Δρ/ρ", self.density_history, self._t["plot_faint2"]))
            secondary = sec or None
        elif self.temperature_history:
            secondary = [("ΔT (K)", self.temperature_history, self._t["plot_temperature"])]
        self._draw_history_plot(
            self.convergence_canvas, history_series, metric_name, secondary_series=secondary
        )


    def _draw_history_plot(
        self,
        canvas: tk.Canvas,
        history_series: list[tuple[str, list[float], str, int]],
        metric_name: str,
        secondary_series: list[tuple[str, list[float], str]] | None = None,
        x_label: str = "Iteration",
    ) -> None:
        canvas.delete("all")
        width = int(canvas.winfo_width() or canvas["width"])
        height = int(canvas.winfo_height() or canvas["height"])
        if width < 160 or height < 120:
            canvas.create_text(
                20,
                20,
                anchor="nw",
                text="Waiting for plot area...",
                fill=self._t["plot_muted"],
            )
            return

        has_secondary = bool(
            secondary_series and any(values for _, values, _ in secondary_series)
        )
        left = 70
        right = width - (70 if has_secondary else 20)
        top = 20
        bottom = height - 45
        if right <= left or bottom <= top:
            return

        all_values = [
            value
            for _label, values, _color, _offset in history_series
            for value in values
            if value > 0.0
        ]
        if not all_values:
            total_iters = sum(len(values) for _label, values, _color, _offset in history_series)
            msg = (
                "Converged immediately — correction below machine precision."
                if total_iters > 0
                else "No convergence data yet."
            )
            canvas.delete("all")
            canvas.create_text(
                20, 20, anchor="nw", text=msg, fill=self._t["plot_muted"],
            )
            return

        min_log = math.log10(min(all_values))
        max_log = math.log10(max(all_values))
        if math.isclose(min_log, max_log):
            min_log -= 1.0
            max_log += 1.0

        max_index = max(
            offset + len(values) - 1
            for _label, values, _color, offset in history_series
        )
        x_den = max(max_index, 1)
        total_iterations = max_index + 1

        canvas.create_line(left, top, left, bottom, fill=self._t["plot_axis"], width=1.5)
        canvas.create_line(left, bottom, right, bottom, fill=self._t["plot_axis"], width=1.5)

        decade_min = math.floor(min_log)
        decade_max = math.ceil(max_log)
        num_decades = decade_max - decade_min
        decade_step = 1 if num_decades <= 6 else 2 if num_decades <= 12 else 3
        decade_ticks = range(decade_min, decade_max + 1, decade_step)
        for decade in decade_ticks:
            if decade < min_log - 0.01 or decade > max_log + 0.01:
                continue
            frac = (max_log - decade) / (max_log - min_log)
            y = top + (bottom - top) * frac
            value = 10 ** decade
            canvas.create_line(left, y, right, y, fill=self._t["plot_grid"])
            canvas.create_text(left - 8, y, text=f"{value:.0e}", anchor="e", font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        tick_count = min(10, total_iterations)
        tick_indices = sorted(
            {
                round(i * (total_iterations - 1) / max(tick_count - 1, 1))
                for i in range(tick_count)
            }
        )
        for i in tick_indices:
            x = left + (right - left) * (i / x_den)
            canvas.create_line(x, bottom, x, bottom + 4, fill=self._t["plot_axis"])
            canvas.create_text(x, bottom + 16, text=str(i + 1), anchor="n", font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        if history_series and len(history_series) == 2:
            transition_index = history_series[1][3]
            if 0 < transition_index <= total_iterations - 1:
                transition_x = left + (right - left) * (transition_index / x_den)
                canvas.create_line(
                    transition_x,
                    top,
                    transition_x,
                    bottom,
                    fill=self._t["plot_faint"],
                    dash=(4, 3),
                )

        canvas.create_text((left + right) / 2, height - 10, text=x_label, anchor="s", fill=self._t["plot_text"])
        canvas.create_text(
            16,
            (top + bottom) / 2,
            text=self._pretty_metric_name(metric_name),
            angle=90,
            fill=self._t["plot_text"],
        )

        for _label, values, color, offset in history_series:
            points: list[float] = []
            for idx, value in enumerate(values):
                safe_value = value if value > 0.0 else min(all_values)
                value_log = math.log10(safe_value)
                x = left + (right - left) * ((offset + idx) / x_den)
                y = top + (max_log - value_log) * (bottom - top) / (max_log - min_log)
                points.extend((x, y))
            if len(points) >= 4:
                canvas.create_line(*points, fill=color, width=2, smooth=False)
                if len(values) <= 20:
                    for idx, value in enumerate(values):
                        safe_value = value if value > 0.0 else min(all_values)
                        value_log = math.log10(safe_value)
                        x = left + (right - left) * ((offset + idx) / x_den)
                        y = top + (max_log - value_log) * (bottom - top) / (max_log - min_log)
                        canvas.create_oval(
                            x - 2,
                            y - 2,
                            x + 2,
                            y + 2,
                            fill=color,
                            outline=color,
                        )
            elif len(points) == 2:
                x, y = points
                canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill=color, outline=color)

        if has_secondary and secondary_series:
            sec_all_values = [
                v for _, values, _ in secondary_series for v in values if v > 0.0
            ]
            if sec_all_values:
                sec_min_log = math.log10(min(sec_all_values))
                sec_max_log = math.log10(max(sec_all_values))
                if math.isclose(sec_min_log, sec_max_log):
                    sec_min_log -= 1.0
                    sec_max_log += 1.0
                sec_decade_min = math.floor(sec_min_log)
                sec_decade_max = math.ceil(sec_max_log)
                sec_num_decades = sec_decade_max - sec_decade_min
                sec_decade_step = 1 if sec_num_decades <= 6 else 2 if sec_num_decades <= 12 else 3

                temp_color = secondary_series[0][2]
                canvas.create_line(right, top, right, bottom, fill=temp_color, width=1.5)

                for decade in range(sec_decade_min, sec_decade_max + 1, sec_decade_step):
                    if decade < sec_min_log - 0.01 or decade > sec_max_log + 0.01:
                        continue
                    frac = (sec_max_log - decade) / (sec_max_log - sec_min_log)
                    y = top + (bottom - top) * frac
                    value = 10 ** decade
                    canvas.create_line(right, y, right + 4, y, fill=temp_color, width=1)
                    canvas.create_text(
                        right + 6, y, text=f"{value:.0e}", anchor="w",
                        font=("TkDefaultFont", 8), fill=temp_color,
                    )

                canvas.create_text(
                    width - 12,
                    (top + bottom) / 2,
                    text="max |ΔT| (K)",
                    angle=90,
                    fill=temp_color,
                    font=("TkDefaultFont", 9),
                )

                for _label, values, color in secondary_series:
                    if not values:
                        continue
                    n = len(values)  # used for x-spacing even if some values are zero
                    # Collect only positive points (skip exact-zero "converged" entries)
                    pos_pts = [
                        (i, v) for i, v in enumerate(values) if v > 0.0
                    ]
                    if not pos_pts:
                        continue
                    pts: list[float] = []
                    for i, value in pos_pts:
                        vlog = math.log10(value)
                        x = right if n == 1 else left + (right - left) * (i / (n - 1))
                        y = top + (sec_max_log - vlog) * (bottom - top) / (sec_max_log - sec_min_log)
                        pts.extend((x, y))
                    if len(pts) >= 4:
                        canvas.create_line(*pts, fill=color, width=2, smooth=False, dash=(6, 3))
                    for i, value in pos_pts:
                        vlog = math.log10(value)
                        x = right if n == 1 else left + (right - left) * (i / (n - 1))
                        y = top + (sec_max_log - vlog) * (bottom - top) / (sec_max_log - sec_min_log)
                        size = 4
                        canvas.create_polygon(
                            x, y - size, x + size, y, x, y + size, x - size, y,
                            fill=color, outline=color,
                        )

        # Canvas legend (bottom-left, inside plot area)
        legend_items: list[tuple[str, str, str]] = []  # (label, color, style)
        for label, _values, color, _offset in history_series:
            legend_items.append((label, color, "line"))
        if has_secondary and secondary_series:
            for label, _values, color in secondary_series:
                legend_items.append((label + " (right axis)", color, "dashed"))
        if legend_items:
            lx = left + 8
            ly = top + 8
            for leg_label, leg_color, leg_style in legend_items:
                if leg_style == "dashed":
                    canvas.create_line(lx, ly + 5, lx + 20, ly + 5, fill=leg_color, width=2, dash=(5, 3))
                else:
                    canvas.create_line(lx, ly + 5, lx + 20, ly + 5, fill=leg_color, width=2)
                canvas.create_text(
                    lx + 24, ly + 5, anchor="w", text=leg_label,
                    fill=self._t["plot_text"], font=("TkDefaultFont", 8),
                )
                ly += 16

    @staticmethod
    def _pretty_metric_name(metric_name: str) -> str:
        if metric_name == "temperature_delta_k":
            return "max |ΔT| (K)"
        if metric_name == "density_rel_change":
            return "max |Δρ/ρ|"
        if metric_name == "outer_residual":
            return "Outer residual"
        labels = {name: label for label, name in NetSimGui.METRIC_OPTIONS}
        return labels.get(metric_name, metric_name)

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    app = NetSimGui()
    app.run()


if __name__ == "__main__":
    main()
