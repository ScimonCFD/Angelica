from __future__ import annotations

import math
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import sv_ttk

from angelica.core.components import FITTING_PRESET_LIBRARY
from angelica.io import export_solve_result_workbook

from .io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file
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
            "plot_laminar":   "#1d3557",
            "plot_turbulent": "#c1121f",
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
            "plot_laminar":   "#3a9fd4",
            "plot_turbulent": "#e8633a",
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
        self.convergence_window: tk.Toplevel | None = None
        self.convergence_canvas: tk.Canvas | None = None
        self._dark = False
        self._unit_system_key = "si"
        self.root = tk.Tk()
        sv_ttk.set_theme("light")
        self.root.title("Angelica GUI")
        self.root.geometry("1100x700")
        self.root.minsize(900, 600)
        self._set_window_icon()

        self.metric_label_to_name = {label: name for label, name in self.METRIC_OPTIONS}
        self.metric_name_to_label = {name: label for label, name in self.METRIC_OPTIONS}
        self.convergence_metric_var = tk.StringVar(
            master=self.root,
            value=self.metric_name_to_label["pressure_correction_abs_pa"],
        )
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
        file_menu.add_command(label="Open", command=self._open_scene)
        file_menu.add_command(label="Export Results Report", command=self._export_results_report)
        file_menu.add_separator()
        file_menu.add_command(label="Close", command=self.root.destroy)
        menu_bar.add_cascade(label="File", menu=file_menu)

        material_menu = tk.Menu(menu_bar, tearoff=False)
        material_menu.add_command(label="Define Material", command=self._open_material_dialog)
        menu_bar.add_cascade(label="Material", menu=material_menu)

        physics_menu = tk.Menu(menu_bar, tearoff=False)
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
        self.tool_var.set("No tool selected")
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
            self.scene = load_scene_from_file(file_path)
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Open failed", f"Could not open case:\n{exc}")
            return

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
        self.tool_var.set("No tool selected")
        self._refresh_global_summaries()
        self._redraw_scene()
        self.status_var.set(f"Opened GUI case: {file_path}")

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
        density_var = tk.StringVar(
            master=dialog,
            value=material.get("density_kg_per_m3", ""),
        )
        viscosity_var = tk.StringVar(
            master=dialog,
            value=material.get("viscosity_pa_s", ""),
        )

        ttk.Label(frame, text="Definition").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        mode_row = ttk.Frame(frame)
        mode_row.grid(row=0, column=1, sticky="w", pady=4)
        ttk.Radiobutton(mode_row, text="Library", variable=mode_var, value="library").pack(side="left")
        ttk.Radiobutton(mode_row, text="Custom", variable=mode_var, value="custom").pack(side="left", padx=(8, 0))

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

        ttk.Label(frame, text="Density (kg/m^3)").grid(row=3, column=0, sticky="w", padx=(0, 8), pady=4)
        density_entry = ttk.Entry(frame, textvariable=density_var, width=26)
        density_entry.grid(row=3, column=1, sticky="ew", pady=4)

        ttk.Label(frame, text="Viscosity (Pa·s)").grid(row=4, column=0, sticky="w", padx=(0, 8), pady=4)
        viscosity_entry = ttk.Entry(frame, textvariable=viscosity_var, width=26)
        viscosity_entry.grid(row=4, column=1, sticky="ew", pady=4)

        frame.columnconfigure(1, weight=1)

        def apply_library_selection(_event: tk.Event | None = None) -> None:
            preset = self.MATERIAL_LIBRARY[library_var.get()]
            name_var.set(preset["name"])
            density_var.set(preset["density_kg_per_m3"])
            viscosity_var.set(preset["viscosity_pa_s"])

        def sync_mode_state(*_args: object) -> None:
            is_library = mode_var.get() == "library"
            library_box.configure(state="readonly" if is_library else "disabled")
            editable_state = "disabled" if is_library else "normal"
            name_entry.configure(state=editable_state)
            density_entry.configure(state=editable_state)
            viscosity_entry.configure(state=editable_state)
            if is_library and library_var.get():
                apply_library_selection()

        library_box.bind("<<ComboboxSelected>>", apply_library_selection)
        mode_var.trace_add("write", sync_mode_state)
        sync_mode_state()

        button_row = ttk.Frame(frame)
        button_row.grid(row=5, column=0, columnspan=2, sticky="e", pady=(10, 0))
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

        current_laminar_iterations = str(
            self.scene.solver_settings.get("laminar_iterations", "")
        )
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
        button_row.grid(row=17, column=0, columnspan=2, sticky="e", pady=(10, 0))
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
    ) -> None:
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
            "definition_mode": mode_var.get().strip(),
            "library_key": library_var.get().strip() if mode_var.get() == "library" else "",
            "name": name_var.get().strip(),
            "density_kg_per_m3": density_var.get().strip(),
            "viscosity_pa_s": viscosity_var.get().strip(),
        }
        if mode_var.get() == "library" and not material["library_key"]:
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
        self.scene.update_material(material)
        self._refresh_global_summaries()
        self.status_var.set(f"Material set to {material['name']}.")
        dialog.destroy()

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
    ) -> None:
        laminar_iterations_text = laminar_iterations_var.get().strip()
        if laminar_iterations_text:
            try:
                laminar_iterations: int | None = int(laminar_iterations_text)
            except ValueError:
                messagebox.showerror(
                    "Invalid numerics",
                    "Laminar iterations must be an integer.",
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
        else:
            laminar_iterations = None

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
        density = self.scene.material.get("density_kg_per_m3", "").strip()
        viscosity = self.scene.material.get("viscosity_pa_s", "").strip()
        lines = [name]
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
            f"f-tol={self._fmt_sci(colebrook_tol)} −\n"
            f"V*-tol={self._fmt_sci(velocity_loop_tol)} m/s\n"
            f"ΔP-tol={self._fmt_sci(dp_tol)} Pa\n"
            f"cont-tol={self._fmt_sci(continuity_tol)} −"
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
        self._prepare_convergence_window()
        result = solver.solve(case, progress_callback=self._on_solver_progress)
        self.latest_result = result
        self.latest_boundary_results = self._build_boundary_results(case, result)
        self.convergence_history = {
            "laminar": list(result.laminar_metrics),
            "turbulent": list(result.turbulent_metrics),
        }
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
        if node_id is None:
            return

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
        self.scene.remove_node(node_id)
        self.selected_node_id = None
        self.latest_boundary_results = {}
        self.latest_result = None
        self._redraw_scene()
        self.status_var.set(f"Deleted {node_type} #{node_id} and its connections.")

    def _on_canvas_release(self, event: tk.Event) -> None:
        if self.moving_node_id is not None:
            moved_node = self.scene.get_node(self.moving_node_id)
            self.moving_node_id = None
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

        try:
            link = self.scene.add_link(source_node_id, target_node_id)
        except ValueError as exc:
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

        self.canvas.create_line(
            start_x,
            start_y,
            end_x,
            end_y,
            fill=self._t["link"],
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

        components_list.bind(
            "<<ListboxSelect>>",
            lambda _event: self._render_link_component_properties(
                link.link_id,
                components_list,
                properties_frame,
            ),
        )

        button_row = ttk.Frame(container)
        button_row.grid(row=2, column=0, columnspan=3, sticky="e", pady=(12, 0))
        ttk.Button(button_row, text="Close", command=dialog.destroy).pack(side="right", padx=(6, 0))
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
        win.geometry("760x440")

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

        win.transient(self.root)
        win.grab_set()
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
                lbl = "pump" if comp.component_type == "pump" else "fitting"
                clr = pump_color if comp.component_type == "pump" else fitting_color
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
            return "\n".join(lines)

        result_data = self.latest_boundary_results.get(node.node_id)
        if result_data and "pressure_pa" in result_data:
            _, p_factor = self._unit_quantities()["pressure"]
            p_unit = self._unit_label("pressure")
            return f"P={self._fmt(result_data['pressure_pa'] * p_factor)} {p_unit}"
        return ""

    def _build_network_case_from_scene(self) -> NetworkCase:
        return build_network_case_from_scene(self.scene)

    def _build_boundary_results(self, case: NetworkCase, result) -> dict[int, dict[str, float]]:
        boundary_results: dict[int, dict[str, float]] = {}
        for node_id, pressure in result.node_pressures_pa.items():
            boundary_results[node_id] = {"pressure_pa": pressure}

        for component, component_result in zip(case.components, result.component_flows):
            start_entry = boundary_results.setdefault(component.start_node, {"pressure_pa": result.node_pressures_pa.get(component.start_node, 0.0)})
            end_entry = boundary_results.setdefault(component.end_node, {"pressure_pa": result.node_pressures_pa.get(component.end_node, 0.0)})
            start_entry["flow_kg_per_s"] = start_entry.get("flow_kg_per_s", 0.0) - component_result.mass_flow_kg_per_s
            end_entry["flow_kg_per_s"] = end_entry.get("flow_kg_per_s", 0.0) + component_result.mass_flow_kg_per_s

        return boundary_results

    def _prepare_convergence_window(self) -> None:
        self.convergence_history = {"laminar": [], "turbulent": []}

        if self.convergence_window is None or not self.convergence_window.winfo_exists():
            self.convergence_window = tk.Toplevel(self.root)
            self.convergence_window.title("Convergence Metrics")
            self.convergence_window.transient(self.root)
            self.convergence_window.geometry("820x460")

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

            self.convergence_canvas = tk.Canvas(
                frame,
                background=self._t["plot_bg"],
                highlightthickness=1,
                highlightbackground=self._t["canvas_hl"],
                width=760,
                height=340,
            )
            self.convergence_canvas.pack(fill="both", expand=True)
            self.convergence_canvas.bind(
                "<Configure>",
                lambda _event: self._redraw_convergence_plot(),
            )

            legend = ttk.Frame(frame)
            legend.pack(fill="x", pady=(8, 0))
            for label, color in (("Laminar", self._t["plot_laminar"]), ("Turbulent", self._t["plot_turbulent"])):
                swatch = tk.Canvas(legend, width=16, height=10, highlightthickness=0)
                swatch.create_line(1, 5, 15, 5, fill=color, width=3)
                swatch.pack(side="left", padx=(0, 4))
                ttk.Label(legend, text=label).pack(side="left", padx=(0, 12))
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

        self._draw_history_plot(self.convergence_canvas, history_series, metric_name)

    def _draw_history_plot(
        self,
        canvas: tk.Canvas,
        history_series: list[tuple[str, list[float], str, int]],
        metric_name: str,
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

        left = 70
        right = width - 20
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

        for i in range(5):
            frac = i / 4 if 4 else 0
            y = top + (bottom - top) * frac
            value_log = max_log - (max_log - min_log) * frac
            value = 10**value_log
            canvas.create_line(left, y, right, y, fill=self._t["plot_grid"])
            canvas.create_text(left - 8, y, text=f"{value:.1e}", anchor="e", font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        tick_count = min(10, total_iterations)
        tick_indices = sorted(
            {
                round(i * max_index / max(tick_count - 1, 1))
                for i in range(tick_count)
            }
        )
        for i in tick_indices:
            x = left + (right - left) * (i / x_den)
            canvas.create_line(x, bottom, x, bottom + 4, fill=self._t["plot_axis"])
            canvas.create_text(x, bottom + 16, text=str(i + 1), anchor="n", font=("TkDefaultFont", 8), fill=self._t["plot_text"])

        if history_series and len(history_series) == 2:
            transition_index = history_series[1][3]
            if 0 < transition_index <= max_index:
                transition_x = left + (right - left) * (transition_index / x_den)
                canvas.create_line(
                    transition_x,
                    top,
                    transition_x,
                    bottom,
                    fill=self._t["plot_faint"],
                    dash=(4, 3),
                )
                canvas.create_text(
                    transition_x + 4,
                    top + 10,
                    text="turbulent",
                    anchor="nw",
                    fill=self._t["plot_faint2"],
                    font=("TkDefaultFont", 8),
                )

        canvas.create_text((left + right) / 2, height - 10, text="Iteration", anchor="s", fill=self._t["plot_text"])
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

    @staticmethod
    def _pretty_metric_name(metric_name: str) -> str:
        labels = {name: label for label, name in NetSimGui.METRIC_OPTIONS}
        return labels.get(metric_name, metric_name)

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    app = NetSimGui()
    app.run()


if __name__ == "__main__":
    main()
