from __future__ import annotations

import math
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from netSim.core.components import FITTING_PRESET_LIBRARY
from netSim.io import export_solve_result_workbook

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
        ("Max abs pressure correction (Pa)", "pressure_correction_abs_pa"),
        ("Mean abs pressure correction (Pa)", "pressure_correction_mean_abs_pa"),
        ("Max rel pressure correction", "pressure_correction_rel"),
        ("Max nodal mass imbalance (kg/s)", "max_nodal_mass_imbalance_kg_per_s"),
        ("Max abs flow (kg/s)", "mass_flow_max_abs_kg_per_s"),
    )

    def __init__(self) -> None:
        self.scene = CanvasScene()
        self.drag_source_node_id: int | None = None
        self.drag_line_id: int | None = None
        self.moving_node_id: int | None = None
        self.middle_pan_anchor: tuple[float, float] | None = None
        self.view_offset_x = 0.0
        self.view_offset_y = 0.0
        self.view_scale = 1.0
        self.latest_result = None
        self.latest_boundary_results: dict[int, dict[str, float]] = {}
        self.convergence_window: tk.Toplevel | None = None
        self.convergence_canvas: tk.Canvas | None = None
        self.root = tk.Tk()
        self.root.title("netSim GUI")
        self.root.geometry("1100x700")
        self.root.minsize(900, 600)

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
            background="#fbfaf4",
            highlightthickness=1,
            highlightbackground="#b8b2a7",
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
            title="Open netSim GUI case",
            filetypes=(
                ("netSim GUI case", "*.gui.json"),
                ("JSON files", "*.json"),
                ("All files", "*.*"),
            ),
        )
        if not file_path:
            return

        try:
            self.scene = load_scene_from_file(file_path)
        except Exception as exc:  # pragma: no cover - UI feedback path
            messagebox.showerror("Open failed", f"Could not open case:\n{exc}")
            return

        self.canvas.delete("all")
        self.drag_source_node_id = None
        self.drag_line_id = None
        self.moving_node_id = None
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
        button_row.grid(row=11, column=0, columnspan=2, sticky="e", pady=(10, 0))
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
            }
        )
        self._refresh_global_summaries()
        self.status_var.set(
            "Numerics updated: "
            f"laminar={'auto' if laminar_iterations is None else laminar_iterations}, "
            f"turbulent={turbulent_iterations}, "
            f"Explicit relaxation (alpha={alpha:g}), "
            f"Colebrook={self.COLEBROOK_FRICTION_STRATEGY_LIBRARY[colebrook_strategy]['name']}, "
            f"{self.FRICTION_FACTOR_METHOD_LIBRARY[friction_method]['name']} "
            f"(max={friction_max_iterations}) for friction, "
            f"{self.VELOCITY_LOOP_METHOD_LIBRARY[velocity_method]['name']} "
            f"(max={velocity_max_iterations})."
        )
        dialog.destroy()

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
        return (
            f"laminar={'auto' if laminar_iterations is None else laminar_iterations}\n"
            f"turbulent={turbulent_iterations}\n"
            f"Explicit\nalpha={alpha}\n"
            f"colebrook={colebrook_strategy_name}\n"
            f"friction={friction_method_name} ({friction_max_iterations})\n"
            f"velocity={velocity_method_name} ({velocity_max_iterations})"
        )

    def _run_simulation(self) -> None:
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
            return

        node = self.scene.get_node(node_id)
        if node is None:
            return

        self.moving_node_id = node_id
        self.status_var.set(
            f"Moving {node.node_type} #{node.node_id}. Drag with left mouse button."
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
            fill="#6c757d",
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
        radius = 24
        canvas_x, canvas_y = self._scene_to_canvas(node.x, node.y)
        x0 = canvas_x - radius
        y0 = canvas_y - radius
        x1 = canvas_x + radius
        y1 = canvas_y + radius

        fill_color = self._node_fill(node.node_type)
        label = self._node_label(node)

        node_tag = f"node_{node.node_id}"

        self.canvas.create_oval(
            x0,
            y0,
            x1,
            y1,
            fill=fill_color,
            outline="#2e2a24",
            width=2,
            tags=(node_tag, "node"),
        )
        self.canvas.create_text(
            canvas_x,
            canvas_y - 2,
            text=label,
            font=("TkDefaultFont", 9, "bold"),
            tags=(node_tag, "node"),
        )
        self.canvas.create_text(
            canvas_x,
            canvas_y + 12,
            text=str(node.node_id),
            font=("TkDefaultFont", 8),
            tags=(node_tag, "node"),
        )

        summary = self._node_summary_text(node)
        if summary:
            self.canvas.create_text(
                canvas_x,
                canvas_y + 40,
                text=summary,
                font=("TkDefaultFont", 8),
                fill="#3d3a35",
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
            fill="#4f5d75",
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

            ttk.Label(container, text="Pressure").grid(row=1, column=0, sticky="w", pady=4)
            pressure_var = tk.StringVar(value=node.properties.get("pressure", ""))
            pressure_entry = ttk.Entry(container, textvariable=pressure_var, width=20)
            pressure_entry.grid(row=1, column=1, sticky="ew", pady=4)
            entries["pressure"] = pressure_var

            ttk.Label(container, text="Flow").grid(row=2, column=0, sticky="w", pady=4)
            flow_var = tk.StringVar(value=node.properties.get("flow", ""))
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
        ttk.Button(button_row, text="Close", command=dialog.destroy).pack(side="right")

        dialog.update_idletasks()
        dialog.grab_set()
        dialog.focus_set()

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
            var = tk.StringVar(value=value)
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
        diameter_var = tk.StringVar(value=component.properties.get("diameter_m", ""))
        mode_var = tk.StringVar(value=component.properties.get("fitting_mode", "manual"))
        preset_var = tk.StringVar(
            value=component.properties.get("fitting_preset", "regular_90_flanged")
        )
        loss_var = tk.StringVar(value=component.properties.get("loss_coefficient", "1.5"))

        entries["diameter_m"] = diameter_var
        entries["fitting_mode"] = mode_var
        entries["fitting_preset"] = preset_var
        entries["loss_coefficient"] = loss_var

        ttk.Label(properties_frame, text="Diameter (m)").grid(
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
        diameter_var = tk.StringVar(value=component.properties.get("diameter_m", ""))

        ttk.Label(properties_frame, text="Diameter (m)").grid(
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
            foreground="#555555",
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
        properties = {key: value.get().strip() for key, value in entries.items()}
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
            "diameter_m": diameter_var.get().strip(),
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
        properties = {key: value.get().strip() for key, value in entries.items()}
        updated_node = self.scene.update_node_properties(node_id, properties)
        self._redraw_scene()
        self.status_var.set(
            f"Updated properties for {updated_node.node_type} #{updated_node.node_id}."
        )
        dialog.destroy()

    @staticmethod
    def _node_fill(node_type: str) -> str:
        colors = {
            "source": "#8ecae6",
            "sink": "#f28482",
            "junction": "#d9d9d9",
        }
        return colors[node_type]

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

    @staticmethod
    def _pretty_field_name(field_name: str) -> str:
        return field_name.replace("_", " ").replace(" m", " (m)").capitalize()

    def _node_summary_text(self, node: CanvasNode) -> str:
        lines: list[str] = []
        if node.node_type in {"source", "sink"}:
            condition_type = node.properties.get("condition_type", "pressure")
            if condition_type == "pressure":
                pressure = node.properties.get("pressure", "").strip()
                if pressure:
                    lines.append(f"P={pressure}")
            else:
                flow = node.properties.get("flow", "").strip()
                if flow:
                    lines.append(f"Q={flow}")

            result_data = self.latest_boundary_results.get(node.node_id)
            if result_data is not None:
                if condition_type == "pressure" and "flow_kg_per_s" in result_data:
                    lines.append(f"Q={result_data['flow_kg_per_s']:.3f} kg/s")
                elif condition_type == "flow" and "pressure_pa" in result_data:
                    lines.append(f"P={result_data['pressure_pa']:.1f} Pa")
            return "\n".join(lines)

        label = node.properties.get("label", "").strip()
        return label

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
                background="white",
                highlightthickness=1,
                highlightbackground="#b8b2a7",
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
            for label, color in (("Laminar", "#1d3557"), ("Turbulent", "#c1121f")):
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
            history_series.append(("Laminar", laminar_values, "#1d3557", 0))
        if turbulent_values:
            history_series.append(
                ("Turbulent", turbulent_values, "#c1121f", len(laminar_values))
            )
        if not history_series:
            self.convergence_canvas.delete("all")
            self.convergence_canvas.create_text(
                20,
                20,
                anchor="nw",
                text="No convergence data yet.",
                fill="#555555",
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
                fill="#555555",
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

        canvas.create_line(left, top, left, bottom, fill="#333333", width=1.5)
        canvas.create_line(left, bottom, right, bottom, fill="#333333", width=1.5)

        for i in range(5):
            frac = i / 4 if 4 else 0
            y = top + (bottom - top) * frac
            value_log = max_log - (max_log - min_log) * frac
            value = 10**value_log
            canvas.create_line(left, y, right, y, fill="#e6e6e6")
            canvas.create_text(left - 8, y, text=f"{value:.1e}", anchor="e", font=("TkDefaultFont", 8))

        tick_count = min(10, total_iterations)
        tick_indices = sorted(
            {
                round(i * max_index / max(tick_count - 1, 1))
                for i in range(tick_count)
            }
        )
        for i in tick_indices:
            x = left + (right - left) * (i / x_den)
            canvas.create_line(x, bottom, x, bottom + 4, fill="#333333")
            canvas.create_text(x, bottom + 16, text=str(i + 1), anchor="n", font=("TkDefaultFont", 8))

        if history_series and len(history_series) == 2:
            transition_index = history_series[1][3]
            if 0 < transition_index <= max_index:
                transition_x = left + (right - left) * (transition_index / x_den)
                canvas.create_line(
                    transition_x,
                    top,
                    transition_x,
                    bottom,
                    fill="#999999",
                    dash=(4, 3),
                )
                canvas.create_text(
                    transition_x + 4,
                    top + 10,
                    text="turbulent",
                    anchor="nw",
                    fill="#666666",
                    font=("TkDefaultFont", 8),
                )

        canvas.create_text((left + right) / 2, height - 10, text="Iteration", anchor="s")
        canvas.create_text(
            16,
            (top + bottom) / 2,
            text=self._pretty_metric_name(metric_name),
            angle=90,
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
