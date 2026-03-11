import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


REQUIRED_COLUMNS = {
    "time_days",
    "layer",
    "i",
    "j",
    "pressure_psia",
    "sw",
    "sg",
}

MAP_FIELDS = {
    "pressure": "pressure_psia",
    "sw": "sw",
    "sg": "sg",
}


class MapViewerApp(tk.Tk):
    def __init__(self, initial_file: str | None = None) -> None:
        super().__init__()
        self.title("BOAST II Map Viewer")
        self.geometry("1100x760")

        self.df: pd.DataFrame | None = None
        self.i_values: list[int] = []
        self.j_values: list[int] = []
        self.layer_values: list[int] = []
        self.time_values: list[float] = []
        self.slice_values: list[int] = []
        self.field_ranges: dict[str, tuple[float, float] | None] = {}
        self.current_grid: pd.DataFrame | None = None
        self.current_field_key: str = "pressure"
        self.time_index: int = 0
        self.animating: bool = False
        self.animation_job: str | None = None
        self.closing: bool = False

        self.file_var = tk.StringVar(value="No file loaded")
        self.map_var = tk.StringVar(value="pressure")
        self.direction_var = tk.StringVar(value="XY")
        self.time_var = tk.StringVar(value="")
        self.slice_var = tk.IntVar(value=1)
        self.slice_label_var = tk.StringVar(value="Layer")
        self.speed_var = tk.StringVar(value="0.5")
        self.hover_var = tk.StringVar(value="Hover over a block to see value")

        self._build_ui()
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.after(100, self._ensure_window_visible)
        if initial_file:
            self.after(150, lambda: self.load_file(initial_file))

    def _ensure_window_visible(self) -> None:
        self.deiconify()
        self.lift()
        self.focus_force()
        self.attributes("-topmost", True)
        self.after(150, lambda: self.attributes("-topmost", False))

    def _build_ui(self) -> None:
        controls = ttk.Frame(self, padding=10)
        controls.pack(side=tk.TOP, fill=tk.X)

        open_button = ttk.Button(controls, text="Open _maps.csv", command=self.open_file)
        open_button.grid(row=0, column=0, padx=(0, 10), pady=4, sticky="w")

        ttk.Label(controls, textvariable=self.file_var).grid(
            row=0, column=1, columnspan=8, pady=4, sticky="w"
        )

        ttk.Label(controls, text="Map").grid(row=1, column=0, pady=4, sticky="w")
        map_combo = ttk.Combobox(
            controls,
            textvariable=self.map_var,
            values=["pressure", "sw", "sg"],
            state="readonly",
            width=12,
        )
        map_combo.grid(row=1, column=1, padx=(0, 12), pady=4, sticky="w")
        map_combo.bind("<<ComboboxSelected>>", lambda _: self.update_plot())

        ttk.Label(controls, text="Direction").grid(row=1, column=2, pady=4, sticky="w")
        direction_combo = ttk.Combobox(
            controls,
            textvariable=self.direction_var,
            values=["XY", "XZ", "YZ"],
            state="readonly",
            width=8,
        )
        direction_combo.grid(row=1, column=3, padx=(0, 12), pady=4, sticky="w")
        direction_combo.bind("<<ComboboxSelected>>", lambda _: self._on_direction_change())

        ttk.Label(controls, text="Time (days)").grid(row=1, column=4, pady=4, sticky="w")
        self.time_combo = ttk.Combobox(
            controls,
            textvariable=self.time_var,
            values=[],
            state="readonly",
            width=12,
        )
        self.time_combo.grid(row=1, column=5, padx=(0, 12), pady=4, sticky="w")
        self.time_combo.bind("<<ComboboxSelected>>", lambda _: self.update_plot())

        self.slice_label = ttk.Label(controls, textvariable=self.slice_label_var)
        self.slice_label.grid(row=1, column=6, pady=4, sticky="w")

        self.slice_scale = tk.Scale(
            controls,
            from_=1,
            to=1,
            orient=tk.HORIZONTAL,
            variable=self.slice_var,
            command=lambda _: self.update_plot(),
            showvalue=True,
            length=220,
        )
        self.slice_scale.grid(row=1, column=7, padx=(0, 8), pady=4, sticky="w")

        self.slice_info = ttk.Label(controls, text="")
        self.slice_info.grid(row=1, column=8, pady=4, sticky="w")

        self.animate_button = ttk.Button(controls, text="Animate", command=self.toggle_animation)
        self.animate_button.grid(row=1, column=9, padx=(8, 0), pady=4, sticky="w")

        ttk.Label(controls, text="Speed (s)").grid(row=1, column=10, padx=(10, 4), pady=4, sticky="w")
        self.speed_combo = ttk.Combobox(
            controls,
            textvariable=self.speed_var,
            values=["0.1","0.2", "0.5", "1.0", "2.0"],
            state="readonly",
            width=6,
        )
        self.speed_combo.grid(row=1, column=11, pady=4, sticky="w")

        ttk.Label(controls, textvariable=self.hover_var).grid(
            row=2, column=0, columnspan=12, pady=(6, 0), sticky="w"
        )

        self.figure, self.ax = plt.subplots(figsize=(9, 6))
        self.colorbar = None
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.mpl_connect("motion_notify_event", self._on_mouse_move)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.figure.tight_layout()

    def _on_mouse_move(self, event) -> None:
        if event.inaxes != self.ax or self.current_grid is None or event.xdata is None or event.ydata is None:
            self.hover_var.set("Hover over a block to see value")
            return

        col_idx = int(round(event.xdata))
        row_idx = int(round(event.ydata))

        if (
            row_idx < 0
            or row_idx >= len(self.current_grid.index)
            or col_idx < 0
            or col_idx >= len(self.current_grid.columns)
        ):
            self.hover_var.set("Hover over a block to see value")
            return

        direction = self.direction_var.get().upper()
        slice_value = self.slice_var.get()
        x_value = self.current_grid.columns[col_idx]
        y_value = self.current_grid.index[row_idx]
        value = self.current_grid.iat[row_idx, col_idx]
        value_text = "nan" if pd.isna(value) else f"{float(value):.6g}"

        if direction == "XY":
            info = (
                f"i={x_value}, j={y_value}, layer={slice_value} | "
                f"{self.current_field_key}={value_text}"
            )
        elif direction == "XZ":
            info = (
                f"i={x_value}, j={slice_value}, layer={y_value} | "
                f"{self.current_field_key}={value_text}"
            )
        else:
            info = (
                f"i={slice_value}, j={x_value}, layer={y_value} | "
                f"{self.current_field_key}={value_text}"
            )

        self.hover_var.set(info)

    def open_file(self) -> None:
        if self.animating:
            self.stop_animation()

        file_path = filedialog.askopenfilename(
            title="Select a _maps.csv file",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialdir=os.getcwd(),
        )
        if not file_path:
            return

        self.load_file(file_path)

    def load_file(self, file_path: str) -> None:
        if self.animating:
            self.stop_animation()

        if not os.path.isfile(file_path):
            messagebox.showerror("File not found", f"Could not find file:\n{file_path}")
            return

        try:
            df = pd.read_csv(file_path)
        except Exception as exc:
            messagebox.showerror("Read error", f"Could not read file:\n{exc}")
            return

        missing = REQUIRED_COLUMNS.difference(df.columns)
        if missing:
            messagebox.showerror("Missing columns", f"Missing columns: {sorted(missing)}")
            return

        try:
            df = df.copy()
            df["time_days"] = pd.to_numeric(df["time_days"], errors="coerce")
            df["layer"] = pd.to_numeric(df["layer"], errors="coerce").astype("Int64")
            df["i"] = pd.to_numeric(df["i"], errors="coerce").astype("Int64")
            df["j"] = pd.to_numeric(df["j"], errors="coerce").astype("Int64")
            df["pressure_psia"] = pd.to_numeric(df["pressure_psia"], errors="coerce")
            df["sw"] = pd.to_numeric(df["sw"], errors="coerce")
            df["sg"] = pd.to_numeric(df["sg"], errors="coerce")
            df = df.dropna(subset=["time_days", "layer", "i", "j"]).reset_index(drop=True)
            df["layer"] = df["layer"].astype(int)
            df["i"] = df["i"].astype(int)
            df["j"] = df["j"].astype(int)
        except Exception as exc:
            messagebox.showerror("Format error", f"Invalid numeric data in file:\n{exc}")
            return

        if df.empty:
            messagebox.showerror("No data", "The selected file has no valid map data.")
            return

        self.df = df
        self.file_var.set(os.path.basename(file_path))

        self.i_values = sorted(df["i"].unique().tolist())
        self.j_values = sorted(df["j"].unique().tolist())
        self.layer_values = sorted(df["layer"].unique().tolist())
        self.time_values = sorted(df["time_days"].unique().tolist())

        self.time_combo["values"] = [f"{t:g}" for t in self.time_values]
        self.time_index = len(self.time_values) - 1
        self.time_var.set(f"{self.time_values[-1]:g}")

        self.field_ranges = {}
        for field_name in MAP_FIELDS.values():
            valid_values = df[field_name].dropna()
            if valid_values.empty:
                self.field_ranges[field_name] = None
                continue
            field_min = float(valid_values.min())
            field_max = float(valid_values.max())
            self.field_ranges[field_name] = (field_min, field_max)

        self._on_direction_change()

    def _on_direction_change(self) -> None:
        direction = self.direction_var.get().upper()
        if direction == "XY":
            self.slice_label_var.set("Layer")
            self.slice_values = self.layer_values
        elif direction == "XZ":
            self.slice_label_var.set("J")
            self.slice_values = self.j_values
        else:
            self.slice_label_var.set("I")
            self.slice_values = self.i_values

        if not self.slice_values:
            self.slice_scale.configure(from_=1, to=1)
            self.slice_var.set(1)
            self.slice_info.configure(text="")
            self.update_plot()
            return

        self.slice_scale.configure(from_=self.slice_values[0], to=self.slice_values[-1])
        current = self.slice_var.get()
        if current not in self.slice_values:
            self.slice_var.set(self.slice_values[0])
        self._snap_slice_to_valid()
        self.update_plot()

    def _snap_slice_to_valid(self) -> None:
        if not self.slice_values:
            return
        selected = self.slice_var.get()
        nearest = min(self.slice_values, key=lambda value: abs(value - selected))
        self.slice_var.set(nearest)
        self.slice_info.configure(text=f"= {nearest}")

    def _selected_time(self) -> float | None:
        if not self.time_var.get():
            return None
        try:
            return float(self.time_var.get())
        except ValueError:
            return None

    def _sync_time_index_from_var(self) -> None:
        selected_time = self._selected_time()
        if selected_time is None or not self.time_values:
            return
        nearest = min(range(len(self.time_values)), key=lambda idx: abs(self.time_values[idx] - selected_time))
        self.time_index = nearest

    def _selected_time_mask(self, selected_time: float) -> tuple[float, pd.Series]:
        """Return canonical selected time and a tolerant boolean mask for dataframe filtering."""
        if not self.time_values:
            return selected_time, pd.Series(dtype=bool)

        canonical_time = self.time_values[self.time_index]
        tolerance = max(1.0e-9, abs(canonical_time) * 1.0e-9)
        mask = (self.df["time_days"] - canonical_time).abs() <= tolerance
        return canonical_time, mask

    def toggle_animation(self) -> None:
        if self.df is None or not self.time_values:
            messagebox.showinfo("No data", "Load a _maps.csv file before starting animation.")
            return

        if self.animating:
            self.stop_animation()
        else:
            self.start_animation()

    def start_animation(self) -> None:
        self._sync_time_index_from_var()
        self.animating = True
        self.animate_button.configure(text="Stop")
        self._animate_next_frame()

    def stop_animation(self) -> None:
        self.animating = False
        self.animate_button.configure(text="Animate")
        if self.animation_job is not None:
            self.after_cancel(self.animation_job)
            self.animation_job = None

    def _animate_next_frame(self) -> None:
        if not self.animating or not self.time_values:
            return

        if self.time_index >= len(self.time_values) - 1:
            self.stop_animation()
            return

        self.time_index += 1

        next_time = self.time_values[self.time_index]
        self.time_var.set(f"{next_time:g}")
        self.update_plot()
        self.animation_job = self.after(self._animation_delay_ms(), self._animate_next_frame)

    def _animation_delay_ms(self) -> int:
        try:
            seconds = float(self.speed_var.get())
        except ValueError:
            seconds = 0.5

        if seconds <= 0.0:
            seconds = 0.5

        return int(seconds * 1000)

    def _on_close(self) -> None:
        if self.closing:
            return
        self.closing = True

        self.after(250, lambda: os._exit(0))

        if self.animation_job is not None:
            self.after_cancel(self.animation_job)
            self.animation_job = None
        self.animating = False

        plt.close(self.figure)
        self.quit()
        self.destroy()

    def update_plot(self) -> None:
        self.ax.clear()
        self.current_grid = None
        self.hover_var.set("Hover over a block to see value")

        if self.df is None:
            self.ax.set_title("Open a _maps.csv file to start")
            self.canvas.draw_idle()
            return

        self._snap_slice_to_valid()
        selected_time = self._selected_time()
        if selected_time is None:
            self.ax.set_title("Select a valid time")
            self.canvas.draw_idle()
            return
        self._sync_time_index_from_var()
        selected_time, time_mask = self._selected_time_mask(selected_time)

        data_t = self.df[time_mask]
        if data_t.empty:
            self.ax.set_title("No data for selected time")
            self.canvas.draw_idle()
            return

        field_key = self.map_var.get().lower()
        field = MAP_FIELDS.get(field_key, "pressure_psia")
        direction = self.direction_var.get().upper()
        slice_value = self.slice_var.get()

        if direction == "XY":
            sliced = data_t[data_t["layer"] == slice_value]
            grid = sliced.pivot_table(index="j", columns="i", values=field, aggfunc="mean")
            xlabel = "i"
            ylabel = "j"
            title = f"{field_key.upper()} map - XY (layer={slice_value}, time={selected_time:g} d)"
            invert_vertical_axis = False
        elif direction == "XZ":
            sliced = data_t[data_t["j"] == slice_value]
            grid = sliced.pivot_table(index="layer", columns="i", values=field, aggfunc="mean")
            xlabel = "i"
            ylabel = "layer"
            title = f"{field_key.upper()} map - XZ (j={slice_value}, time={selected_time:g} d)"
            invert_vertical_axis = True
        else:
            sliced = data_t[data_t["i"] == slice_value]
            grid = sliced.pivot_table(index="layer", columns="j", values=field, aggfunc="mean")
            xlabel = "j"
            ylabel = "layer"
            title = f"{field_key.upper()} map - YZ (i={slice_value}, time={selected_time:g} d)"
            invert_vertical_axis = True

        if grid.empty:
            self.ax.set_title("No data for selected slice")
            self.canvas.draw_idle()
            return

        self.current_grid = grid
        self.current_field_key = field_key

        color_limits = self.field_ranges.get(field)
        if color_limits is not None:
            vmin, vmax = color_limits
            if vmin == vmax:
                vmin -= 0.5
                vmax += 0.5
        else:
            vmin = None
            vmax = None

        image = self.ax.imshow(
            grid.values,
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

        x_tick_positions = list(range(len(grid.columns)))
        y_tick_positions = list(range(len(grid.index)))
        self.ax.set_xticks(x_tick_positions)
        self.ax.set_xticklabels([str(value) for value in grid.columns], rotation=45, ha="right")
        self.ax.set_yticks(y_tick_positions)
        self.ax.set_yticklabels([str(value) for value in grid.index])

        if invert_vertical_axis:
            self.ax.invert_yaxis()

        if self.colorbar is None:
            self.colorbar = self.figure.colorbar(image, ax=self.ax, shrink=0.88)
        else:
            self.colorbar.update_normal(image)
        self.colorbar.set_label(field)

        self.canvas.draw_idle()


def main() -> None:
    initial_file = sys.argv[1] if len(sys.argv) > 1 else None
    app = MapViewerApp(initial_file=initial_file)
    app.mainloop()


if __name__ == "__main__":
    main()
