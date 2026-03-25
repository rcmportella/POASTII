import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Dict, List

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.ticker import MaxNLocator


REQUIRED_COLS = {
    "time_days",
    "well_index",
    "well_id",
    "well_name",
    "qo_stb_d",
    "qw_stb_d",
    "qg_mscf_d",
    "bhp_psia",
}

METRICS = [
    ("qo_stb_d", "Oil rate (STB/D)"),
    ("qw_stb_d", "Water rate (STB/D)"),
    ("qg_mscf_d", "Gas rate (MSCF/D)"),
    ("rgo_mscf_stb", "RGO (MSCF/STB)"),
    ("wor", "WOR (STB/STB)"),
    ("bhp_psia", "BHP (PSIA)"),
]


def pick_files(parent: tk.Tk) -> List[str]:
    file_paths = filedialog.askopenfilenames(
        title="Select a _wells.csv file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        initialdir=os.getcwd(),
        parent=parent,
    )
    return list(file_paths)


def _build_case_name(file_path: str) -> str:
    base = os.path.basename(file_path)
    stem, _ = os.path.splitext(base)
    return stem


def prepare_case_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    missing = REQUIRED_COLS.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    case_name = _build_case_name(file_path)
    df = df.copy()
    oil_rate = df["qo_stb_d"].astype(float)
    safe_oil_rate = oil_rate.where(oil_rate != 0)
    df["rgo_mscf_stb"] = (df["qg_mscf_d"].astype(float) / safe_oil_rate)
    df["wor"] = (df["qw_stb_d"].astype(float) / safe_oil_rate)

    df["well_label"] = (
        df["well_name"].fillna("")
        + " (id=" + df["well_id"].astype(str)
        + ", idx=" + df["well_index"].astype(str) + ")"
    )
    df["case_name"] = case_name
    df["plot_label"] = "[" + case_name + "] " + df["well_label"]
    return df


class WellViewerApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.geometry("1100x800")
        self.root.title("Well Viewer")

        self.cases: List[Dict[str, object]] = []
        self.notebook: ttk.Notebook | None = None

        self._build_layout()
        self._render_empty_state()

    def _build_layout(self) -> None:
        top = ttk.Frame(self.root)
        top.pack(fill="x", padx=8, pady=(8, 4))

        ttk.Button(top, text="Add case(s)", command=self.add_cases).pack(side="left")
        ttk.Button(top, text="New session (clear all)", command=self.reset_session).pack(side="left", padx=(8, 0))

        self.case_info_var = tk.StringVar(value="No cases loaded")
        ttk.Label(top, textvariable=self.case_info_var).pack(side="left", padx=(12, 0))

        self.content_frame = ttk.Frame(self.root)
        self.content_frame.pack(fill="both", expand=True)

    def _update_header(self) -> None:
        if not self.cases:
            self.case_info_var.set("No cases loaded")
            self.root.title("Well Viewer")
            return

        names = [case["name"] for case in self.cases]
        self.case_info_var.set(f"Loaded cases ({len(names)}): " + ", ".join(names))
        self.root.title(f"Well Viewer - {len(names)} case(s)")

    def _clear_notebook(self) -> None:
        if self.notebook is not None:
            self.notebook.destroy()
            self.notebook = None

    def _render_empty_state(self) -> None:
        self._clear_notebook()
        empty = ttk.Frame(self.content_frame)
        empty.pack(fill="both", expand=True)
        ttk.Label(
            empty,
            text="Load one or more *_wells.csv files to compare well curves.",
            anchor="center",
        ).place(relx=0.5, rely=0.5, anchor="center")
        self.notebook = None
        self._update_header()

    def add_cases(self) -> None:
        file_paths = pick_files(self.root)
        if not file_paths:
            return

        loaded = 0
        failures: List[str] = []
        existing_paths = {case["path"] for case in self.cases}

        for file_path in file_paths:
            if file_path in existing_paths:
                continue
            try:
                df = prepare_case_dataframe(file_path)
            except Exception as exc:
                failures.append(f"{os.path.basename(file_path)}: {exc}")
                continue

            self.cases.append({
                "path": file_path,
                "name": _build_case_name(file_path),
                "df": df,
            })
            loaded += 1

        if loaded > 0:
            self.render_plots()

        if failures:
            messagebox.showwarning(
                "Some files were not loaded",
                "\n".join(failures),
                parent=self.root,
            )

    def reset_session(self) -> None:
        self.cases.clear()
        for child in self.content_frame.winfo_children():
            child.destroy()
        self._render_empty_state()

    def render_plots(self) -> None:
        for child in self.content_frame.winfo_children():
            child.destroy()

        self.notebook = ttk.Notebook(self.content_frame)
        self.notebook.pack(fill="both", expand=True)

        for col, title in METRICS:
            tab = ttk.Frame(self.notebook)
            self.notebook.add(tab, text=title)

            fig = plt.Figure(figsize=(10, 7), dpi=100)
            ax = fig.add_subplot(111)

            merged_time_values: List[float] = []
            merged_metric_values: List[float] = []

            for case in self.cases:
                df = case["df"]
                labels = sorted(df["plot_label"].unique())
                for label in labels:
                    sub = df[df["plot_label"] == label].sort_values("time_days")
                    ax.plot(sub["time_days"], sub[col], label=label)
                    merged_time_values.extend(pd.to_numeric(sub["time_days"], errors="coerce").tolist())
                    merged_metric_values.extend(pd.to_numeric(sub[col], errors="coerce").tolist())

            ax.set_title(title)
            ax.set_xlabel("Time (days)")
            ax.set_ylabel(title)
            ax.grid(True, alpha=0.3)

            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=8))
            ax.ticklabel_format(axis="both", style="plain", useOffset=False)
            ax.tick_params(axis="both", which="major", labelsize=9)

            ax.legend(loc="best", fontsize=8)
            fig.tight_layout()

            toolbar_frame = ttk.Frame(tab)
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas_widget = canvas.get_tk_widget()

            controls_frame = ttk.Frame(tab)
            controls_frame.pack(fill="x", padx=6, pady=(6, 2))

            ttk.Label(controls_frame, text="X scale:").pack(side="left")
            x_scale_var = tk.StringVar(value="linear")
            ttk.Combobox(
                controls_frame,
                textvariable=x_scale_var,
                values=("linear", "log"),
                width=8,
                state="readonly",
            ).pack(side="left", padx=(4, 10))

            ttk.Label(controls_frame, text="X min:").pack(side="left")
            x_min_var = tk.StringVar(value="")
            ttk.Entry(controls_frame, textvariable=x_min_var, width=10).pack(side="left", padx=(4, 10))

            ttk.Label(controls_frame, text="X max:").pack(side="left")
            x_max_var = tk.StringVar(value="")
            ttk.Entry(controls_frame, textvariable=x_max_var, width=10).pack(side="left", padx=(4, 10))

            ttk.Label(controls_frame, text="Y scale:").pack(side="left")
            y_scale_var = tk.StringVar(value="linear")
            ttk.Combobox(
                controls_frame,
                textvariable=y_scale_var,
                values=("linear", "log"),
                width=8,
                state="readonly",
            ).pack(side="left", padx=(4, 10))

            ttk.Label(controls_frame, text="Y min:").pack(side="left")
            y_min_var = tk.StringVar(value="")
            ttk.Entry(controls_frame, textvariable=y_min_var, width=10).pack(side="left", padx=(4, 10))

            ttk.Label(controls_frame, text="Y max:").pack(side="left")
            y_max_var = tk.StringVar(value="")
            ttk.Entry(controls_frame, textvariable=y_max_var, width=10).pack(side="left", padx=(4, 10))

            time_values = pd.to_numeric(pd.Series(merged_time_values), errors="coerce")
            positive_time_values = time_values[time_values > 0]
            metric_values = pd.to_numeric(pd.Series(merged_metric_values), errors="coerce")
            positive_metric_values = metric_values[metric_values > 0]

            def apply_x_axis(
                ax=ax,
                canvas=canvas,
                x_scale_var=x_scale_var,
                x_min_var=x_min_var,
                x_max_var=x_max_var,
                positive_time_values=positive_time_values,
            ) -> None:
                x_scale = x_scale_var.get()
                x_min_text = x_min_var.get().strip()
                x_max_text = x_max_var.get().strip()

                try:
                    x_min = float(x_min_text) if x_min_text else None
                    x_max = float(x_max_text) if x_max_text else None
                except ValueError:
                    messagebox.showerror("Invalid axis limits", "X min and X max must be numeric values.")
                    return

                if x_min is not None and x_max is not None and x_min >= x_max:
                    messagebox.showerror("Invalid axis limits", "X min must be smaller than X max.")
                    return

                if x_scale == "log":
                    if x_min is not None and x_min <= 0:
                        messagebox.showerror("Invalid log scale", "X min must be greater than zero for log scale.")
                        return
                    if x_max is not None and x_max <= 0:
                        messagebox.showerror("Invalid log scale", "X max must be greater than zero for log scale.")
                        return
                    if positive_time_values.empty:
                        messagebox.showerror("Invalid log scale", "No positive X values are available for log scale.")
                        return

                    auto_min = float(positive_time_values.min())
                    auto_max = float(positive_time_values.max())
                    left = x_min if x_min is not None else auto_min
                    right = x_max if x_max is not None else auto_max
                    if left >= right:
                        messagebox.showerror("Invalid axis limits", "X min must be smaller than X max.")
                        return

                    ax.set_xscale("log")
                    ax.set_xlim(left=left, right=right)
                else:
                    ax.set_xscale("linear")
                    if x_min is None and x_max is None:
                        ax.relim()
                        ax.autoscale_view(scalex=True, scaley=False)
                    else:
                        current_left, current_right = ax.get_xlim()
                        left = x_min if x_min is not None else current_left
                        right = x_max if x_max is not None else current_right
                        if left >= right:
                            messagebox.showerror("Invalid axis limits", "X min must be smaller than X max.")
                            return
                        ax.set_xlim(left=left, right=right)

                canvas.draw_idle()

            def apply_y_axis(
                ax=ax,
                canvas=canvas,
                y_scale_var=y_scale_var,
                y_min_var=y_min_var,
                y_max_var=y_max_var,
                positive_metric_values=positive_metric_values,
            ) -> None:
                y_scale = y_scale_var.get()
                y_min_text = y_min_var.get().strip()
                y_max_text = y_max_var.get().strip()

                try:
                    y_min = float(y_min_text) if y_min_text else None
                    y_max = float(y_max_text) if y_max_text else None
                except ValueError:
                    messagebox.showerror("Invalid axis limits", "Y min and Y max must be numeric values.")
                    return

                if y_min is not None and y_max is not None and y_min >= y_max:
                    messagebox.showerror("Invalid axis limits", "Y min must be smaller than Y max.")
                    return

                if y_scale == "log":
                    if y_min is not None and y_min <= 0:
                        messagebox.showerror("Invalid log scale", "Y min must be greater than zero for log scale.")
                        return
                    if y_max is not None and y_max <= 0:
                        messagebox.showerror("Invalid log scale", "Y max must be greater than zero for log scale.")
                        return
                    if positive_metric_values.empty:
                        messagebox.showerror("Invalid log scale", "No positive values are available for log scale.")
                        return

                    auto_min = float(positive_metric_values.min())
                    auto_max = float(positive_metric_values.max())
                    bottom = y_min if y_min is not None else auto_min
                    top = y_max if y_max is not None else auto_max
                    if bottom >= top:
                        messagebox.showerror("Invalid axis limits", "Y min must be smaller than Y max.")
                        return

                    ax.set_yscale("log")
                    ax.set_ylim(bottom=bottom, top=top)
                else:
                    ax.set_yscale("linear")
                    if y_min is None and y_max is None:
                        ax.relim()
                        ax.autoscale_view(scalex=False, scaley=True)
                    else:
                        current_bottom, current_top = ax.get_ylim()
                        bottom = y_min if y_min is not None else current_bottom
                        top = y_max if y_max is not None else current_top
                        if bottom >= top:
                            messagebox.showerror("Invalid axis limits", "Y min must be smaller than Y max.")
                            return
                        ax.set_ylim(bottom=bottom, top=top)

                canvas.draw_idle()

            def reset_axes(
                ax=ax,
                canvas=canvas,
                x_scale_var=x_scale_var,
                x_min_var=x_min_var,
                x_max_var=x_max_var,
                y_scale_var=y_scale_var,
                y_min_var=y_min_var,
                y_max_var=y_max_var,
            ) -> None:
                x_scale_var.set("linear")
                x_min_var.set("")
                x_max_var.set("")
                y_scale_var.set("linear")
                y_min_var.set("")
                y_max_var.set("")
                ax.set_xscale("linear")
                ax.set_yscale("linear")
                ax.relim()
                ax.autoscale_view()
                canvas.draw_idle()

            ttk.Button(controls_frame, text="Apply X axis", command=apply_x_axis).pack(side="left", padx=(4, 6))
            ttk.Button(controls_frame, text="Apply Y axis", command=apply_y_axis).pack(side="left", padx=(0, 6))
            ttk.Button(controls_frame, text="Reset", command=reset_axes).pack(side="left")

            toolbar_frame.pack(fill="x")
            canvas_widget.pack(fill="both", expand=True)

            toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
            toolbar.update()
            canvas.draw()

        self._update_header()


def main() -> None:
    root = tk.Tk()
    app = WellViewerApp(root)
    app.add_cases()

    root.mainloop()


if __name__ == "__main__":
    main()
