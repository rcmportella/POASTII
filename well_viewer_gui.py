import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

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


def pick_file() -> str:
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select a _wells.csv file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        initialdir=os.getcwd(),
    )
    root.destroy()
    return file_path


def main() -> None:
    file_path = pick_file()
    if not file_path:
        print("No file selected.")
        return

    df = pd.read_csv(file_path)
    missing = REQUIRED_COLS.difference(df.columns)
    if missing:
        messagebox.showerror("Missing columns", f"Missing columns: {sorted(missing)}")
        return

    # Build labels
    df = df.copy()
    oil_rate = df["qo_stb_d"].astype(float)
    safe_oil_rate = oil_rate.where(oil_rate != 0)
    df["rgo_mscf_stb"] = (df["qg_mscf_d"].astype(float) / safe_oil_rate)
    df["wor"] = (df["qw_stb_d"].astype(float) / safe_oil_rate)

    df["label"] = (
        df["well_name"].fillna("")
        + " (id=" + df["well_id"].astype(str)
        + ", idx=" + df["well_index"].astype(str) + ")"
    )

    metrics = [
        ("qo_stb_d", "Oil rate (STB/D)"),
        ("qw_stb_d", "Water rate (STB/D)"),
        ("qg_mscf_d", "Gas rate (MSCF/D)"),
        ("rgo_mscf_stb", "RGO (MSCF/STB)"),
        ("wor", "WOR (STB/STB)"),
        ("bhp_psia", "BHP (PSIA)"),
    ]

    labels = sorted(df["label"].unique())

    root = tk.Tk()
    root.title(f"Well Viewer - {os.path.basename(file_path)}")
    root.geometry("1100x800")

    notebook = ttk.Notebook(root)
    notebook.pack(fill="both", expand=True)

    for col, title in metrics:
        tab = ttk.Frame(notebook)
        notebook.add(tab, text=title)

        fig = plt.Figure(figsize=(10, 7), dpi=100)
        ax = fig.add_subplot(111)

        for label in labels:
            sub = df[df["label"] == label].sort_values("time_days")
            ax.plot(sub["time_days"], sub[col], label=label)

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

        time_values = pd.to_numeric(df["time_days"], errors="coerce")
        positive_time_values = time_values[time_values > 0]
        metric_values = pd.to_numeric(df[col], errors="coerce")
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

    root.mainloop()


if __name__ == "__main__":
    main()
