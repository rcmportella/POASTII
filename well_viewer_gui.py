import os
import tkinter as tk
from tkinter import filedialog, messagebox

import pandas as pd
import matplotlib.pyplot as plt


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
    df["label"] = (
        df["well_name"].fillna("")
        + " (id=" + df["well_id"].astype(str)
        + ", idx=" + df["well_index"].astype(str) + ")"
    )

    metrics = [
        ("qo_stb_d", "Oil rate (STB/D)"),
        ("qw_stb_d", "Water rate (STB/D)"),
        ("qg_mscf_d", "Gas rate (MSCF/D)"),
        ("bhp_psia", "BHP (PSIA)"),
    ]

    labels = sorted(df["label"].unique())
    nrows = len(metrics)
    fig, axes = plt.subplots(nrows=nrows, ncols=1, sharex=True, figsize=(10, 8))
    if nrows == 1:
        axes = [axes]

    for ax, (col, title) in zip(axes, metrics):
        for label in labels:
            sub = df[df["label"] == label]
            ax.plot(sub["time_days"], sub[col], label=label)
        ax.set_ylabel(title)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (days)")
    axes[0].legend(loc="best", fontsize=8)
    fig.suptitle(os.path.basename(file_path))
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


if __name__ == "__main__":
    main()
