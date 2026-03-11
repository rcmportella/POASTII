import os
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


HEADER_PATTERN = re.compile(r"\bSAT\b.*\bKRO\b.*\bKRW\b.*\bKRG\b", re.IGNORECASE)
REPEAT_PATTERN = re.compile(r"^(\d+)\*(.+)$")


class RelPermViewerApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("BOAST II Relative Permeability Viewer")
        self.geometry("980x760")
        self.closing = False

        self.file_var = tk.StringVar(value="No input file loaded")
        self.region_var = tk.StringVar(value="")

        self.tables: list[dict] = []
        self.fig1 = None
        self.fig2 = None
        self.ax1 = None
        self.ax2 = None
        self.canvas1: FigureCanvasTkAgg | None = None
        self.canvas2: FigureCanvasTkAgg | None = None

        self._build_ui()
        self.protocol("WM_DELETE_WINDOW", self._close_app)

    def _build_ui(self) -> None:
        frame = ttk.Frame(self, padding=12)
        frame.pack(fill=tk.X)

        ttk.Button(frame, text="Open BOASTII input", command=self.open_file).grid(
            row=0, column=0, padx=(0, 10), pady=(0, 8), sticky="w"
        )

        ttk.Label(frame, textvariable=self.file_var).grid(
            row=0, column=1, columnspan=3, sticky="w", pady=(0, 8)
        )

        ttk.Label(frame, text="Rock table:").grid(row=1, column=0, sticky="w")
        self.region_combo = ttk.Combobox(
            frame,
            textvariable=self.region_var,
            state="readonly",
            values=[],
            width=30,
        )
        self.region_combo.grid(row=1, column=1, padx=(0, 10), sticky="w")
        self.region_combo.bind("<<ComboboxSelected>>", lambda _: self.update_plots())

        ttk.Button(frame, text="Plot", command=self.update_plots).grid(row=1, column=2, sticky="w")

        ttk.Label(
            frame,
            text=(
                "Tab 1: KRO and KRW vs SW\n"
                "Tab 2: KRG and KRW vs SW\n"
                "Axes are fixed from 0 to 1."
            ),
        ).grid(row=2, column=0, columnspan=4, pady=(12, 0), sticky="w")

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=12, pady=(6, 12))

        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab1, text="KRO & KRW vs SW")
        self.notebook.add(self.tab2, text="KRG & KRW vs SW")

        self.fig1, self.ax1 = plt.subplots(figsize=(8.6, 6.0), dpi=100)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.tab1)
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar1 = NavigationToolbar2Tk(self.canvas1, self.tab1)
        toolbar1.update()

        self.fig2, self.ax2 = plt.subplots(figsize=(8.6, 6.0), dpi=100)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.tab2)
        self.canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar2 = NavigationToolbar2Tk(self.canvas2, self.tab2)
        toolbar2.update()

        self._draw_empty_plots()

    def open_file(self) -> None:
        file_path = filedialog.askopenfilename(
            title="Select BOASTII input file",
            filetypes=[("Data files", "*.DAT *.dat *.txt"), ("All files", "*.*")],
            initialdir=os.getcwd(),
        )
        if not file_path:
            return

        try:
            tables = self.parse_relperm_tables(file_path)
        except Exception as exc:
            messagebox.showerror("Parse error", f"Could not parse file:\n{exc}")
            return

        if not tables:
            messagebox.showerror(
                "No relative permeability data",
                "Could not find SAT/KRO/KRW/KRG table in selected file.",
            )
            return

        self.tables = tables
        self.file_var.set(os.path.basename(file_path))

        labels = [f"Table {idx + 1} (line {tbl['line_start']})" for idx, tbl in enumerate(tables)]
        self.region_combo["values"] = labels
        self.region_var.set(labels[0])
        self.update_plots()

    @staticmethod
    def _expand_fortran_repeats(tokens: list[str]) -> list[str]:
        expanded: list[str] = []
        for token in tokens:
            token = token.strip()
            if not token:
                continue

            match = REPEAT_PATTERN.match(token)
            if match:
                repeat = int(match.group(1))
                value = match.group(2).strip()
                expanded.extend([value] * repeat)
            else:
                expanded.append(token)

        return expanded

    @staticmethod
    def _to_float(token: str) -> float:
        token = token.replace("D", "E").replace("d", "e")
        return float(token)

    def parse_relperm_tables(self, file_path: str) -> list[dict]:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        tables: list[dict] = []
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not HEADER_PATTERN.search(line):
                i += 1
                continue

            sat_vals: list[float] = []
            kro_vals: list[float] = []
            krw_vals: list[float] = []
            krg_vals: list[float] = []

            j = i + 1
            while j < len(lines):
                raw = lines[j].strip()
                if not raw:
                    break

                first = raw.split()[0]
                if first.upper().startswith("ITHREE"):
                    break

                tokens = self._expand_fortran_repeats(raw.split())
                floats: list[float] = []
                for tok in tokens:
                    try:
                        floats.append(self._to_float(tok))
                    except ValueError:
                        # Stop row parse once a non-numeric token appears.
                        break

                if len(floats) < 4:
                    break

                sat_vals.append(floats[0])
                kro_vals.append(floats[1])
                krw_vals.append(floats[2])
                krg_vals.append(floats[3])
                j += 1

            if sat_vals:
                tables.append(
                    {
                        "line_start": i + 1,
                        "sat": sat_vals,
                        "kro": kro_vals,
                        "krw": krw_vals,
                        "krg": krg_vals,
                    }
                )

            i = j + 1

        return tables

    def _draw_empty_plots(self) -> None:
        self.ax1.clear()
        self.ax1.set_xlim(0.0, 1.0)
        self.ax1.set_ylim(0.0, 1.0)
        self.ax1.set_xlabel("SW")
        self.ax1.set_ylabel("Relative Permeability")
        self.ax1.set_title("KRO and KRW vs SW")
        self.ax1.grid(True, alpha=0.3)
        self.canvas1.draw_idle()

        self.ax2.clear()
        self.ax2.set_xlim(0.0, 1.0)
        self.ax2.set_ylim(0.0, 1.0)
        self.ax2.set_xlabel("SW")
        self.ax2.set_ylabel("Relative Permeability")
        self.ax2.set_title("KRG and KRW vs SW")
        self.ax2.grid(True, alpha=0.3)
        self.canvas2.draw_idle()

    def _selected_table(self) -> dict | None:
        if not self.tables:
            return None

        labels = list(self.region_combo["values"])
        if not labels:
            return self.tables[0]

        try:
            idx = labels.index(self.region_var.get())
        except ValueError:
            idx = 0

        return self.tables[idx]

    def update_plots(self) -> None:
        table = self._selected_table()
        if table is None:
            messagebox.showinfo("No data", "Load a BOASTII input file first.")
            return

        sat = table["sat"]
        kro = table["kro"]
        krw = table["krw"]
        krg = table["krg"]

        # Table saturation is phase saturation: So for KRO, Sw for KRW, and Sg for KRG.
        # Convert KRO and KRG x-axis to equivalent Sw scale to show both requested plots vs Sw.
        sw_for_kro = [1.0 - s for s in sat]
        sw_for_krw = sat
        sw_for_krg = [1.0 - s for s in sat]

        self.ax1.clear()
        self.ax1.plot(sw_for_kro, kro, marker="o", label="KRO vs SW (SW=1-SO)")
        self.ax1.plot(sw_for_krw, krw, marker="s", label="KRW vs SW")
        self.ax1.set_xlim(0.0, 1.0)
        self.ax1.set_ylim(0.0, 1.0)
        self.ax1.set_xlabel("SW")
        self.ax1.set_ylabel("Relative Permeability")
        self.ax1.set_title("KRO and KRW vs SW")
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend(loc="best")
        self.canvas1.draw_idle()

        self.ax2.clear()
        self.ax2.plot(sw_for_krg, krg, marker="^", label="KRG vs SW (SW=1-SG)")
        self.ax2.plot(sw_for_krw, krw, marker="s", label="KRW vs SW")
        self.ax2.set_xlim(0.0, 1.0)
        self.ax2.set_ylim(0.0, 1.0)
        self.ax2.set_xlabel("SW")
        self.ax2.set_ylabel("Relative Permeability")
        self.ax2.set_title("KRG and KRW vs SW")
        self.ax2.grid(True, alpha=0.3)
        self.ax2.legend(loc="best")
        self.canvas2.draw_idle()

    def _close_app(self) -> None:
        if self.closing:
            return
        self.closing = True

        if self.fig1 is not None:
            plt.close(self.fig1)
            self.fig1 = None
        if self.fig2 is not None:
            plt.close(self.fig2)
            self.fig2 = None
        plt.close("all")
        self.quit()
        self.destroy()


def main() -> None:
    app = RelPermViewerApp()
    app.mainloop()


if __name__ == "__main__":
    main()
