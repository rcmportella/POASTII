#!/usr/bin/env python3
"""Read a text file interval, parse 4 numeric columns, run calculations, and write results.

Examples
--------
python relperm_conversion.py \
    --input SPE1CASE1.DATA \
    --output interval_results_case1.txt
"""

from __future__ import annotations

import argparse
import bisect
import re
from pathlib import Path

NUMBER_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
SAT_ROUND_DECIMALS = 8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read a line interval from a text file, extract the first 4 numeric columns, "
            "run calculations, and write a report."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input text file path (example: SPE1CASE1.DATA)",
    )
    parser.add_argument(
        "--output",
        default="interval_results.txt",
        help="Output text file path (default: interval_results.txt)",
    )
    parser.add_argument(
        "--start-line",
        type=int,
        default=144,
        help="Start line (1-based, inclusive)",
    )
    parser.add_argument(
        "--end-line",
        type=int,
        default=158,
        help="End line (1-based, inclusive)",
    )
    parser.add_argument(
        "--start-line-2",
        type=int,
        default=175,
        help="Second interval start line (1-based, inclusive)",
    )
    parser.add_argument(
        "--end-line-2",
        type=int,
        default=189,
        help="Second interval end line (1-based, inclusive)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if a line in the interval has fewer than 4 numeric values",
    )
    return parser.parse_args()


def prompt_line_interval(default_start: int, default_end: int, label: str) -> tuple[int, int]:
    start_text = input(f"{label} start line [{default_start}]: ").strip()
    end_text = input(f"{label} end line [{default_end}]: ").strip()

    start_line = int(start_text) if start_text else default_start
    end_line = int(end_text) if end_text else default_end
    return start_line, end_line


def to_float(token: str) -> float:
    # Fortran outputs often use D exponent (e.g. 1.0D+03), normalize to E.
    return float(token.replace("D", "E").replace("d", "e"))


def extract_four_columns(line: str) -> list[float]:
    matches = NUMBER_RE.findall(line)
    if len(matches) < 4:
        raise ValueError("Line has fewer than 4 numeric values")
    return [to_float(value) for value in matches[:4]]


def is_comment_line(line: str) -> bool:
    return line.lstrip().startswith("--")


def summarize(columns: list[list[float]], labels: list[str]) -> list[str]:
    lines: list[str] = []
    for label, values in zip(labels, columns):
        count = len(values)
        total = sum(values)
        mean = total / count if count else float("nan")
        min_v = min(values) if count else float("nan")
        max_v = max(values) if count else float("nan")
        lines.append(
            f"{label}: count={count} sum={total:.8g} mean={mean:.8g} min={min_v:.8g} max={max_v:.8g}"
        )
    return lines


def unique_xy_sorted(x_values: list[float], y_values: list[float]) -> tuple[list[float], list[float]]:
    """Sort x/y pairs by x and collapse duplicate x values using y average."""
    pairs = sorted(zip(x_values, y_values), key=lambda p: p[0])
    if not pairs:
        return [], []

    unique_x: list[float] = []
    unique_y: list[float] = []

    i = 0
    n = len(pairs)
    while i < n:
        x = pairs[i][0]
        total = 0.0
        count = 0
        while i < n and pairs[i][0] == x:
            total += pairs[i][1]
            count += 1
            i += 1
        unique_x.append(x)
        unique_y.append(total / count)

    return unique_x, unique_y


def interpolate_linear(x_known: list[float], y_known: list[float], x_new: list[float]) -> list[float]:
    """Piecewise-linear interpolation with edge clamping outside known domain."""
    if len(x_known) != len(y_known):
        raise ValueError("Interpolation arrays must have the same length")
    if not x_known:
        raise ValueError("Interpolation requires at least one known point")
    if len(x_known) == 1:
        return [y_known[0] for _ in x_new]

    result: list[float] = []
    for x in x_new:
        if x <= x_known[0]:
            result.append(y_known[0])
            continue
        if x >= x_known[-1]:
            result.append(y_known[-1])
            continue

        hi = bisect.bisect_right(x_known, x)
        lo = hi - 1
        x0, y0 = x_known[lo], y_known[lo]
        x1, y1 = x_known[hi], y_known[hi]

        if x1 == x0:
            result.append(y0)
            continue

        t = (x - x0) / (x1 - x0)
        result.append(y0 + t * (y1 - y0))

    return result


def read_interval_rows(
    all_lines: list[str], start_line: int, end_line: int, strict: bool
) -> tuple[list[tuple[int, list[float]]], list[tuple[int, str]], int]:
    start_idx = start_line - 1
    end_idx = min(end_line, len(all_lines))

    if start_idx >= len(all_lines):
        raise ValueError(f"Start line {start_line} is beyond end of file")

    rows: list[tuple[int, list[float]]] = []
    invalid_lines: list[tuple[int, str]] = []

    for line_no in range(start_idx + 1, end_idx + 1):
        raw_line = all_lines[line_no - 1]
        try:
            if is_comment_line(raw_line):
                raise ValueError("Comment line")
            values = extract_four_columns(raw_line)
            rows.append((line_no, values))
        except ValueError:
            invalid_lines.append((line_no, raw_line.strip()))
            if strict:
                raise ValueError(f"Line {line_no} does not contain 4 numeric values")

    return rows, invalid_lines, end_idx


def main() -> int:
    args = parse_args()

    # Always ask for interval so users can quickly adjust range each run.
    args.start_line, args.end_line = prompt_line_interval(args.start_line, args.end_line, "Interval 1")
    args.start_line_2, args.end_line_2 = prompt_line_interval(
        args.start_line_2, args.end_line_2, "Interval 2"
    )

    if args.start_line < 1 or args.end_line < args.start_line:
        raise ValueError("Line interval 1 is invalid. Use start >= 1 and end >= start.")
    if args.start_line_2 < 1 or args.end_line_2 < args.start_line_2:
        raise ValueError("Line interval 2 is invalid. Use start >= 1 and end >= start.")

    input_path = Path(args.input)
    output_path = Path(args.output)

    all_lines = input_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    rows, invalid_lines, end_idx = read_interval_rows(
        all_lines, args.start_line, args.end_line, args.strict
    )

    if not rows:
        raise ValueError("No valid lines with 4 numeric columns were found in interval 1")

    rows2, invalid_lines2, end_idx2 = read_interval_rows(
        all_lines, args.start_line_2, args.end_line_2, args.strict
    )

    if not rows2:
        raise ValueError("No valid lines with 4 numeric columns were found in interval 2")

    sw_values: list[float] = []
    so_values: list[float] = []
    krw_values: list[float] = []
    kro_values: list[float] = []
    cp_values: list[float] = []
    for line_no, values in rows:
        sw, krw, kro, cp = values
        so = 1.0 - sw
        sw_values.append(sw)
        so_values.append(so)
        krw_values.append(krw)
        kro_values.append(kro)
        cp_values.append(cp)

    sg_values: list[float] = []
    sog_values: list[float] = []
    krg_values: list[float] = []
    krog_values: list[float] = []
    cpog_values: list[float] = []
    for _line_no, values in rows2:
        sg, krg, krog, cpog = values
        sog = 1.0 - sg
        sg_values.append(sg)
        sog_values.append(sog)
        krg_values.append(krg)
        krog_values.append(krog)
        cpog_values.append(cpog)

    # Build a global merged saturation grid from all sets: sw/so/sg/sog.
    # Round to remove floating-point artifacts (e.g., 0.3 vs 0.30000000000000004).
    saturation_values = sorted({round(v, SAT_ROUND_DECIMALS) for v in (sw_values + so_values + sg_values + sog_values)})

    sw_x, krw_y = unique_xy_sorted(sw_values, krw_values)
    so_x, kro_y = unique_xy_sorted(so_values, kro_values)
    sg_x, krg_y = unique_xy_sorted(sg_values, krg_values)
    sog_x, krog_y = unique_xy_sorted(sog_values, krog_values)
    sw_x_cp, cp_y = unique_xy_sorted(sw_values, cp_values)
    sg_x_cpog, cpog_y = unique_xy_sorted(sg_values, cpog_values)

    kro_interp = interpolate_linear(so_x, kro_y, saturation_values)
    krw_interp = interpolate_linear(sw_x, krw_y, saturation_values)
    krg_interp = interpolate_linear(sg_x, krg_y, saturation_values)
    krog_interp = interpolate_linear(sog_x, krog_y, saturation_values)
    cp_interp = interpolate_linear(sw_x_cp, cp_y, saturation_values)
    cpog_interp = interpolate_linear(sg_x_cpog, cpog_y, saturation_values)

    calculated_rows: list[list[float]] = []
    for sat, kro_i, krw_i, krg_i, krog_i, cp_i, cpog_i in zip(
        saturation_values,
        kro_interp,
        krw_interp,
        krg_interp,
        krog_interp,
        cp_interp,
        cpog_interp,
    ):
        calculated_rows.append([sat, kro_i, krw_i, krg_i, krog_i, cp_i, cpog_i])

    original_values_only = [values for _, values in rows]
    col_values = [list(col) for col in zip(*original_values_only)]

    header = ["saturation", "kro", "krw", "krg", "krog", "cp", "cpog"]

    report_lines: list[str] = []
    report_lines.append(f"Input file: {input_path}")
    report_lines.append(f"Interval 1 (requested): {args.start_line}..{args.end_line}")
    report_lines.append(f"Interval 1 (processed): {args.start_line}..{end_idx}")
    report_lines.append(f"Interval 1 valid rows: {len(rows)}")
    report_lines.append(f"Interval 1 skipped rows: {len(invalid_lines)}")
    report_lines.append(f"Interval 2 (requested): {args.start_line_2}..{args.end_line_2}")
    report_lines.append(f"Interval 2 (processed): {args.start_line_2}..{end_idx2}")
    report_lines.append(f"Interval 2 valid rows: {len(rows2)}")
    report_lines.append(f"Interval 2 skipped rows: {len(invalid_lines2)}")

    report_lines.append("")
    report_lines.append("DATA (global merged saturation with interpolated kro, krw, krg, krog, cp, cpog)")
    report_lines.append(" ".join(f"{name:>14}" for name in header))

    for values in calculated_rows:
        formatted = [f"{v:>14.8g}" for v in values]
        report_lines.append(" ".join(formatted))

    report_lines.append("")
    report_lines.append("SUMMARY (original columns only)")
    report_lines.extend(summarize(col_values, ["sw", "krw", "kro", "cp"]))

    report_lines.append("")
    report_lines.append("SUMMARY (merged/interpolated columns)")
    sat_values = [values[0] for values in calculated_rows]
    kro_out = [values[1] for values in calculated_rows]
    krw_out = [values[2] for values in calculated_rows]
    krg_out = [values[3] for values in calculated_rows]
    krog_out = [values[4] for values in calculated_rows]
    cp_out = [values[5] for values in calculated_rows]
    cpog_out = [values[6] for values in calculated_rows]
    report_lines.extend(
        summarize(
            [sat_values, kro_out, krw_out, krg_out, krog_out, cp_out, cpog_out],
            ["saturation", "kro", "krw", "krg", "krog", "cp", "cpog"],
        )
    )

    if invalid_lines:
        report_lines.append("")
        report_lines.append("SKIPPED LINES (interval 1)")
        for line_no, text in invalid_lines:
            report_lines.append(f"line {line_no}: {text}")

    if invalid_lines2:
        report_lines.append("")
        report_lines.append("SKIPPED LINES (interval 2)")
        for line_no, text in invalid_lines2:
            report_lines.append(f"line {line_no}: {text}")

    output_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Wrote results to {output_path}")
    print(
        f"Processed interval1 rows={len(rows)} (skipped {len(invalid_lines)}), "
        f"interval2 rows={len(rows2)} (skipped {len(invalid_lines2)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
