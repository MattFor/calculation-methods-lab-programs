from __future__ import annotations

import csv
import math
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

import matplotlib.pyplot as plt

METHODS = [
	("forward_2", "2-point forward"),
	("forward_3", "3-point forward"),
	("backward_2", "2-point backward"),
	("backward_3", "3-point backward"),
	("central_2", "2-point central"),
]

POINTS = [
	("0", r"$x=0$"),
	("pi_over_4", r"$x=\pi/4$"),
	("pi_over_2", r"$x=\pi/2$"),
]

TYPES = [
	("double", "double", "solid"),
	("long_double", "long double", "dashed"),
]


@dataclass
class Series:
	h: list[float]
	errors: dict[str, list[float]]


def read_csv(path: Path) -> Series:
	h: list[float] = []
	errors: dict[str, list[float]] = {key: [] for key, _ in METHODS}

	with path.open("r", newline="", encoding="utf-8-sig") as f:
		reader = csv.DictReader(f)
		for row in reader:
			if not row:
				continue
			h_val = float(row["h"])
			h.append(h_val)
			for key, _ in METHODS:
				errors[key].append(float(row[key]))

	return Series(h=h, errors=errors)


def machine_threshold_index(errors: list[float]) -> int:
	if len(errors) < 3:
		return len(errors)

	best = 0
	for i in range(1, len(errors)):
		if errors[i] < errors[best]:
			best = i
	return max(2, best + 1)


def fit_slope_loglog(h: list[float], errors: list[float]) -> Optional[float]:
	xs: list[float] = []
	ys: list[float] = []

	stop = machine_threshold_index(errors)
	for i in range(min(stop, len(h))):
		if h[i] > 0.0 and errors[i] > 0.0 and math.isfinite(h[i]) and math.isfinite(errors[i]):
			xs.append(math.log10(h[i]))
			ys.append(math.log10(errors[i]))

	if len(xs) < 2:
		return None

	n = float(len(xs))
	sx = sum(xs)
	sy = sum(ys)
	sxx = sum(x * x for x in xs)
	sxy = sum(x * y for x, y in zip(xs, ys))

	denom = n * sxx - sx * sx
	if abs(denom) < 1e-30:
		return None

	return (n * sxy - sx * sy) / denom


def find_input_file(base_dir: Path, type_name: str, point_key: str) -> Path:
	return base_dir / f"{type_name}_{point_key}.csv"


def plot_one_panel(ax, base_dir: Path, point_key: str, point_label: str) -> None:
	method_colors = {
		"forward_2": "C0",
		"forward_3": "C1",
		"backward_2": "C2",
		"backward_3": "C3",
		"central_2": "C4",
	}

	for type_name, type_label, linestyle in TYPES:
		path = find_input_file(base_dir, type_name, point_key)
		if not path.exists():
			ax.text(0.5, 0.5, f"Missing:\n{path.name}", ha="center", va="center", transform=ax.transAxes)
			continue

		series = read_csv(path)

		for key, method_label in METHODS:
			color = method_colors[key]

			# Keep only finite values and the same length as x/h
			xs: list[float] = []
			ys: list[float] = []

			for hv, ev in zip(series.h, series.errors[key]):
				if hv > 0.0 and ev > 0.0 and math.isfinite(hv) and math.isfinite(ev):
					xs.append(math.log10(hv))
					ys.append(math.log10(ev))

			label = f"{type_label}: {method_label}"
			ax.plot(xs, ys, marker="o", markersize=3.5, linewidth=1.2, linestyle=linestyle, color=color, label=label)

			slope = fit_slope_loglog(series.h, series.errors[key])
			if slope is not None:
				# Annotate only once per method/type pair near the first point
				ax.annotate(
					f"{slope:.2f}",
					xy=(xs[0], ys[0]),
					xytext=(4, 4),
					textcoords="offset points",
					fontsize=7,
					color=color,
					alpha=0.8,
				)

	ax.set_title(point_label)
	ax.set_xlabel(r"$\log_{10}(h)$")
	ax.grid(True, which="both", linestyle=":", linewidth=0.7, alpha=0.8)


def main() -> int:
	parser = argparse.ArgumentParser(description="Create a single comparison graph from CSV files.")
	parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent.parent, help="Directory with CSV files (default: data)")
	parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent.parent / "output" / "graphs.png", help="Output image path (default: data/graphs.png)")
	parser.add_argument("--show", action="store_true", help="Show the figure interactively")
	args = parser.parse_args()

	graph, axes = plt.subplots(
		1, 3,
		figsize=(19, 9.5),
		sharey=True
	)

	graph.subplots_adjust(
		top=0.80,
		bottom=0.12,
		left=0.05,
		right=0.98,
		wspace=0.04
	)

	for ax, (point_key, point_label) in zip(axes, POINTS):
		plot_one_panel(ax, args.data_dir, point_key, point_label)

	axes[0].set_ylabel(r"$\log_{10}(|\mathrm{error}|)$")

	handles, labels = axes[0].get_legend_handles_labels()
	seen: set[str] = set()
	filtered_handles = []
	filtered_labels = []
	for handle, label in zip(handles, labels):
		if label not in seen:
			seen.add(label)
			filtered_handles.append(handle)
			filtered_labels.append(label)

	graph.legend(
		filtered_handles,
		filtered_labels,
		loc="lower center",
		ncol=2,
		fontsize=8,
		frameon=True,
		bbox_to_anchor=(0.5, -0.02)
	)

	graph.subplots_adjust(bottom=0.13)

	graph.suptitle(
		"Absolute error of finite-difference approximations for f(x) = cos(x)",
		y=0.87,
		fontsize=30
	)

	args.output.parent.mkdir(parents=True, exist_ok=True)
	graph.savefig(args.output, dpi=500, bbox_inches="tight")
	print(f"Saved graphs: {args.output}")

	if args.show:
		plt.show()

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
