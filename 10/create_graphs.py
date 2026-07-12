import csv
import numpy as np
import matplotlib.pyplot as plt

path = "../output/"


def read_solution(path):
	t, exact, numerical = [], [], []

	with open(path, newline="") as f:
		reader = csv.DictReader(f)
		for row in reader:
			t.append(float(row["t"]))
			exact.append(float(row["exact"]))
			numerical.append(float(row["numerical"]))

	return np.array(t), np.array(exact), np.array(numerical)


def plot_solution(csv_file, png_file, title):
	t, ye, yn = read_solution(csv_file)

	plt.figure(figsize=(8, 5))
	plt.plot(t, ye, "-", linewidth=2, label="analityczne")
	plt.plot(t, yn, "o", markersize=3, label="numeryczne")
	plt.xlabel("t")
	plt.ylabel("y(t)")
	plt.title(title)
	plt.grid(True, alpha=0.3)
	plt.legend()
	plt.tight_layout()
	plt.savefig(path + png_file, dpi=200)
	plt.close()


def read_convergence(path):
	h, e1, e2, e3 = [], [], [], []

	with open(path, newline="") as f:
		reader = csv.DictReader(f)
		for row in reader:
			h.append(float(row["h"]))
			e1.append(float(row["explicit_euler"]))
			e2.append(float(row["implicit_euler"]))
			e3.append(float(row["trapezoid"]))

	return np.array(h), np.array(e1), np.array(e2), np.array(e3)


def order_fit(h, err, mask=None):
	if mask is None:
		mask = np.ones_like(h, dtype=bool)

	x = np.log10(h[mask])
	y = np.log10(err[mask])
	slope, intercept = np.polyfit(x, y, 1)

	return slope, intercept


def plot_convergence(csv_file, png_file):
	h, e1, e2, e3 = read_convergence(csv_file)

	plt.figure(figsize=(8, 5))
	plt.plot(np.log10(h), np.log10(e1), "o-", markersize=3, label="Euler jawny")
	plt.plot(np.log10(h), np.log10(e2), "o-", markersize=3, label="Euler pośredni")
	plt.plot(np.log10(h), np.log10(e3), "o-", markersize=3, label="Metoda trapezów")
	plt.xlabel("log10(dt)")
	plt.ylabel("log10(max |błędu|)")
	plt.title("Zależność błędu od kroku czasowego")
	plt.grid(True, alpha=0.3)
	plt.legend()
	plt.tight_layout()
	plt.savefig(path + png_file, dpi=200)
	plt.close()

	# Oszacowanie rzędów z najmniejszych kroków
	mask = h <= 2 ** -10
	for name, err in [("Euler jawny", e1), ("Euler pośredni", e2), ("Metoda trapezów", e3)]:
		p, _ = order_fit(h, err, mask=mask)
		print(f"{name}: oszacowany rząd ~ {p:.3f}")


if __name__ == "__main__":
	plot_solution("../output/csv/explicit_stable.csv", "explicit_stable.png", "Euler jawny — przypadek stabilny (h = 0.01)")
	plot_solution("../output/csv/explicit_border.csv", "explicit_border.png", "Euler jawny — przypadek graniczny (h = 0.02)")
	plot_solution("../output/csv/explicit_unstable.csv", "explicit_unstable.png", "Euler jawny — przypadek niestabilny (h = 0.03)")
	plot_solution("../output/csv/implicit_euler.csv", "implicit_euler.png", "Euler pośredni (h = 0.02)")
	plot_solution("../output/csv/trapezoid.csv", "trapezoid.png", "Metoda trapezów (h = 0.02)")
	plot_convergence("../output/csv/convergence.csv", "convergence.png")
