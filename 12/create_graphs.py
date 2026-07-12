import csv
from pathlib import Path

import matplotlib.pyplot as plt

"""
Czyta dane z pliku CSV i zwraca kolumny jako listy
"""


def load_main_data(path: str):
	xs = []
	exact = []
	equidistant = []
	chebyshev = []

	with open(path, newline="", encoding="utf-8") as file:
		reader = csv.DictReader(file)
		for row in reader:
			xs.append(float(row["x"]))
			exact.append(float(row["exact"]))
			equidistant.append(float(row["equidistant"]))
			chebyshev.append(float(row["chebyshev"]))

	return xs, exact, equidistant, chebyshev


"""
Czyta węzły z pliku CSV i zwraca punkty do zaznaczenia na wykresie
"""


def load_nodes(path: str):
	xs = []
	ys = []

	with open(path, newline="", encoding="utf-8") as file:
		reader = csv.DictReader(file)
		for row in reader:
			xs.append(float(row["x"]))
			ys.append(float(row["f"]))

	return xs, ys


"""
Tworzy wykres porównujący funkcję dokładną i oba wielomiany interpolacyjne
"""


def make_plots():
	data_path = Path("../output/csv/runge_data.csv")
	eq_nodes_path = Path("../output/csv/equidistant_nodes.csv")
	ch_nodes_path = Path("../output/csv/chebyshev_nodes.csv")

	if not data_path.exists():
		raise FileNotFoundError("Brak pliku runge_data.csv")

	xs, exact, equidistant, chebyshev = load_main_data(str(data_path))
	eq_x, eq_y = load_nodes(str(eq_nodes_path))
	ch_x, ch_y = load_nodes(str(ch_nodes_path))

	plt.figure(figsize=(11, 6))
	plt.plot(xs, exact, label="f(x) = x / (1 + a x^4)")
	plt.plot(xs, equidistant, label="Interpolacja na węzłach równoodległych")
	plt.plot(xs, chebyshev, label="Interpolacja na węzłach Czebyszewa")
	plt.scatter(eq_x, eq_y, s=20, label="Węzły równoodległe")
	plt.scatter(ch_x, ch_y, s=20, label="Węzły Czebyszewa")
	plt.title("Zjawisko Rungego w interpolacji wielomianowej Newtona")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.grid(True)
	plt.legend()
	plt.tight_layout()
	plt.savefig("../output/runge_full.png", dpi=200)

	plt.figure(figsize=(11, 6))
	plt.plot(xs, exact, label="f(x) = x / (1 + a x^4)")
	plt.plot(xs, equidistant, label="Interpolacja na węzłach równoodległych")
	plt.plot(xs, chebyshev, label="Interpolacja na węzłach Czebyszewa")
	plt.scatter(eq_x, eq_y, s=20, label="Węzły równoodległe")
	plt.scatter(ch_x, ch_y, s=20, label="Węzły Czebyszewa")
	plt.xlim(-1.0, -0.6)
	plt.title("Powiększenie lewego brzegu przedziału")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.grid(True)
	plt.legend()
	plt.tight_layout()
	plt.savefig("../output/runge_zoom_left.png", dpi=200)


if __name__ == "__main__":
	make_plots()
