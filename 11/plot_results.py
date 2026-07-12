from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
results_dir = BASE_DIR / "results"

convergence = pd.read_csv(results_dir / "convergence.csv")
time_error = pd.read_csv(results_dir / "time_error.csv")
snapshots = pd.read_csv(results_dir / "snapshots.csv")


def savefig(path: Path):
	plt.tight_layout()
	plt.savefig(path, dpi=200)
	plt.close()


METHOD_NAMES = {
	"explicit_ftcs": "Klasyczna Metoda Bezpośrednia",
	"crank_nicolson_thomas": "Crank–Nicolson (Thomas)",
	"crank_nicolson_jacobi": "Crank–Nicolson (Jacobi)",
}


# 1. Wykres zbieżności: maksymalny błąd dla t_max w funkcji h
plt.figure(figsize=(8, 5))

for method, group in convergence.groupby("method"):
	group = group.sort_values("h")

	plt.loglog(
		group["h"],
		group["max_abs_error_at_t_max"],
		marker="o",
		label=METHOD_NAMES.get(method, method),
	)

plt.xlabel("Krok przestrzenny h")
plt.ylabel("Maksymalny błąd bezwzględny dla t_max")
plt.title("Badanie zbieżności metod")
plt.grid(True, which="both", linestyle=":")
plt.legend()

savefig(results_dir / "convergence.png")


# 2. Błąd w funkcji czasu
plt.figure(figsize=(8, 5))

for method, group in time_error.groupby("method"):
	group = group.sort_values("time")

	plt.plot(
		group["time"],
		group["max_abs_error"],
		label=METHOD_NAMES.get(method, method),
	)

plt.xlabel("Czas t")
plt.ylabel("Maksymalny błąd bezwzględny")
plt.title("Zmiana błędu w czasie")
plt.grid(True, linestyle=":")
plt.legend()

savefig(results_dir / "time_error.png")


# 3. Snapshoty z rozwiązaniem analitycznym + numerycznym
snapshot_ids = sorted(snapshots["snapshot_id"].unique())

for snapshot_id in snapshot_ids:
	subset = snapshots[snapshots["snapshot_id"] == snapshot_id].copy()
	time_value = subset["time"].iloc[0]

	plt.figure(figsize=(8, 5))

	reference = (
		subset.groupby("x", as_index=False)["exact"]
		.mean()
		.sort_values("x")
	)

	plt.plot(
		reference["x"],
		reference["exact"],
		linewidth=2,
		label="Rozwiązanie analityczne",
	)

	for method, group in subset.groupby("method"):
		group = group.sort_values("x")

		plt.scatter(
			group["x"],
			group["numerical"],
			s=14,
			label=METHOD_NAMES.get(method, method),
		)

	plt.xlabel("Położenie x")
	plt.ylabel("U(x,t)")
	plt.title(f"Porównanie rozwiązań dla t = {time_value:.6f}")
	plt.grid(True, linestyle=":")
	plt.legend()

	savefig(results_dir / f"snapshot_{snapshot_id:02d}.png")


print("Wykresy zostały zapisane w katalogu results/")