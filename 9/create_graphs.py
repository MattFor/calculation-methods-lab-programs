import math
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

base = Path(__file__).resolve().parent
csv_path = base / "data.csv"

df = pd.read_csv(csv_path)

# --- wykres porównania rozwiązania ---
sol = df[df["typ"] == "solution"].copy()

plt.figure(figsize=(8, 5))
plt.plot(sol["x"], sol["exact"], label="analityczne")
plt.plot(sol["x"], sol["fd"], label="różnice skończone")
plt.plot(sol["x"], sol["shoot"], label="metoda strzałów")
plt.xlabel("x")
plt.ylabel("U(x)")
plt.title("Porównanie rozwiązania zagadnienia brzegowego")
plt.legend()
plt.tight_layout()
plt.savefig(base / "bvp_comparison.png", dpi=200)
plt.show()

# --- wykres zbieżności ---
conv = df[df["typ"] == "convergence"].copy()

plt.figure(figsize=(8, 5))
plt.plot([math.log10(h) for h in conv["h"]], [math.log10(e) for e in conv["err_fd"]], marker="o", label="różnice skończone")
plt.plot([math.log10(h) for h in conv["h"]], [math.log10(e) for e in conv["err_sh"]], marker="o", label="metoda strzałów")
plt.xlabel("log10(h)")
plt.ylabel("log10(maks. |błąd|)")
plt.title("Zbieżność metod numerycznych")
plt.legend()
plt.tight_layout()
plt.savefig(base / "bvp_convergence.png", dpi=200)
plt.show()
