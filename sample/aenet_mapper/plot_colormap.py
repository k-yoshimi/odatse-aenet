"""Plot distance-energy curve from mapper ColorMap.txt output."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = np.loadtxt("output/ColorMap.txt")
z = data[:, 0]
fx = data[:, 1]

plt.figure(figsize=(6, 4))
plt.plot(z, fx, "o-", markersize=4)
plt.xlabel("N-N distance (Angstrom)")
plt.ylabel(r"Energy (eV/$N_{\rm atom}$)")
plt.title("N2 dimer: Energy vs N-N distance")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("distance_energy.png", dpi=150)
print("Saved: distance_energy.png")

# Print minimum
imin = np.argmin(fx)
print(f"Energy minimum: z = {z[imin]:.4f} Angstrom, E = {fx[imin]:.6f} eV/atom")
