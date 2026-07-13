"""Cleaner gate: dig map of the m=4 critical field (real chaos) — island=high,
chaotic sea=low. Complements the single-resonance point probes."""
import os, sys, time
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4
from wba import dig_map_rk4
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
f = paul_m4(); N = 80
rho = np.linspace(0, 1, N); th = np.linspace(0, 2*np.pi, N, endpoint=False)
t0 = time.time()
dig = dig_map_rk4(f, rho, th, n_periods=1000, steps_per_period=32, which="dig_psi")
print(f"m4 dig map N={N} M=32 in {time.time()-t0:.0f}s")
fig, ax = plt.subplots(figsize=(7,5))
im = ax.imshow(dig, origin="lower", extent=[0,2*np.pi,0,1], aspect="auto",
               cmap="viridis", vmin=0, vmax=12)
fig.colorbar(im, ax=ax, label="dig (WBA, h=$\\psi$)")
for r in (0.125,0.875): ax.axhline(r, color="w", lw=0.5, ls=":")
ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\rho$")
ax.set_title("WBA gate (m=4 critical field): island chains = high dig,\nchaotic sea = low dig")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"wba_gate_m4.png"), dpi=140)
print("wrote wba_gate_m4.png")
