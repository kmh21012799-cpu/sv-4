"""
C3a P-2 GATE: validate WBA on the integrable single-resonance Paul field.
Expectation: high dig (regular) in island + circulating regions, low dig along
the separatrix (divergent period -> slow WBA convergence at finite T).
Cheap (single mode). Also prints the regular-vs-chaotic contrast on the real
m=4/m=36 fields (from point probes) for the record.
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import PaulField                                  # noqa: E402
from wba import dig_map                                           # noqa: E402
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def main(N=72, n_periods=300):
    f = PaulField([(4, 2, 0.01)])       # integrable single resonance, island ~[0.4,0.6]
    rho = np.linspace(0.0, 1.0, N)
    th = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    t0 = time.time()
    dig = dig_map(f, rho, th, n_periods=n_periods, which="dig_psi")
    dt = time.time() - t0
    np.savez(os.path.join(DATA, "wba_gate.npz"), dig=dig, rho=rho, th=th)
    print(f"single-resonance gate: N={N}, T={n_periods}, {dt:.0f}s")
    v = dig[~np.isnan(dig)]
    print(f"  dig_psi: median={np.median(v):.1f}  10th pct={np.percentile(v,10):.1f}"
          f"  90th={np.percentile(v,90):.1f}  frac(dig<5)={100*np.mean(v<5):.1f}%")

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(dig, origin="lower", extent=[0, 2 * np.pi, 0, 1],
                   aspect="auto", cmap="viridis", vmin=0, vmax=12)
    fig.colorbar(im, ax=ax, label="dig (WBA, h=$\\psi$)")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\rho$")
    ax.set_title(f"WBA gate: single-resonance Paul (integrable), T={n_periods}\n"
                 "high dig = regular (island + circulating); low = separatrix")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "wba_gate.png"), dpi=140)
    plt.close(fig)
    print("wrote wba_gate.png")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=72)
    ap.add_argument("--T", type=int, default=300)
    args = ap.parse_args()
    main(N=args.N, n_periods=args.T)
