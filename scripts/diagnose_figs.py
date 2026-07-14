"""Two diagnostic figures for the C2 checkpoint (data already computed)."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1, KMMField                     # noqa: E402
from kmm_island import island_area, _fixed_points_21         # noqa: E402
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")


def fig_grid_noise():
    f = example1(0.004); S_I = island_area(f)
    psiO, psiX = _fixed_points_21(f); PsiX = f.invariant_single(psiX, np.pi / 2, 0.0)

    def area(N, L=0.9, noff=1):
        dy = 2 * L / (N - 1); tot = 0.0; cnt = 0
        for a in range(noff):
            for b in range(noff):
                yg = np.linspace(-L, L, N) + (a / noff) * dy
                zg = np.linspace(-L, L, N) + (b / noff) * dy
                Y, Z = np.meshgrid(yg, zg)
                psi, th = KMMField.from_symplectic(Y, Z, B0=f.B0)
                inside = (psi > 1e-4) & (f.invariant_single(psi, th, 0.0) < PsiX)
                tot += inside.sum() * dy * dy; cnt += 1
        return tot / cnt

    Ns = list(range(80, 241, 5))
    e1 = [100 * (area(N, noff=1) - S_I) / S_I for N in Ns]
    e4 = [100 * (area(N, noff=4) - S_I) / S_I for N in Ns]
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.plot(Ns, e1, "o-", ms=3, label="naive cell-count (pure discretisation)")
    ax.plot(Ns, e4, "s-", ms=3, label="4x4 offset-averaged")
    ax.axhline(0, color="gray", lw=0.6)
    ax.axhspan(-0.1, 0.1, color="green", alpha=0.08, label="+-0.1% band")
    ax.set_xlabel("grid N"); ax.set_ylabel("area error vs $S_I$  (%)")
    ax.set_title("Grid noise is pure discretisation (no converse-KAM here):\n"
                 "true-island cell-count wobbles +-0.5%, offset-averaging removes it")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "diag_grid_noise.png"), dpi=140)
    plt.close(fig); print("wrote diag_grid_noise.png")


def fig_m36_span():
    # values measured in diagnose_c2.py (recorded here to avoid a 15-min re-run)
    turns = np.array([1000, 10000, 100000])
    width = np.array([0.186, 0.533, 0.756])
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.semilogx(turns, width, "o-", ms=6)
    ax.axhline(0.875 - 0.125, color="r", ls="--", label="full span 7/8-1/8 = 0.75")
    ax.set_xlabel("toroidal turns (single m=36 orbit)")
    ax.set_ylabel(r"cumulative $\psi$-span width")
    ax.set_title("m=36 chaos is CONNECTED but SLOW (cantori):\n"
                 "one orbit's span grows to [1/8,7/8] given enough time")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "diag_m36_span.png"), dpi=140)
    plt.close(fig); print("wrote diag_m36_span.png")


if __name__ == "__main__":
    fig_grid_noise()
    fig_m36_span()
