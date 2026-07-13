"""
C2-B GATE: converse-KAM on a SINGLE-resonance (integrable) Paul field must
reproduce the analytic island area (Paul-coordinate answer key).  Only after
this passes do we touch the four chaotic fields.
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import PaulField                                       # noqa: E402
from paul_island import island_area_gridcount, island_area_pendulum    # noqa: E402
from converse_kam import detection_map_uv, nonexistence_area_uv        # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def main(N=140, t_f=200.0, mode=(4, 2, 0.01)):
    f = PaulField([mode])
    S_grid = island_area_gridcount(f)
    S_pend = island_area_pendulum(f)
    rho = np.linspace(0.0, 1.0, N)
    th = np.linspace(0.0, 2.0 * np.pi, N)
    t0 = time.time()
    tc = detection_map_uv(f, rho, th, t_f=t_f, symmetry=False)
    dt = time.time() - t0
    S_ck = nonexistence_area_uv(tc, rho, th, t_f)
    np.savez(os.path.join(DATA, "paul_single_tc.npz"), tc=tc, rho=rho, th=th)
    print(f"single resonance {mode}:  N={N} t_f={t_f}  ({dt:.0f}s)")
    print(f"  answer key (gridcount) S_I = {S_grid:.4f}   (pendulum 8W = {S_pend:.4f})")
    print(f"  converse-KAM area      S   = {S_ck:.4f}   ratio = {S_ck/S_grid:.3f}")
    # timeout/undetected inside the island band (rho in [0.4,0.6])
    band = (rho >= 0.35) & (rho <= 0.65)
    det = np.isfinite(tc) & (tc <= t_f)
    print(f"  detected fraction overall = {det.mean():.3f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    m = np.where(np.isfinite(tc), tc, np.nan)
    im = ax.imshow(m, origin="lower", extent=[0, 2 * np.pi, 0, 1],
                   aspect="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label=r"$t_c$")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\psi=\rho$")
    ax.set_title(f"Paul single resonance {mode}: converse-KAM detection\n"
                 f"area {S_ck:.3f} vs analytic island {S_grid:.3f} (ratio {S_ck/S_grid:.3f})")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "paul_single_ckam.png"), dpi=140)
    plt.close(fig)
    print("wrote paul_single_ckam.png")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=140)
    ap.add_argument("--tf", type=float, default=200.0)
    args = ap.parse_args()
    main(N=args.N, t_f=args.tf)
