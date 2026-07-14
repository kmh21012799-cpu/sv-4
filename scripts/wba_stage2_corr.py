"""
C3a Stage 2 (focused): converse-KAM t_c (C2) vs WBA dig (C3a) at the SAME points.
Do the two 'topological' axes measure the same thing?
  corr ~ +-1  -> tautology (measuring the same thing twice)
  corr ~ 0    -> independent
  intermediate-> related but different (cf. B4 dig<->FTLE = -0.53)

Sample random core points, look up C2 t_c on its N=160 grid (nearest cell),
compute dig via RK4 (M=96, T=1000). Spearman rank correlation (robust to the
skewed distributions); reported on DETECTED points (finite t_c) and noted for
the undetected (regular) set.
"""
import os
import sys
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36        # noqa: E402
from wba import wba_dig_rk4                                          # noqa: E402
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
FIELDS = [("m4", paul_m4), ("m12", paul_m12), ("m20", paul_m20), ("m36", paul_m36)]


def lookup_tc(key, r, t):
    d = np.load(os.path.join(DATA, f"paul_{key}_tc.npz"))
    tc = d["tc"]; rho = d["rho"]; th = d["th"]
    ir = np.clip(np.searchsorted(rho, r), 0, len(rho) - 1)
    jt = np.clip(np.round(t / (2 * np.pi) * len(th)).astype(int), 0, len(th) - 1)
    return tc[ir, jt]


def main(n_ic=400, T=1000, M=96, seed=99):
    rng = np.random.default_rng(seed)
    r = rng.uniform(0.25, 0.75, n_ic)
    t = rng.uniform(0.0, 2.0 * np.pi, n_ic)
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    print(f"{n_ic} random core ICs, dig RK4 M={M} T={T}, t_c from C2 (t_f=200)")
    print(f"{'field':6s} {'n_det':>6s} {'Spearman(tc,dig)':>17s} {'undet dig med':>14s}")
    for ax, (key, ctor) in zip(axes.ravel(), FIELDS):
        f = ctor()
        dig = wba_dig_rk4(f, r, t, n_periods=T, steps_per_period=M)["dig_psi"]
        tc = lookup_tc(key, r, t)
        det = np.isfinite(tc) & (tc <= 200.0)
        und = ~det
        if det.sum() > 5:
            rho_s, p = spearmanr(tc[det], dig[det])
        else:
            rho_s = np.nan
        undmed = np.median(dig[und]) if und.sum() else np.nan
        print(f"{key:6s} {det.sum():6d} {rho_s:17.3f} {undmed:14.2f}")
        ax.scatter(tc[det], dig[det], s=8, alpha=0.5, label=f"detected (n={det.sum()})")
        if und.sum():
            ax.scatter(np.full(und.sum(), 205), dig[und], s=8, alpha=0.5, c="crimson",
                       label=f"undetected t_c=inf (n={und.sum()})")
        ax.set_title(f"{key}: Spearman(t_c,dig)={rho_s:.2f}")
        ax.set_xlabel("converse-KAM $t_c$ (C2)"); ax.set_ylabel("WBA dig (h=$\\psi$)")
        ax.legend(fontsize=7)
    fig.suptitle("converse-KAM t_c vs WBA dig at the same core points")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "wba_tc_dig_corr.png"), dpi=140)
    print("wrote wba_tc_dig_corr.png")


if __name__ == "__main__":
    main()
