"""
C0 validation figures for the Paul (2022) critical-overlap chaotic fields
(paper Sec. 7 / Fig. 6):  m = 4, 12, 20, 36.

Each is a chain of n/m resonances with amplitudes chosen so island half-width =
half the resonance spacing (critical Chirikov overlap).  Expected: visually
chaotic sea in the interior with surviving KAM surfaces near the boundaries
(perturbation ~ psi(psi-1) vanishes at psi=0,1), plus secondary island chains.

Section: zeta = 0 (mod 2pi).  Plotted in (theta, psi=rho).
"""
import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36     # noqa: E402

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")

FIELDS = {"m4": paul_m4, "m12": paul_m12, "m20": paul_m20, "m36": paul_m36}


def make_fig(key, n_seed=60, n_cross=300, n_chaos=3000, rtol=1e-8, atol=1e-10):
    f = FIELDS[key]()
    m = f.modes[0][0]
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("turbo")

    # (a) regular seeds along theta=0 -> KAM curves + island O-points
    seeds = np.linspace(0.01, 0.99, n_seed)
    for i, p0 in enumerate(seeds):
        psis, ths = f.poincare(p0, 0.0, n_cross=n_cross, rtol=rtol, atol=atol)
        ax.scatter(np.mod(ths, 2 * np.pi), psis, s=0.12,
                   c=[cmap(i / n_seed)], marker=".", linewidths=0)

    # (b) chaos-filling orbits seeded at the midpoints between adjacent
    #     resonances (psi = (n+n')/(2m)), where neighbouring separatrices meet.
    #     At critical overlap these layers connect into one sea; a long orbit
    #     fills it and reveals whether domain-spanning KAM surfaces survive.
    res_psi = sorted(n / mm for (mm, n, _) in f.modes)
    mids = [(a + b) / 2.0 for a, b in zip(res_psi[:-1], res_psi[1:])]
    # a handful of seeds suffice since the connected layers merge into one sea
    n_orbit = min(4, len(mids))
    pick = [mids[int(round(i))] for i in np.linspace(0, len(mids) - 1, n_orbit)]
    for p0 in pick:
        psis, ths = f.poincare(p0, 1e-3, n_cross=n_chaos, rtol=rtol, atol=atol)
        good = (psis > 0) & (psis < 1)
        ax.scatter(np.mod(ths[good], 2 * np.pi), psis[good], s=0.05,
                   c="k", marker=".", linewidths=0, alpha=0.5)

    for (mm, n, _) in f.modes:
        ax.axhline(n / mm, color="gray", lw=0.25, alpha=0.2)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\psi=\rho$")
    ax.set_title(f"Paul critical-overlap field  $m={m}$  (cf. paper Fig. 6)")
    fig.tight_layout()
    out = os.path.join(FIGDIR, f"paul_{key}_poincare.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("keys", nargs="*", default=["m4", "m12", "m36"])
    ap.add_argument("--nseed", type=int, default=60)
    ap.add_argument("--ncross", type=int, default=300)
    ap.add_argument("--nchaos", type=int, default=4000)
    args = ap.parse_args()
    for k in args.keys:
        make_fig(k, n_seed=args.nseed, n_cross=args.ncross, n_chaos=args.nchaos)
