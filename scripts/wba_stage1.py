"""
C3a Stage 1 — T-convergence of the WBA dig distribution (the decisive check).

100 RANDOM ICs (not grid/line -> avoid the symmetry trap) in the core
rho in [0.25,0.75], SAME ICs across all fields and all T. dig via the vectorised
fixed-step RK4 (M=96, validated vs adaptive on the core: median matches,
corr 0.96). Question: as T grows, do the four fields' dig distributions
converge (finite-time artifact -> hypothesis A rejected) or stay separated
(real difference)?
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36        # noqa: E402
from wba import wba_dig_rk4                                          # noqa: E402
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
FIELDS = [("m4", paul_m4), ("m12", paul_m12), ("m20", paul_m20), ("m36", paul_m36)]
TS = [500, 1000, 2000, 5000]
M = 96


def main(n_ic=100, seed=20240713):
    rng = np.random.default_rng(seed)
    r = rng.uniform(0.25, 0.75, n_ic)
    th = rng.uniform(0.0, 2.0 * np.pi, n_ic)
    results = {}
    print(f"Stage 1: {n_ic} random core ICs, RK4 M={M}, h=psi")
    print(f"{'field':6s} " + " ".join(f"T={T:<5d}" for T in TS))
    for key, ctor in FIELDS:
        f = ctor()
        meds, iqrs, arrs = [], [], []
        for T in TS:
            t0 = time.time()
            dig = wba_dig_rk4(f, r, th, n_periods=T, steps_per_period=M)["dig_psi"]
            arrs.append(dig)
            meds.append(float(np.median(dig)))
            iqrs.append((float(np.percentile(dig, 25)), float(np.percentile(dig, 75))))
            print(f"    {key} T={T}: median={meds[-1]:.2f} "
                  f"IQR=[{iqrs[-1][0]:.2f},{iqrs[-1][1]:.2f}] frac<5={100*np.mean(dig<5):.0f}% "
                  f"({time.time()-t0:.0f}s)", flush=True)
        results[key] = dict(med=meds, iqr=iqrs, arr=np.array(arrs))
    np.savez(os.path.join(DATA, "wba_stage1.npz"),
             r=r, th=th, TS=TS,
             **{f"{k}_arr": v["arr"] for k, v in results.items()})

    # summary table
    print("\n=== dig median vs T (core) ===")
    print(f"{'field':6s} " + " ".join(f"{T:>7d}" for T in TS))
    for key, _ in FIELDS:
        print(f"{key:6s} " + " ".join(f"{m:7.2f}" for m in results[key]["med"]))

    # figure: median dig vs T with IQR band
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for key, _ in FIELDS:
        m = np.array(results[key]["med"])
        lo = np.array([q[0] for q in results[key]["iqr"]])
        hi = np.array([q[1] for q in results[key]["iqr"]])
        ax.plot(TS, m, "o-", label=key)
        ax.fill_between(TS, lo, hi, alpha=0.12)
    ax.set_xscale("log"); ax.set_xlabel("T (toroidal periods)")
    ax.set_ylabel("core dig (median, band=IQR)")
    ax.set_title("Stage 1: T-convergence of WBA dig, core rho in [0.25,0.75]\n"
                 "converge with T => finite-time artifact; stay apart => real")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "wba_stage1_Tconv.png"), dpi=140)
    plt.close(fig); print("wrote wba_stage1_Tconv.png")

    # histograms at largest T
    fig, ax = plt.subplots(figsize=(7, 4.6))
    bins = np.linspace(0, 15, 46)
    for key, _ in FIELDS:
        ax.hist(results[key]["arr"][-1], bins=bins, histtype="step", density=True,
                lw=1.6, label=f"{key} (med {results[key]['med'][-1]:.1f})")
    ax.set_xlabel("dig (h=psi)"); ax.set_ylabel("density")
    ax.set_title(f"core dig distribution at T={TS[-1]}")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "wba_stage1_hist.png"), dpi=140)
    plt.close(fig); print("wrote wba_stage1_hist.png")


if __name__ == "__main__":
    main()
