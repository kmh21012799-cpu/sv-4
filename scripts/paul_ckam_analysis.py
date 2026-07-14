"""
C2-C analysis: read the four Paul converse-KAM t_c maps and produce
  - comparison table (area full + rho[0.25,0.75], undetected%, t_c median/IQR)
  - overlaid t_c histograms
  - 2x2 detection maps (undetected shown white)
Primary signal is the t_c distribution; area is carried with a discretisation
error bar (measured on the single-resonance gate ~1-1.5%).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from converse_kam import nonexistence_area_uv                       # noqa: E402
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
KEYS = ["m4", "m12", "m20", "m36"]
LABEL = {"m4": "m=4", "m12": "m=12", "m20": "m=20", "m36": "m=36"}


def load(key, tag="tc"):
    d = np.load(os.path.join(DATA, f"paul_{key}_{tag}.npz"))
    return d["tc"], d["rho"], d["th"]


def stats(key, tag="tc", t_f=200.0):
    tc, rho, th = load(key, tag)
    valid = ~np.isnan(tc)
    det = valid & np.isfinite(tc) & (tc <= t_f)
    full_dom = 1.0 * 2.0 * np.pi
    band = (rho >= 0.25) & (rho <= 0.75)
    S_full = nonexistence_area_uv(np.where(det, tc, np.inf), rho, th, t_f)
    subfrac = 100.0 * det[band, :].sum() / valid[band, :].sum()
    tcv = tc[det]
    med = float(np.median(tcv)); q1, q3 = np.percentile(tcv, [25, 75])
    undet = 100 * (valid & ~det).sum() / valid.sum()
    return dict(key=key, S_full=S_full, frac_full=100 * S_full / full_dom,
                frac_sub=subfrac, undet=undet, med=med, q1=q1, q3=q3, tcv=tcv)


def table(tag="tc"):
    print(f"\n{'field':6s} {'area_full':>9s} {'%dom':>6s} {'det[.25,.75]%':>13s} "
          f"{'undet%':>7s} {'t_c med':>8s} {'IQR':>14s}")
    rows = []
    for k in KEYS:
        s = stats(k, tag)
        rows.append(s)
        print(f"{LABEL[k]:6s} {s['S_full']:9.3f} {s['frac_full']:6.1f} "
              f"{s['frac_sub']:13.1f} {s['undet']:7.1f} "
              f"{s['med']:8.1f} [{s['q1']:5.1f},{s['q3']:5.1f}]")
    return rows


def fig_histograms(tag="tc", t_f=200.0):
    fig, ax = plt.subplots(figsize=(7, 4.6))
    bins = np.linspace(0, 120, 61)
    for k in KEYS:
        s = stats(k, tag, t_f)
        ax.hist(s["tcv"], bins=bins, histtype="step", density=True, lw=1.6,
                label=f"{LABEL[k]} (med {s['med']:.0f})")
    ax.set_xlabel(r"detection time $t_c$")
    ax.set_ylabel("density (detected points)")
    ax.set_title("converse-KAM $t_c$ distribution, four Paul fields (t_f=200)")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "paul_tc_hist.png"), dpi=140)
    plt.close(fig); print("wrote paul_tc_hist.png")


def fig_maps(tag="tc"):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, k in zip(axes.ravel(), KEYS):
        tc, rho, th = load(k, tag)
        m = np.where(np.isfinite(tc), tc, np.nan)
        im = ax.imshow(m, origin="lower", extent=[0, 2 * np.pi, 0, 1],
                       aspect="auto", cmap="viridis", vmin=0, vmax=80)
        ax.axhline(0.125, color="w", lw=0.5, ls=":")
        ax.axhline(0.875, color="w", lw=0.5, ls=":")
        ax.set_title(LABEL[k]); ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\rho$")
        fig.colorbar(im, ax=ax, fraction=0.046, label=r"$t_c$")
    fig.suptitle("converse-KAM detection maps (white = undetected)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "paul_ckam_maps.png"), dpi=130)
    plt.close(fig); print("wrote paul_ckam_maps.png")


def convergence_table():
    print("\n=== t_f convergence (N=60, tag conv60): area_full at t_f=200/400/800 ===")
    for k in KEYS:
        tc, rho, th = load(k, "conv60")
        vals = []
        for tf in [200.0, 400.0, 800.0]:
            det = (~np.isnan(tc)) & np.isfinite(tc) & (tc <= tf)
            vals.append(nonexistence_area_uv(np.where(det, tc, np.inf), rho, th, tf))
        undet = 100 * ((~np.isnan(tc)) & ~((np.isfinite(tc)) & (tc <= 800))).sum() / (~np.isnan(tc)).sum()
        print(f"  {LABEL[k]:5s}: {vals[0]:.3f} / {vals[1]:.3f} / {vals[2]:.3f}"
              f"  (undetected@800 = {undet:.1f}%)  saturated={abs(vals[2]-vals[0])<1e-6}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="*", default=["table", "hist", "maps", "conv"])
    ap.add_argument("--tag", default="tc")
    args = ap.parse_args()
    if "conv" in args.which:
        convergence_table()
    if "table" in args.which:
        table(args.tag)
    if "hist" in args.which:
        fig_histograms(args.tag)
    if "maps" in args.which:
        fig_maps(args.tag)
