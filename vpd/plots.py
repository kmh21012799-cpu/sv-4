"""Figure generation for the C3b record."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fig_validation1(path="results/validation1.json", out="figures/fig1_island.png"):
    with open(path) as fh:
        data = json.load(fh)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, d in data.items():
        rows = d["rows"]
        eps = np.array([r["eps"] for r in rows])
        vpd = np.array([r["V_PD"] for r in rows])
        ec = d["eps_crit"]
        ax1.plot(eps / ec, vpd, "o-", label=f"$\\kappa_\\perp$={key}")
    ax1.axvline(1.0, color="k", ls="--", lw=1, label="$\\epsilon_{crit}=\\sqrt{\\kappa_\\perp}/2$")
    ax1.set_xscale("log"); ax1.set_xlabel("$\\epsilon_{21}/\\epsilon_{crit}$")
    ax1.set_ylabel("$V_{PD}$"); ax1.set_title("(a) transition at $\\epsilon_{crit}$")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    # panel b: sqrt(eps) scaling for smallest kperp
    key = sorted(data.keys())[-1]
    rows = data[key]["rows"]
    eps = np.array([r["eps"] for r in rows]); vpd = np.array([r["V_PD"] for r in rows])
    m = vpd > 0.02
    ax2.loglog(eps[m], vpd[m], "o-", label=f"$V_{{PD}}$ ($\\kappa_\\perp$={key})")
    ref = eps[m]
    ax2.loglog(ref, vpd[m][0] * np.sqrt(ref / ref[0]), "k--", label="$\\propto\\sqrt{\\epsilon}$")
    ax2.set_xlabel("$\\epsilon_{21}$"); ax2.set_ylabel("$V_{PD}$")
    ax2.set_title("(b) large-$\\epsilon$ scaling"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print("wrote", out)


def fig_validation2(path="results/validation2.json", out="figures/fig2_chaos.png"):
    with open(path) as fh:
        data = json.load(fh)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for key, d in sorted(data.items()):
        rows = d["rows"]
        kp = np.array([r["kperp"] for r in rows])
        vp = np.array([r["V_PD"] for r in rows])
        ax.semilogx(kp, vp, "o-", label=f"s={d['chirikov']:.2f}")
        ax.axvline(d["D_QL"], color=ax.lines[-1].get_color(), ls=":", lw=1, alpha=0.6)
    ax.axhline(0.5, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_xlabel("$\\kappa_\\perp$"); ax.set_ylabel("$V_{PD}$")
    ax.set_title("Chaotic layer: $V_{PD}$ vs $\\kappa_\\perp$\n(dotted = quasilinear $D_{QL}$)")
    ax.invert_xaxis(); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print("wrote", out)


def fig_stage2(path="results/stage2.json", out="figures/fig3_fourfields.png"):
    with open(path) as fh:
        data = json.load(fh)
    ms = sorted(int(k) for k in data)
    kperps = [r["kperp"] for r in data[str(ms[0])]["rows"]]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for m in ms:
        rows = data[str(m)]["rows"]
        kp = [r["kperp"] for r in rows]
        ax1.semilogx(kp, [r["V_PD"] for r in rows], "o-", label=f"m={m}")
        ax2.semilogx(kp, [r["DeltaT"] for r in rows], "s-", label=f"m={m}")
    for ax in (ax1, ax2):
        ax.invert_xaxis(); ax.set_xlabel("$\\kappa_\\perp$"); ax.legend(); ax.grid(alpha=0.3)
    ax1.set_ylabel("$V_{PD}$"); ax1.set_title("(a) $V_{PD}$: largest for m=4")
    ax2.set_ylabel("$\\Delta T$"); ax2.set_title("(b) $\\Delta T$: largest for m=36 (insulates)")
    fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)


def fig_stage3(path="results/stage3.json", out="figures/fig4_correlation.png"):
    with open(path) as fh:
        data = json.load(fh)
    ms = sorted(int(k) for k in data)
    kperps = sorted({kp for m in ms for kp in
                     [v["kperp"] for v in data[str(m)]["per_kperp"].values()]},
                    reverse=True)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for m in ms:
        pk = data[str(m)]["per_kperp"]
        kps = sorted(pk.values(), key=lambda d: -d["kperp"])
        ax.semilogx([d["kperp"] for d in kps], [d["r_ck_vpd"] for d in kps],
                    "o-", label=f"m={m}")
    ax.axhline(0, color="k", lw=0.6)
    ax.invert_xaxis(); ax.set_xlabel("$\\kappa_\\perp$")
    ax.set_ylabel("Spearman $r$(converse-KAM, $V_{PD}$)")
    ax.set_title("Does converse-KAM $\\leftrightarrow$ $V_{PD}$ correlation\ngrow as $\\kappa_\\perp\\to0$?")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "1"):
        fig_validation1()
    if which in ("all", "2"):
        fig_validation2()
    if which in ("all", "3"):
        fig_stage2()
    if which in ("all", "4"):
        fig_stage3()
