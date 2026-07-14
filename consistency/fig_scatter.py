"""Scatter: original vs reimplemented t_c and dig, on the same points."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, figsize=(10, 9))
for col, m in enumerate([12, 4]):
    d = np.load(f"results/consistency_m{m}.npz")
    tco, tcm = d["tc_orig"], d["tc_mine"]
    dgo, dgm = d["dig_orig"], d["dig_mine"]
    both = np.isfinite(tco) & (tco <= 200) & np.isfinite(tcm) & (tcm <= 200)
    ax = axs[0, col]
    ax.scatter(tco[both], tcm[both], s=8, alpha=0.5)
    lim = [0, max(tco[both].max(), tcm[both].max()) * 1.05]
    ax.plot(lim, lim, "k--", lw=1, label="y=x")
    ax.set_xlabel("original $t_c$ (zeta)"); ax.set_ylabel("C3b $t_c$ (zeta)")
    ax.set_title(f"converse-KAM $t_c$  (m={m})\nSpearman "
                 f"{__import__('scipy.stats', fromlist=['spearmanr']).spearmanr(tco[both], tcm[both])[0]:.3f}")
    ax.legend(); ax.grid(alpha=0.3)
    ax = axs[1, col]
    fin = np.isfinite(dgo) & np.isfinite(dgm)
    ax.scatter(dgo[fin], dgm[fin], s=8, alpha=0.5, color="C1")
    ax.axvline(5, color="r", ls=":", lw=1); ax.axhline(5, color="r", ls=":", lw=1)
    ax.set_xlabel("original dig"); ax.set_ylabel("C3b dig")
    ax.set_title(f"WBA dig  (m={m})\nSpearman "
                 f"{__import__('scipy.stats', fromlist=['spearmanr']).spearmanr(dgo[fin], dgm[fin])[0]:.3f} "
                 f"(dashed=chaos threshold)")
    ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("figures/figA_consistency.png", dpi=120)
print("wrote figures/figA_consistency.png")
