"""C2 summary figure: converse-KAM metrics vs m, with Paul's transport trend."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from paul_ckam_analysis import stats                             # noqa: E402
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")

ms = [4, 12, 20, 36]
keys = ["m4", "m12", "m20", "m36"]
S = [stats(k) for k in keys]
area = [s["frac_full"] for s in S]
core = [s["frac_sub"] for s in S]
undet = [s["undet"] for s in S]
med = [s["med"] for s in S]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

# left: converse-KAM metrics vs m
err = 1.0   # ~+-1% discretisation error bar (single-resonance gate / KMM scan)
ax1.errorbar(ms, area, yerr=err, fmt="o-", capsize=3, label="non-existence area (% domain)")
ax1.plot(ms, core, "s--", label=r"core $\rho\in[0.25,0.75]$ detected (%)")
ax1.plot(ms, [100 - u for u in undet], "^:", label="detected (% domain)")
ax1.set_xlabel("m (resonance count 3/9/15/27)")
ax1.set_ylabel("percent")
ax1.set_title("converse-KAM: core is 100% for all;\nfull-domain area saturates by m=12")
ax1.set_ylim(50, 105)
ax1.legend(fontsize=8)

# right: the direction relative to Paul's transport
ax2.plot(ms, area, "o-", label="converse-KAM area (this work)")
ax2.annotate("Paul V_PD / heat transport\nDECREASES with m\n(m=4 leaks most, m=36 insulates)",
             xy=(36, 74), xytext=(12, 60), fontsize=9,
             arrowprops=dict(arrowstyle="->", color="crimson"), color="crimson")
ax2.annotate("", xy=(36, 57), xytext=(4, 68),
             arrowprops=dict(arrowstyle="->", color="crimson", lw=2))
ax2.set_xlabel("m")
ax2.set_ylabel("converse-KAM area (% domain)")
ax2.set_title("converse-KAM area is ANTI-correlated with transport\n"
              "(and t_c medians ~19-21 are flat) -> does not track V_PD")
ax2.set_ylim(50, 80)
ax2.legend(fontsize=8, loc="lower right")

fig.tight_layout()
fig.savefig(os.path.join(FIG, "paul_ckam_summary.png"), dpi=140)
print("wrote paul_ckam_summary.png")
print("area %dom:", [round(a, 1) for a in area])
print("core %   :", [round(c, 1) for c in core])
print("t_c med  :", [round(m, 1) for m in med])
