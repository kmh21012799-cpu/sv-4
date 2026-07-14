"""
C0 cross-check: Greene residue of the KMM example-1 (2/1) O- and X-points vs eps.

The (2,1) chain is a period-2 orbit of the Poincare map (theta advances by ~pi
per toroidal turn at the resonance).  O-point lies on the symmetry line theta=0,
X-point at theta=pi/2 (alpha=pi).  Compare the trend to paper Fig. 5:
  - R_O grows from 0 (elliptic, small perturbation) and crosses 1 when the O-point
    period-doubles;
  - R_X < 0 (hyperbolic) throughout.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1                              # noqa: E402
from kmm_island import _fixed_points_21                     # noqa: E402
from residue import greene_residue, refine_periodic         # noqa: E402

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def residue_curve(eps_list):
    RO, RX = [], []
    for eps in eps_list:
        f = example1(eps)
        psiO, psiX = _fixed_points_21(f)
        # O-point on symmetry line theta=0 (Newton-refine along psi for safety)
        pO, _ = refine_periodic(f, psiO, 0.0, q=2, fix_theta=True)
        RO.append(greene_residue(f, pO, 0.0, 2))
        # X-point at theta=pi/2
        pX, _ = refine_periodic(f, psiX, np.pi / 2.0, q=2, fix_theta=True)
        RX.append(greene_residue(f, pX, np.pi / 2.0, 2))
    return np.array(RO), np.array(RX)


def main():
    eps = np.linspace(0.0005, 0.02, 30)
    RO, RX = residue_curve(eps)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(eps, RO, "o-", ms=3, label=r"$R_O$ (O-point, elliptic)")
    ax.plot(eps, RX, "s-", ms=3, label=r"$R_X$ (X-point, hyperbolic)")
    ax.axhline(0.0, color="gray", lw=0.6)
    ax.axhline(1.0, color="gray", lw=0.6, ls="--")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("Greene residue $R$")
    ax.set_title("KMM example 1 (2/1): residue vs $\\varepsilon$  (cf. paper Fig. 5)")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(FIGDIR, "kmm_ex1_residue.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)

    # small-eps scaling check: R_O ~ eps
    print("eps      R_O        R_O/eps")
    for e, r in zip(eps[:6], RO[:6]):
        print(f"{e:.5f}  {r:.5f}   {r/e:.3f}")


if __name__ == "__main__":
    main()
