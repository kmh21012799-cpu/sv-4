"""
C0 validation figures for the KMM field.

  fig 1 : example 1 (integrable, single 2/1) Poincare  +  invariant Psi contours
          -> the two must coincide  (integrability check, paper Fig. 2-ish)
  fig 2 : example 2 (2/1 + 3/2) Poincare
  fig 3 : example 3 (2/1 + 5/4) Poincare

Section: phi = 0 (mod 2 pi).  Plotted in symplectic coords (ytil, ztil) where
area = toroidal flux.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1, example2, example3, KMMField  # noqa: E402

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIGDIR, exist_ok=True)


def poincare_scatter(ax, field, seeds, n_cross=300, s=0.25, color="k"):
    for (p0, t0) in seeds:
        psis, ths = field.poincare(p0, t0, n_cross=n_cross)
        y, z = field.to_symplectic(psis, ths)
        ax.scatter(y, z, s=s, c=color, marker=".", linewidths=0)


def radial_seeds(psi_lo, psi_hi, n, theta0=0.0):
    return [(p, theta0) for p in np.linspace(psi_lo, psi_hi, n)]


def fig_example1():
    eps = 0.004
    f = example1(eps)
    fig, ax = plt.subplots(figsize=(6, 6))

    # Poincare
    seeds = radial_seeds(0.02, 0.42, 26, theta0=0.05)
    poincare_scatter(ax, f, seeds, n_cross=250, s=0.4)

    # invariant contours (should coincide with the section, since integrable)
    g = np.linspace(-0.95, 0.95, 600)
    Y, Z = np.meshgrid(g, g)
    psi, th = KMMField.from_symplectic(Y, Z, B0=f.B0)
    Psi = f.invariant_single(psi, th, 0.0)
    ax.contour(Y, Z, Psi, levels=35, colors="tab:red", linewidths=0.4, alpha=0.7)

    ax.set_aspect("equal")
    ax.set_xlabel(r"$\tilde y=\sqrt{2\psi}\cos\vartheta$")
    ax.set_ylabel(r"$\tilde z=\sqrt{2\psi}\sin\vartheta$")
    ax.set_title(f"KMM example 1 (integrable, 2/1), $\\varepsilon={eps}$\n"
                 "black = Poincare, red = invariant $\\Psi$ contours")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "kmm_ex1_poincare_invariant.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)


def fig_example2():
    eps = 0.003
    f = example2(eps)
    fig, ax = plt.subplots(figsize=(6, 6))
    seeds = radial_seeds(0.02, 0.45, 40, theta0=0.05)
    poincare_scatter(ax, f, seeds, n_cross=400, s=0.3)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\tilde y$"); ax.set_ylabel(r"$\tilde z$")
    ax.set_title(f"KMM example 2 (2/1 + 3/2), $\\varepsilon={eps}$")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "kmm_ex2_poincare.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)


def fig_example3():
    eps21, eps54 = 0.001, 0.01
    f = example3(eps21, eps54)
    fig, ax = plt.subplots(figsize=(6, 6))
    seeds = radial_seeds(0.02, 0.5, 45, theta0=0.05)
    poincare_scatter(ax, f, seeds, n_cross=400, s=0.3)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\tilde y$"); ax.set_ylabel(r"$\tilde z$")
    ax.set_title(f"KMM example 3 (2/1 + 5/4), "
                 f"$\\varepsilon_{{21}}={eps21},\\ \\varepsilon_{{54}}={eps54}$")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "kmm_ex3_poincare.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    fig_example1()
    fig_example2()
    fig_example3()
