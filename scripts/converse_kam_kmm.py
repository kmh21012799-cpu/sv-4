"""
C1 validation of converse-KAM 3D on the KMM fields.

Validations (answer keys):
  (1) area convergence : example 1 (integrable) non-existence area S(t_f) rises
      monotonically to the analytic island area S_I  (paper Fig. 4).
  (2) first-detection time : min t_c over the island ~ (pi/2) T / sqrt(R),
      T = 4 pi (period-2 island-centre period), R = Greene residue  (Fig. 5).
  (3) symmetry line vs 2-D grid : Thm 3.2 (symmetry line) detects at ~half the
      time of Thm 3.1 but over-samples islands and cannot estimate area.
  (4) examples 2, 3 : both island chains are detected as non-existence regions.

Grid: symplectic (ytil, ztil); area in these coords = toroidal flux (eq. 12).
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1, example2, example3, KMMField           # noqa: E402
from kmm_island import island_area, _fixed_points_21                    # noqa: E402
from residue import greene_residue                                      # noqa: E402
from converse_kam import detect_tc, detection_map, nonexistence_area    # noqa: E402

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATADIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATADIR, exist_ok=True)


def compute_map(field, N=160, L=0.9, t_f=200.0, symmetry=False, tag=""):
    yg = np.linspace(-L, L, N)
    zg = np.linspace(-L, L, N)
    t0 = time.time()
    tc = detection_map(field, yg, zg, t_f=t_f, symmetry=symmetry)
    print(f"[{tag}] map {N}x{N} t_f={t_f} sym={symmetry} in {time.time()-t0:.1f}s")
    return tc, yg, zg


# ---------------------------------------------------------------------------
def validate_example1(N=160, t_f=200.0):
    f = example1(0.004)
    S_I = island_area(f)
    tc, yg, zg = compute_map(f, N=N, L=0.9, t_f=t_f, symmetry=False, tag="ex1")
    np.savez(os.path.join(DATADIR, "ex1_tc.npz"), tc=tc, yg=yg, zg=zg, S_I=S_I)

    # (1) area convergence
    tfs = np.linspace(20, t_f, 60)
    S = [nonexistence_area(tc, yg, zg, tf) for tf in tfs]
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(tfs, S, "-o", ms=3, label=r"$S(t_f)$ converse-KAM non-existence area")
    ax.axhline(S_I, color="r", ls="--", label=fr"$S_I={S_I:.4f}$ (analytic island area)")
    ax.set_xlabel(r"timeout $t_f$ (in $\varphi$)")
    ax.set_ylabel("non-existence area (toroidal flux)")
    ax.set_title("KMM ex 1: converse-KAM area $\\to$ island area (cf. Fig. 4)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "ck_ex1_area_convergence.png"), dpi=140)
    plt.close(fig)
    print(f"[ex1] S(t_f={t_f})={S[-1]:.4f}  S_I={S_I:.4f}  ratio={S[-1]/S_I:.3f}")

    # detection map figure
    plot_map(tc, yg, zg, "KMM ex 1 (2/1): converse-KAM detection time $t_c$",
             "ck_ex1_map.png", overlay_field=f)


def plot_map(tc, yg, zg, title, fname, overlay_field=None):
    fig, ax = plt.subplots(figsize=(6, 5.4))
    masked = np.where(np.isfinite(tc), tc, np.nan)
    im = ax.imshow(masked, origin="lower", extent=[yg[0], yg[-1], zg[0], zg[-1]],
                   aspect="equal", cmap="viridis")
    fig.colorbar(im, ax=ax, label=r"$t_c$ (detected non-existence)")
    ax.set_xlabel(r"$\tilde y$"); ax.set_ylabel(r"$\tilde z$")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, fname), dpi=140)
    plt.close(fig)
    print("wrote", fname)


# ---------------------------------------------------------------------------
def validate_tc_residue():
    """min t_c near the O-point vs eps, compared to (pi/2) T / sqrt(R), T=4pi."""
    eps_list = [0.002, 0.003, 0.004, 0.006, 0.008, 0.01]
    T = 4.0 * np.pi
    rows = []
    for eps in eps_list:
        f = example1(eps)
        psiO, _ = _fixed_points_21(f)
        R = greene_residue(f, psiO, 0.0, 2)
        # scan a small neighbourhood of the O-line for the minimum t_c
        best = np.inf
        for dp in np.linspace(-0.02, 0.02, 21):
            for th in [0.0, 0.3, 0.6]:
                tc = detect_tc(f, psiO + dp, th, t_f=300.0)
                if tc < best:
                    best = tc
        pred = (np.pi / 2.0) * T / np.sqrt(R)
        rows.append((eps, R, best, pred))
        print(f"eps={eps:.4f} R_O={R:.4f}  t_c(min)={best:.2f}  (pi/2)T/sqrt(R)={pred:.2f}")

    eps_a = np.array([r[0] for r in rows])
    tcmin = np.array([r[2] for r in rows])
    pred = np.array([r[3] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(eps_a, tcmin, "o-", ms=4, label=r"measured min $t_c$")
    ax.plot(eps_a, pred, "s--", ms=4, label=r"$(\pi/2)\,T/\sqrt{R}$,  $T=4\pi$")
    ax.set_xlabel(r"$\varepsilon$"); ax.set_ylabel(r"first-detection time $t_c$")
    ax.set_title("KMM ex 1: converse-KAM $t_c$ vs residue prediction (cf. Fig. 5)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "ck_ex1_tc_residue.png"), dpi=140)
    plt.close(fig)
    print("wrote ck_ex1_tc_residue.png")


# ---------------------------------------------------------------------------
def validate_symmetry_vs_grid(N=120, t_f=200.0):
    """Two things the paper warns about (and the D1/B2 trap):

      (a) the theta=0 symmetry line passes through EVERY primary island's
          elliptic point -> it over-samples islands;
      (b) the Thm 3.2 short-cut (condition i only) is valid ONLY on the symmetry
          line; applied to a 2-D grid it OVER-detects, so area must come from the
          full Thm 3.1 grid.

    We show (b) directly by mapping both criteria on the same grid, and quantify
    the paper's 'factor of 2' on the symmetry line via the lambda-condition
    timing.
    """
    from converse_kam import _aug_rhs, _lambda
    from scipy.integrate import solve_ivp
    f = example1(0.004)
    tc_g, yg, zg = compute_map(f, N=N, L=0.9, t_f=t_f, symmetry=False, tag="ex1-grid(3.1)")
    tc_s, _, _ = compute_map(f, N=N, L=0.9, t_f=t_f, symmetry=True, tag="ex1-shortcut(3.2)")
    np.savez(os.path.join(DATADIR, "ex1_sym.npz"), tc_g=tc_g, tc_s=tc_s,
             yg=yg, zg=zg)
    Sg = nonexistence_area(tc_g, yg, zg, t_f)
    Ss = nonexistence_area(tc_s, yg, zg, t_f)
    print(f"area: Thm3.1 grid = {Sg:.4f}   Thm3.2 short-cut on 2-D grid = {Ss:.4f}"
          f"   (difference {100*(Ss-Sg)/Sg:+.1f}%)")

    # (b') on the theta=0 line: first beta sign-change (Thm 3.2 detection) vs
    # first lambda<0 (the extra Thm 3.1 condition) -> ratio ~ 1/2.
    print("theta=0 line:   psi   t_c(beta sign change)   t(first lambda<0)   ratio")
    for p0 in [0.10, 0.13, 0.16]:
        sol = solve_ivp(_aug_rhs(f), (0, 80), [p0, 0.0, 1.0, 0.0],
                        rtol=1e-9, atol=1e-11, dense_output=True, max_step=0.4)
        phi = np.linspace(0, 80, 4000); Y = sol.sol(phi)
        et = Y[3]
        lam = np.array([_lambda(f, Y[0, i], Y[1, i], phi[i], [Y[2, i], Y[3, i]])
                        for i in range(len(phi))])
        nz = np.nonzero(np.abs(et) > 1e-6)[0][0]
        cr = [phi[k] for k in range(nz + 1, len(et))
              if np.sign(et[k]) != np.sign(et[nz])]
        tbeta = cr[0] if cr else np.nan
        tlam = phi[np.nonzero(lam < 0)[0][0]] if np.any(lam < 0) else np.nan
        print(f"              {p0:.3f}      {tbeta:8.2f}            {tlam:8.2f}"
              f"          {tlam/tbeta:.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, tc, ttl in ((axes[0], tc_g, "Thm 3.1  (cond. i + ii)"),
                        (axes[1], tc_s, "Thm 3.2 short-cut  (cond. i only)")):
        m = np.where(np.isfinite(tc), tc, np.nan)
        im = ax.imshow(m, origin="lower", extent=[yg[0], yg[-1], zg[0], zg[-1]],
                       aspect="equal", cmap="viridis", vmin=0, vmax=t_f)
        ax.set_title(ttl); ax.set_xlabel(r"$\tilde y$"); ax.set_ylabel(r"$\tilde z$")
        fig.colorbar(im, ax=ax, fraction=0.046, label=r"$t_c$")
    fig.suptitle(f"KMM ex 1 (integrable): short-cut and full criterion coincide  "
                 f"(area {Ss:.3f} vs {Sg:.3f}); factor-of-2 is in the $\\lambda$-timing")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "ck_ex1_sym_vs_grid.png"), dpi=140)
    plt.close(fig)
    print("wrote ck_ex1_sym_vs_grid.png")


# ---------------------------------------------------------------------------
def validate_examples_23(N=140, t_f=200.0):
    for name, f in (("ex2", example2(0.003)), ("ex3", example3(0.001, 0.01))):
        tc, yg, zg = compute_map(f, N=N, L=0.95, t_f=t_f, symmetry=False, tag=name)
        np.savez(os.path.join(DATADIR, f"{name}_tc.npz"), tc=tc, yg=yg, zg=zg)
        plot_map(tc, yg, zg, f"KMM {name}: converse-KAM detection $t_c$",
                 f"ck_{name}_map.png")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="*",
                    default=["ex1", "tcres", "sym", "ex23"])
    ap.add_argument("--N", type=int, default=160)
    ap.add_argument("--tf", type=float, default=200.0)
    args = ap.parse_args()
    if "ex1" in args.which:
        validate_example1(N=args.N, t_f=args.tf)
    if "tcres" in args.which:
        validate_tc_residue()
    if "sym" in args.which:
        validate_symmetry_vs_grid(N=min(args.N, 120), t_f=args.tf)
    if "ex23" in args.which:
        validate_examples_23(N=min(args.N, 140), t_f=args.tf)
