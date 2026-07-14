"""
C2 diagnostics (KMM example 1 for the grid noise; Paul m=36 for connectivity).
No C2 body, no Paul converse-KAM, no grid-size decision -- diagnosis only.

Grid noise (decision 1):
  The converse-KAM area error vs S_I is decomposed into
     (ckam area) - S_I  =  [ (ideal grid count) - S_I ]      <- pure discretisation
                         + [ (ckam area) - (ideal grid count) ] <- converse-KAM physics
  where "ideal grid count" = cell-count of the TRUE island {Psi < Psi_X} on the
  same grid (analytic, instant).  This isolates whether the +-0.5% non-monotone
  scatter is discretisation (case B) or converse-KAM timeout/physics (case A/C).
  Plus: timeout fraction and boundary fraction from the saved t_c maps.

Connectivity (decision 2):
  Chirikov overlap S and outer-separatrix location for all four Paul fields
  (code check), and long single-orbit + multi-IC span for m=36 (cantori vs
  real barrier vs code bug).
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1, KMMField                         # noqa: E402
from kmm_island import island_area, _fixed_points_21             # noqa: E402
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36     # noqa: E402
DATA = os.path.join(os.path.dirname(__file__), "..", "data")


# ---------------------------------------------------------------------------
# GRID NOISE
# ---------------------------------------------------------------------------
def ideal_grid_area(field, N, L=0.9):
    """Cell-count of the true island {Psi < Psi_X} on the (ytil,ztil) grid."""
    psiO, psiX = _fixed_points_21(field)
    PsiX = field.invariant_single(psiX, np.pi / 2.0, 0.0)
    yg = np.linspace(-L, L, N); zg = np.linspace(-L, L, N)
    Y, Z = np.meshgrid(yg, zg)
    psi, th = KMMField.from_symplectic(Y, Z, B0=field.B0)
    inside = (psi > 1e-4) & (field.invariant_single(psi, th, 0.0) < PsiX)
    dy = (2 * L) / (N - 1); dz = (2 * L) / (N - 1)
    return float(inside.sum()) * dy * dz, inside


def grid_discretisation_scan():
    f = example1(0.004)
    S_I = island_area(f)
    print(f"\n=== GRID DIAG 1: discretisation-only area (true island cell-count) ===")
    print(f"    S_I(analytic) = {S_I:.5f}")
    print("     N    ideal_grid   err_geom(%)")
    for N in [80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240]:
        A, _ = ideal_grid_area(f, N)
        print(f"   {N:4d}   {A:.5f}    {100*(A-S_I)/S_I:+.3f}")
    print("  -> if this reproduces the +-0.5% non-monotone scatter, the wobble is")
    print("     pure discretisation (case B), independent of converse-KAM.")


def timeout_boundary_from_map(npz, key="tc", t_f=200.0):
    """From a saved t_c map, split island(true) vs detected and count timeouts."""
    d = np.load(os.path.join(DATA, npz))
    tc = d[key]; yg = d["yg"]; zg = d["zg"]
    f = example1(0.004)
    psiO, psiX = _fixed_points_21(f)
    PsiX = f.invariant_single(psiX, np.pi / 2.0, 0.0)
    Y, Z = np.meshgrid(yg, zg)
    psi, th = KMMField.from_symplectic(Y, Z, B0=f.B0)
    inside = (psi > 1e-4) & (f.invariant_single(psi, th, 0.0) < PsiX)
    detected = np.isfinite(tc) & (tc <= t_f)
    N = len(yg)
    npts = inside.sum()
    timeout = inside & ~detected           # island point that never detected
    falsepos = (~inside) & detected        # outside point flagged (should be ~0)
    # boundary cells of the true island (4-neighbour edge)
    edge = np.zeros_like(inside)
    edge[1:, :] |= inside[1:, :] ^ inside[:-1, :]
    edge[:, 1:] |= inside[:, 1:] ^ inside[:, :-1]
    print(f"   {npz:16s} N={N:4d}: island cells={npts:5d}  "
          f"timeout(inside&missed)={timeout.sum():4d} ({100*timeout.sum()/npts:.2f}%)  "
          f"falsepos={falsepos.sum():3d}  boundarycells={edge.sum():4d} "
          f"({100*edge.sum()/npts:.1f}% of island)")


def grid_timeout_boundary():
    print(f"\n=== GRID DIAG 2/3: timeout + boundary from saved t_c maps ===")
    for npz, key in [("ex1_tc.npz", "tc"), ("ex1_sym.npz", "tc_g")]:
        try:
            timeout_boundary_from_map(npz, key)
        except Exception as e:
            print("   (skip", npz, ":", e, ")")


# ---------------------------------------------------------------------------
# CONNECTIVITY
# ---------------------------------------------------------------------------
def chirikov_and_separatrix():
    print(f"\n=== CONN DIAG 3/4: critical-overlap check (all four Paul fields) ===")
    for name, f in [("m4", paul_m4()), ("m12", paul_m12()),
                    ("m20", paul_m20()), ("m36", paul_m36())]:
        m = f.modes[0][0]
        dpsi = 1.0 / m
        res = sorted((n / mm, mm, n, e) for (mm, n, e) in f.modes)
        Ws = [f.island_half_width(mm, n, e) for (_, mm, n, e) in res]
        # W/(dpsi/2): should be ~1 (half-width = half spacing)
        wr = np.mean([2 * w / dpsi for w in Ws])
        # Chirikov S between neighbours
        S = np.mean([(Ws[i] + Ws[i + 1]) / dpsi for i in range(len(Ws) - 1)])
        psi_in = res[0][0]; psi_out = res[-1][0]
        sep_in = psi_in - Ws[0]; sep_out = psi_out + Ws[-1]
        print(f"  {name:4s}: dpsi=1/{m}={dpsi:.4f}  mean W/(dpsi/2)={wr:.3f}  "
              f"Chirikov S={S:.3f}  inner sep={sep_in:.4f} outer sep={sep_out:.4f} "
              f"(want 0.125, 0.875)")


def m36_long_orbit(turns=(1000, 10000, 100000)):
    print(f"\n=== CONN DIAG 1: m=36 single orbit, cumulative span vs #turns ===")
    f = paul_m36()
    res = sorted(n / m for (m, n, e) in f.modes)
    mids = [(a + b) / 2 for a, b in zip(res[:-1], res[1:])]
    p0 = mids[len(mids) // 2]
    t0 = time.time()
    psis, ths = f.poincare(p0, 1e-3, n_cross=max(turns), rtol=1e-8, atol=1e-10)
    for T in turns:
        seg = psis[:T + 1]
        print(f"   after {T:7d} turns: psi span [{seg.min():.3f}, {seg.max():.3f}] "
              f"width={seg.max()-seg.min():.3f}")
    print(f"   (seed psi0={p0:.3f}, {time.time()-t0:.0f}s)  want [0.125,0.875] if connected")


def m36_multi_ic(seeds=(0.20, 0.35, 0.50, 0.65, 0.80), turns=20000):
    print(f"\n=== CONN DIAG 2: m=36 multi-IC spans ({turns} turns each) ===")
    f = paul_m36()
    spans = []
    for p0 in seeds:
        psis, _ = f.poincare(p0, 1e-3, n_cross=turns, rtol=1e-8, atol=1e-10)
        spans.append((p0, psis.min(), psis.max()))
        print(f"   seed {p0:.2f} -> [{psis.min():.3f}, {psis.max():.3f}]")
    # do the spans overlap into one connected range?
    lo = min(s[1] for s in spans); hi = max(s[2] for s in spans)
    print(f"   union span = [{lo:.3f}, {hi:.3f}]  (overlap => one sea; gaps => barrier)")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="*",
                    default=["griddisc", "gridto", "chirikov", "multi", "long"])
    args = ap.parse_args()
    if "griddisc" in args.which:
        grid_discretisation_scan()
    if "gridto" in args.which:
        grid_timeout_boundary()
    if "chirikov" in args.which:
        chirikov_and_separatrix()
    if "multi" in args.which:
        m36_multi_ic()
    if "long" in args.which:
        m36_long_orbit()
