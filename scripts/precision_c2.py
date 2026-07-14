"""
C2 preparation:
  P-1  grid-resolution convergence of the converse-KAM area on the KMM example 1
       (integrable, analytic answer key S_I).  Which N drives the error down, and
       how far can cell-counting realistically go?
  P-2  per-configuration cost for the four Paul fields, and a total budget.

Also a cheap sanity pass on the Paul critical-overlap fields: outermost
separatrix location (should be rho = 1/8, 7/8) and connected-chaos span.
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_kmm import example1                                   # noqa: E402
from kmm_island import island_area                               # noqa: E402
from converse_kam import detection_map, nonexistence_area        # noqa: E402
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36     # noqa: E402


def p1_grid_convergence(Ns=(90, 120, 160), t_f=200.0):
    f = example1(0.004)
    S_I = island_area(f)
    print(f"\n=== P-1  KMM ex1 area convergence  (S_I = {S_I:.5f}, t_f={t_f}) ===")
    print("   N     seconds     S(t_f)     abs.err     rel.err")
    for N in Ns:
        yg = np.linspace(-0.9, 0.9, N); zg = np.linspace(-0.9, 0.9, N)
        t0 = time.time()
        tc = detection_map(f, yg, zg, t_f=t_f, symmetry=False)
        dt = time.time() - t0
        S = nonexistence_area(tc, yg, zg, t_f)
        print(f"  {N:4d}   {dt:7.1f}   {S:.5f}   {S-S_I:+.5f}   {100*(S-S_I)/S_I:+.2f}%")


def p1_tf_convergence(N=140, t_fs=(150.0, 200.0, 300.0, 450.0)):
    """At fixed N, how much of the residual is the finite-t_f separatrix cut-off?"""
    f = example1(0.004)
    S_I = island_area(f)
    print(f"\n=== P-1b  t_f convergence at N={N}  (S_I={S_I:.5f}) ===")
    yg = np.linspace(-0.9, 0.9, N); zg = np.linspace(-0.9, 0.9, N)
    tf_max = max(t_fs)
    t0 = time.time()
    tc = detection_map(f, yg, zg, t_f=tf_max, symmetry=False)
    print(f"  (one map to t_f={tf_max} in {time.time()-t0:.1f}s; area is cumulative)")
    print("   t_f      S(t_f)     rel.err")
    for tf in t_fs:
        S = nonexistence_area(tc, yg, zg, tf)
        print(f"  {tf:5.0f}   {S:.5f}   {100*(S-S_I)/S_I:+.2f}%")


def p2_paul_cost(t_f=200.0, n_probe=24):
    """Per-point field-line cost for each Paul field; estimate converse-KAM cost.

    converse-KAM integrates state + tangent (variational) => ~2.5x the plain
    field-line RHS work, but detected points terminate early (Paul: KAM
    destroyed => most points detect), so plain-full-t_f is a CONSERVATIVE proxy.
    """
    print(f"\n=== P-2  Paul per-configuration cost  (t_f={t_f}) ===")
    n_cross = int(round(t_f / (2 * np.pi)))
    fields = [("m4", paul_m4()), ("m12", paul_m12()),
              ("m20", paul_m20()), ("m36", paul_m36())]
    var_overhead = 2.5
    print(f"  probing {n_probe} points, n_cross={n_cross}")
    print("  field   #modes   s/point(fieldline)   est s/point(ckam)")
    per_pt = {}
    for name, fld in fields:
        seeds = np.linspace(0.1, 0.9, n_probe)
        t0 = time.time()
        for p0 in seeds:
            fld.poincare(p0, 1e-3, n_cross=n_cross, rtol=1e-7, atol=1e-9)
        dt = (time.time() - t0) / n_probe
        est = dt * var_overhead
        per_pt[name] = est
        print(f"  {name:5s}   {len(fld.modes):5d}    {dt:12.4f}       {est:12.4f}")
    return per_pt


def p2_budget(per_pt, Ns=(120, 160), ncores=4):
    print(f"\n=== P-2  total budget  ({ncores} cores, 4 fields) ===")
    print("   N     points/field    est wall-time/field    est total (4 fields)")
    for N in Ns:
        pts = N * N
        for name in ["m36"]:  # worst case reported; others cheaper
            pass
        # sum over the four fields
        total = 0.0
        for name, s in per_pt.items():
            total += pts * s / ncores
        worst = pts * per_pt["m36"] / ncores
        print(f"  {N:4d}   {pts:9d}      m36 ~ {worst/60:6.1f} min      "
              f"all4 ~ {total/60:6.1f} min  (conservative)")


def sanity_paul_boundary():
    print("\n=== sanity: Paul outermost separatrix + chaos span ===")
    for name, fld in [("m4", paul_m4()), ("m12", paul_m12()),
                      ("m20", paul_m20()), ("m36", paul_m36())]:
        res = sorted(n / m for (m, n, e) in fld.modes)
        mids = [(a + b) / 2 for a, b in zip(res[:-1], res[1:])]
        p0 = mids[len(mids) // 2] if mids else res[0]
        psis, ths = fld.poincare(p0, 1e-3, n_cross=4000, rtol=1e-8, atol=1e-10)
        m = fld.modes[0][0]
        # half-width of outermost chains
        W = fld.island_half_width(m, fld.modes[0][1], fld.modes[0][2])
        print(f"  {name:5s}: {len(fld.modes)} resonances n/m in [{res[0]:.3f},{res[-1]:.3f}]"
              f"  half-width(outer)={W:.4f}  chaos span psi=[{psis.min():.3f},{psis.max():.3f}]"
              f"  (want [0.125,0.875])")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="*", default=["sanity", "p2cost", "p1"])
    args = ap.parse_args()
    if "sanity" in args.which:
        sanity_paul_boundary()
    if "p2cost" in args.which:
        pp = p2_paul_cost()
        p2_budget(pp)
    if "p1" in args.which:
        p1_grid_convergence()
    if "p1tf" in args.which:
        p1_tf_convergence()
