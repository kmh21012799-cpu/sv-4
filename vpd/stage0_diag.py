"""
Stage 0 -- diagnose why the C3b validation passed on the wrong (envelope-less)
field, and show the boundary parallel-flux contamination.

Compares the angle-averaged local V_PD indicator chi(rho) over the FULL domain
[0,1] for the OLD constant-amplitude field vs the NEW Paul-envelope field.
Envelope => perturbation vanishes at rho=0,1 => chi -> 0 there (no parallel
flux across the Dirichlet boundaries).  Constant amplitude => chi stays large
near the boundaries (spurious parallel flux).
"""
import numpy as np
from vpd.solver3d import solve3d
from vpd.field import paul_field as old_field       # OLD constant-amplitude
from vpd.paulsolve import paul_critical             # NEW Paul envelope


def grid_for_field(field, m):
    n_max = max(n for (mm, n, e) in field.modes)
    Nzeta = int(2 ** np.ceil(np.log2(max(32, 2.6 * n_max))))
    return dict(Nrho=129, Ntheta=32, Nzeta=Nzeta)


def chi_profile(field, m, kperp=1e-6):
    g = grid_for_field(field, m)
    r = solve3d(field, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                theta_period=2 * np.pi / m, core=(0.25, 0.75), method="amg", **g)
    chi = r["chi"]                    # (Nrho, Ntheta, Nzeta)
    rho = r["rho"]
    prof = chi.mean(axis=(1, 2))      # angle-averaged chi(rho)
    return rho, prof, r


if __name__ == "__main__":
    m = 12
    print(f"Stage 0 diagnostic (m={m}, kperp=1e-6): angle-averaged chi(rho)\n")
    ro, po, _ = chi_profile(old_field(m, chirikov=1.0), m)
    rn, pn, _ = chi_profile(paul_critical(m), m)
    print(f"{'rho':>6} {'OLD chi':>10} {'NEW chi':>10}")
    for rr in [0.02, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 0.85, 0.90, 0.95, 0.98]:
        io = int(np.argmin(np.abs(ro - rr))); ino = int(np.argmin(np.abs(rn - rr)))
        print(f"{rr:6.2f} {po[io]:10.3f} {pn[ino]:10.3f}")
    # boundary-band means
    def band(rho, prof, lo, hi):
        m_ = (rho >= lo) & (rho <= hi)
        return float(prof[m_].mean())
    print("\nnear-boundary chi (mean over rho in [0,0.1] and [0.9,1]):")
    print(f"  OLD: lower={band(ro,po,0,0.1):.3f}  upper={band(ro,po,0.9,1):.3f}")
    print(f"  NEW: lower={band(rn,pn,0,0.1):.3f}  upper={band(rn,pn,0.9,1):.3f}")
    print("  (NEW should be ~0 at the boundaries; OLD should not.)")
    np.savez("results/stage0_chiprofile.npz", ro=ro, po=po, rn=rn, pn=pn)
