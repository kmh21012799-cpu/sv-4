"""
Decisive robustness test: recompute C3b's correlations using the ORIGINAL
converse-KAM / WBA (regenerated on the C3b field, saved in results/consistency_m*.npz)
and compare to the correlations using the C3b reimplementation.

If r(cKAM,VPD), r(WBA,VPD), r(cKAM,WBA) and the kappa_perp trend are the same
with the originals as with the reimplementation, the C3b conclusions (and key
question #8) are robust regardless of the diagnostic-fidelity Spearman.

The field-line diagnostics (t_c, dig) are kappa_perp-INDEPENDENT, so we reuse
them across the three kappa_perp values; only V_PD's chi is re-solved per kappa_perp.
"""
import json
import numpy as np
from scipy.stats import spearmanr

from vpd.field import paul_field
from vpd.solver3d import solve3d
from vpd.stage2 import grid_for

MS = [12, 4]
KPERPS = [1e-4, 1e-5, 1e-6]


def chi_on_grid(m, kperp, rho_grid, theta_grid):
    """Solve V_PD PDE, return chi(zeta=0) nearest-sampled onto (rho_grid,theta_grid)."""
    f = paul_field(m, chirikov=1.0)
    g = grid_for(m)
    r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                theta_period=2 * np.pi / m, core=(0.25, 0.75), method="amg", **g)
    chi = r["chi"][:, :, 0]                 # (Nrho_full, Ntheta_solve)
    rho_full = r["rho"]
    th_solve = np.arange(g["Ntheta"]) * (2 * np.pi / m) / g["Ntheta"]
    out = np.zeros((len(rho_grid), len(theta_grid)))
    for i, rr in enumerate(rho_grid):
        ir = int(np.argmin(np.abs(rho_full - rr)))
        for j, tt in enumerate(theta_grid):
            jt = int(np.argmin(np.abs(th_solve - tt)))
            out[i, j] = chi[ir, jt]
    return out.ravel()


def sp(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if np.std(a[ok]) < 1e-9 or np.std(b[ok]) < 1e-9:
        return float("nan")
    return float(spearmanr(a[ok], b[ok])[0])


def run():
    out = {}
    for m in MS:
        d = np.load(f"results/consistency_m{m}.npz")
        rho = d["rho"]; theta = d["theta"]
        tc_orig = d["tc_orig"]; tc_mine = d["tc_mine"]
        dig_orig = d["dig_orig"]; dig_mine = d["dig_mine"]
        det_orig = d["det_orig"].astype(float)
        det_mine = (d["det_mine"] > 0.5).astype(float)
        # converse-KAM vs WBA (kappa_perp independent)
        r_ckwba_orig = sp(det_orig, dig_orig)
        r_ckwba_mine = sp(det_mine, dig_mine)
        rows = {}
        for kperp in KPERPS:
            chi = chi_on_grid(m, kperp, rho, theta)
            rows[f"{kperp:.0e}"] = dict(
                r_ck_vpd_orig=sp(det_orig, chi), r_ck_vpd_mine=sp(det_mine, chi),
                # also using graded t_c (detected points), sign as -t_c so "more
                # destroyed" ~ earlier detection lines up with larger chi
                r_tc_vpd_orig=sp(-tc_orig, chi), r_tc_vpd_mine=sp(-tc_mine, chi),
                r_wba_vpd_orig=sp(dig_orig, chi), r_wba_vpd_mine=sp(dig_mine, chi),
                chi_mean=float(np.mean(chi)))
            print(f" m={m} kperp={kperp:.0e}: "
                  f"r(cKAM,VPD) orig={rows[f'{kperp:.0e}']['r_ck_vpd_orig']:.3f} "
                  f"mine={rows[f'{kperp:.0e}']['r_ck_vpd_mine']:.3f} | "
                  f"r(WBA,VPD) orig={rows[f'{kperp:.0e}']['r_wba_vpd_orig']:.3f} "
                  f"mine={rows[f'{kperp:.0e}']['r_wba_vpd_mine']:.3f}")
        out[str(m)] = dict(r_ckwba_orig=r_ckwba_orig, r_ckwba_mine=r_ckwba_mine,
                           per_kperp=rows)
        print(f" m={m}: r(cKAM,WBA) orig={r_ckwba_orig:.3f} mine={r_ckwba_mine:.3f}")
    return out


if __name__ == "__main__":
    out = run()
    json.dump(out, open("results/consistency_correlations.json", "w"), indent=2)
    print("\nwrote results/consistency_correlations.json")
