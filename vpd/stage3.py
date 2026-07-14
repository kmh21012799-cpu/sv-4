"""
Stage 3 -- the three-axis comparison.

On the SAME field, SAME (rho, theta) grid (zeta=0 slice), SAME core
[0.25, 0.75], compute all three indicators point by point:

  [topology] converse-KAM : detected(rho,theta), t_c
  [dynamics] WBA          : dig(rho,theta)
  [transport] V_PD        : local indicator chi(rho,theta) = Theta(kpar|grad_par T|^2
                                                              - kperp|grad_perp T|^2)

and correlate the three pairs (Spearman).  Central question: does the
converse-KAM <-> V_PD correlation grow as kperp -> 0 (Paul's expectation)?
"""
import json
import numpy as np
from scipy.stats import spearmanr

from vpd.field import paul_field, core_resonances
from vpd.solver3d import solve3d
from vpd.diagnostics import run_diagnostics
from vpd.stage2 import grid_for

MS = [4, 12, 20, 36]
KPERPS = [1e-4, 1e-5, 1e-6]
CORE = (0.25, 0.75)


def chi_zeta0_on_core(solve_result, core):
    """Local V_PD indicator on the zeta=0 slice, restricted to core rho."""
    chi = solve_result["chi"]           # (Nrho, Ntheta, Nzeta)
    rho = solve_result["rho"]
    cmask = (rho >= core[0]) & (rho <= core[1])
    return chi[cmask, :, 0], rho[cmask], cmask


def safe_spearman(a, b):
    a = np.asarray(a).ravel().astype(float)
    b = np.asarray(b).ravel().astype(float)
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan"), float("nan")
    r, p = spearmanr(a, b)
    return float(r), float(p)


def run():
    out = {}
    for m in MS:
        f = paul_field(m, chirikov=1.0)
        g = grid_for(m)
        # first solve at coarsest kperp to get the core rho grid & theta grid
        # (theta grid is fixed by g)
        theta = np.arange(g["Ntheta"]) * (2 * np.pi / m) / g["Ntheta"]
        rho_full = np.linspace(0.0, 1.0, g["Nrho"])
        cmask = (rho_full >= CORE[0]) & (rho_full <= CORE[1])
        rho_core = rho_full[cmask]
        # subsample rho for the (expensive) field-line diagnostics
        sub = max(1, len(rho_core) // 48)
        rho_diag = rho_core[::sub]

        print(f"[m={m}] diagnostics on {len(rho_diag)}x{len(theta)} grid ...")
        diag = run_diagnostics(f, rho_diag, theta, n_periods=1200,
                               substeps=28, ck_periods=150)

        per_kperp = {}
        x0 = None
        for kperp in KPERPS:
            r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / m, core=CORE, method="amg",
                        x0=x0, **g)
            x0 = r["T"]
            chi_core, rc, _ = chi_zeta0_on_core(r, CORE)
            chi_sub = chi_core[::sub, :]        # align to rho_diag
            # correlations (flatten over the core rho x theta grid)
            det = diag["detected"].astype(float)
            dig = diag["dig"]
            r_ck_vpd, p1 = safe_spearman(det, chi_sub)
            r_wba_vpd, p2 = safe_spearman(dig, chi_sub)
            r_ck_wba, p3 = safe_spearman(det, dig)
            agree_ck_vpd = float(np.mean((det > 0.5) == (chi_sub > 0.5)))
            per_kperp[f"{kperp:.0e}"] = dict(
                kperp=kperp, V_PD=r["V_PD"], DeltaT=r["DeltaT"],
                chi_mean=float(chi_sub.mean()),
                r_ck_vpd=r_ck_vpd, r_wba_vpd=r_wba_vpd, r_ck_wba=r_ck_wba,
                agree_ck_vpd=agree_ck_vpd)
            print(f"   kperp={kperp:.0e}: V_PD={r['V_PD']:.3f} "
                  f"chi_mean(z=0core)={chi_sub.mean():.3f} "
                  f"r(cKAM,VPD)={r_ck_vpd:.3f} r(WBA,VPD)={r_wba_vpd:.3f} "
                  f"r(cKAM,WBA)={r_ck_wba:.3f} agree(cKAM,VPD)={agree_ck_vpd:.2f}")

        # representative maps (finest kperp) for the visualisation figure
        r_last = solve3d(f, kpar=1.0, kperp=KPERPS[-1], rho_min=0.0, rho_max=1.0,
                         theta_period=2 * np.pi / m, core=CORE, method="amg",
                         x0=x0, **g)
        chi_core, rc, _ = chi_zeta0_on_core(r_last, CORE)
        Tslice = r_last["T"][:, :, 0]
        cmask = (r_last["rho"] >= CORE[0]) & (r_last["rho"] <= CORE[1])

        out[str(m)] = dict(
            m=m, grid=g,
            ck_nonexist=float(diag["detected"].mean()),
            dig_median=float(np.median(diag["dig"])),
            chaos_frac=float(np.mean(diag["dig"] < 5)),
            tc_median=float(np.nanmedian(diag["t_c"])),
            per_kperp=per_kperp,
            # store fields for figures
            dig=diag["dig"].tolist(), detected=diag["detected"].astype(int).tolist(),
            rho_diag=rho_diag.tolist(), theta=theta.tolist(),
            chi_map=chi_core.tolist(), T_map=Tslice[cmask].tolist(),
            rho_core=r_last["rho"][cmask].tolist())
    return out


if __name__ == "__main__":
    out = run()
    with open("results/stage3.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n=== Stage 3: does r(cKAM,VPD) grow as kperp -> 0? ===")
    for m in MS:
        rs = [out[str(m)]["per_kperp"][f"{kp:.0e}"]["r_ck_vpd"] for kp in KPERPS]
        print(f" m={m}: r(cKAM,VPD) at kperp={KPERPS} = "
              f"{[f'{x:.3f}' for x in rs]}")
    print("wrote results/stage3.json")
