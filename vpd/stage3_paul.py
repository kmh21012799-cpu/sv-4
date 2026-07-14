"""
Stage 3 on PAUL's field, three axes with the ORIGINAL C2/C3a code:
  converse-KAM : consistency/orig/converse_kam.detection_map_uv  (C2, 393328e)
  WBA          : consistency/orig/wba.dig_map_rk4                 (C3a, bf5ca06)
  V_PD         : this repo's solver on Paul's field

Same field (original PaulField), same (rho,theta) grid, core [0.25,0.75].
Saves all grid data to results/grid/ (committed, not gitignored).
"""
import os, sys, json, time
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "consistency", "orig"))
import field_paul as fp           # ORIGINAL Paul field
import converse_kam as ck         # ORIGINAL C2
import wba as wba_mod             # ORIGINAL C3a

from vpd.paulsolve import PaulSolveField
from vpd.solver3d import solve3d

MS = [4, 12, 20, 36]
KPERPS = [1e-4, 1e-5, 1e-6]
CORE = (0.25, 0.75)
os.makedirs("results/grid", exist_ok=True)


def grid_for(field, m):
    n_max = max(n for (mm, n, e) in field.modes)
    Nzeta = int(2 ** np.ceil(np.log2(max(32, 2.6 * n_max))))
    return dict(Nrho=129, Ntheta=32, Nzeta=Nzeta)


def sp(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if np.std(a[ok]) < 1e-12 or np.std(b[ok]) < 1e-12:
        return float("nan")
    return float(spearmanr(a[ok], b[ok])[0])


def _sample(chi, rho_full, th_solve, rho_grid, theta_grid):
    out = np.zeros((len(rho_grid), len(theta_grid)))
    for i, rr in enumerate(rho_grid):
        ir = int(np.argmin(np.abs(rho_full - rr)))
        for j, tt in enumerate(theta_grid):
            jt = int(np.argmin(np.abs(th_solve - tt)))
            out[i, j] = chi[ir, jt]
    return out


def chi_scan(m, rho_grid, theta_grid):
    """Warm-started kperp scan; returns {kperp: (chi_grid, V_PD, DeltaT)}."""
    f = PaulSolveField({4: fp.paul_m4, 12: fp.paul_m12, 20: fp.paul_m20,
                        36: fp.paul_m36}[m]())
    g = grid_for(f, m)
    th_solve = np.arange(g["Ntheta"]) * (2 * np.pi / m) / g["Ntheta"]
    res = {}
    x0 = None
    for kperp in KPERPS:
        r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                    theta_period=2 * np.pi / m, core=CORE, method="amg",
                    x0=x0, **g)
        x0 = r["T"]
        res[kperp] = (_sample(r["chi"][:, :, 0], r["rho"], th_solve, rho_grid, theta_grid),
                      r["V_PD"], r["DeltaT"])
    return res


def run(n_rho=24, n_theta=24, t_f=200.0, wba_periods=1000):
    out = {}
    for m in MS:
        pf = {4: fp.paul_m4, 12: fp.paul_m12, 20: fp.paul_m20, 36: fp.paul_m36}[m]()
        rho = np.linspace(0.25, 0.75, n_rho)
        theta = np.linspace(0, 2 * np.pi / m, n_theta, endpoint=False)
        # --- ORIGINAL converse-KAM t_c map ---
        t0 = time.time()
        tc = ck.detection_map_uv(pf, rho, theta, t_f=t_f, symmetry=False)
        det = np.isfinite(tc) & (tc <= t_f)
        # --- ORIGINAL WBA dig map (vectorised, native to Paul field) ---
        dig = wba_mod.dig_map_rk4(pf, rho, theta, n_periods=wba_periods,
                                  steps_per_period=16, which="dig_psi")
        print(f"[m={m}] cKAM+WBA maps in {time.time()-t0:.0f}s "
              f"(nonexist={np.nanmean(det):.3f}, dig_med={np.nanmedian(dig):.2f})")
        scan = chi_scan(m, rho, theta)
        per = {}
        for kperp in KPERPS:
            chi, vpd, dT = scan[kperp]
            per[f"{kperp:.0e}"] = dict(
                V_PD=vpd, DeltaT=dT,
                r_ck_vpd=sp(det.astype(float), chi),
                r_wba_vpd=sp(dig, chi),
                chi_mean=float(np.nanmean(chi)))
            if kperp == KPERPS[-1]:
                np.savez_compressed(f"results/grid/stage3_m{m}.npz",
                                    tc=tc, detected=det.astype(int), dig=dig,
                                    chi=chi, rho=rho, theta=theta)
            print(f"   kperp={kperp:.0e}: V_PD={vpd:.3f} DeltaT={dT:.3f} "
                  f"r(cKAM,VPD)={per[f'{kperp:.0e}']['r_ck_vpd']:.3f} "
                  f"r(WBA,VPD)={per[f'{kperp:.0e}']['r_wba_vpd']:.3f}")
        r_ck_wba = sp(det.astype(float), dig)
        out[str(m)] = dict(
            m=m, ck_nonexist=float(np.nanmean(det)),
            tc_median=float(np.nanmedian(tc[det])) if det.any() else float("nan"),
            dig_median=float(np.nanmedian(dig)),
            chaos_frac=float(np.nanmean(dig < 5)),
            r_ck_wba=r_ck_wba, per_kperp=per)
        print(f" m={m}: nonexist={out[str(m)]['ck_nonexist']:.3f} "
              f"t_c_med={out[str(m)]['tc_median']:.1f} dig_med={out[str(m)]['dig_median']:.2f} "
              f"chaos={out[str(m)]['chaos_frac']:.2f} r(cKAM,WBA)={r_ck_wba:.3f}")
    return out


if __name__ == "__main__":
    out = run()
    json.dump(out, open("results/stage3_paul.json", "w"), indent=2)
    print("\n=== Stage 3 (Paul field): key question #8 (does r(cKAM,VPD) grow as kperp->0?) ===")
    for m in MS:
        rs = [out[str(m)]["per_kperp"][f"{kp:.0e}"]["r_ck_vpd"] for kp in KPERPS]
        print(f" m={m}: {[f'{x:.3f}' for x in rs]}")
    print("wrote results/stage3_paul.json + results/grid/stage3_m*.npz")
