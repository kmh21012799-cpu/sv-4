"""Graded-t_c correlations on Paul's field (converse-KAM is fully saturated, so
the binary detected<->chi correlation is undefined for m>=12; the graded t_c
still carries spatial variance). Reuses saved grid data; re-solves only the
fast kperp=1e-4/1e-5 chi maps (1e-6 loaded from Stage 3)."""
import json
import numpy as np
from scipy.stats import spearmanr
from vpd.paulsolve import paul_critical
from vpd.solver3d import solve3d
from vpd.stage3_paul import grid_for, _sample

MS = [4, 12, 20, 36]
KPERPS = [1e-4, 1e-5, 1e-6]


def sp(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5 or np.std(a[ok]) < 1e-12 or np.std(b[ok]) < 1e-12:
        return float("nan")
    return float(spearmanr(a[ok], b[ok])[0])


def run():
    out = {}
    for m in MS:
        d = np.load(f"results/grid/stage3_m{m}.npz")
        tc, dig, chi_1e6 = d["tc"], d["dig"], d["chi"]
        rho, theta = d["rho"], d["theta"]
        f = paul_critical(m)
        g = grid_for(f, m)
        th_solve = np.arange(g["Ntheta"]) * (2 * np.pi / m) / g["Ntheta"]
        chis = {1e-6: chi_1e6}
        x0 = None
        for kp in [1e-4, 1e-5]:                 # re-solve the cheap ones
            r = solve3d(f, kpar=1.0, kperp=kp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / m, core=(0.25, 0.75),
                        method="amg", x0=x0, **g)
            x0 = r["T"]
            chis[kp] = _sample(r["chi"][:, :, 0], r["rho"], th_solve, rho, theta)
        per = {}
        for kp in KPERPS:
            chi = chis[kp]
            # graded: earlier detection (smaller t_c) ~ "more destroyed"; use -t_c
            per[f"{kp:.0e}"] = dict(
                r_tc_vpd=sp(-tc, chi),           # graded converse-KAM <-> V_PD
                r_wba_vpd=sp(dig, chi),
                chi_mean=float(np.nanmean(chi)))
            print(f" m={m} kperp={kp:.0e}: r(-t_c,VPD)={per[f'{kp:.0e}']['r_tc_vpd']:.3f} "
                  f"r(WBA,VPD)={per[f'{kp:.0e}']['r_wba_vpd']:.3f}")
        out[str(m)] = dict(m=m, r_tc_wba=sp(-tc, dig), per_kperp=per)
        print(f" m={m}: r(-t_c,WBA)={out[str(m)]['r_tc_wba']:.3f}")
    return out


if __name__ == "__main__":
    out = run()
    json.dump(out, open("results/correlations_paul.json", "w"), indent=2)
    print("\n=== key question #8: does graded r(-t_c,VPD) grow as kperp->0? ===")
    for m in MS:
        rs = [out[str(m)]["per_kperp"][f"{kp:.0e}"]["r_tc_vpd"] for kp in KPERPS]
        print(f" m={m}: {[f'{x:.3f}' for x in rs]}")
    print("wrote results/correlations_paul.json")
