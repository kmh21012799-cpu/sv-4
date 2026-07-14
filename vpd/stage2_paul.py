"""Stage 2 on PAUL's field: V_PD and DeltaT for the four original fields
(paul_m4/12/20/36, envelope included), core [0.25,0.75], kperp=1e-4/1e-5/1e-6.
Saves grid data (chi, T at zeta=0) to results/grid/."""
import json
import os
import numpy as np
from vpd.paulsolve import paul_critical
from vpd.solver3d import solve3d

MS = [4, 12, 20, 36]
KPERPS = [1e-4, 1e-5, 1e-6]
CORE = (0.25, 0.75)
os.makedirs("results/grid", exist_ok=True)


def grid_for(field, m):
    n_max = max(n for (mm, n, e) in field.modes)
    Nzeta = int(2 ** np.ceil(np.log2(max(32, 2.6 * n_max))))
    return dict(Nrho=129, Ntheta=32, Nzeta=Nzeta)


def run():
    out = {}
    for m in MS:
        f = paul_critical(m)
        g = grid_for(f, m)
        x0 = None; rows = []
        for kperp in KPERPS:
            r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / m, core=CORE, method="amg",
                        x0=x0, **g)
            x0 = r["T"]
            rows.append(dict(kperp=float(kperp), V_PD=r["V_PD"], DeltaT=r["DeltaT"],
                             residual=r["residual"], n_iter=r["n_iter"],
                             mp_min=r["mp_min"], mp_max=r["mp_max"]))
            print(f"  m={m:2d} kperp={kperp:.0e} V_PD={r['V_PD']:.4f} "
                  f"DeltaT={r['DeltaT']:.4f} it={r['n_iter']} "
                  f"resid={r['residual']:.1e} mp=({r['mp_min']:.3f},{r['mp_max']:.3f})")
            # save grid chi(zeta=0) for the finest kperp (for Stage 3 + archive)
            if kperp == KPERPS[-1]:
                np.savez_compressed(f"results/grid/vpd_m{m}_kperp{kperp:.0e}.npz",
                                    chi_zeta0=r["chi"][:, :, 0], rho=r["rho"],
                                    theta=np.arange(g["Ntheta"]) * (2*np.pi/m)/g["Ntheta"],
                                    T_zeta0=r["T"][:, :, 0])
        out[str(m)] = dict(m=m, modes=[[int(a),int(b),float(c)] for a,b,c in f.modes],
                           grid=g, D_QL=f.quasilinear_D(), rows=rows)
    return out


if __name__ == "__main__":
    out = run()
    json.dump(out, open("results/stage2_paul.json", "w"), indent=2)
    print("\n=== Stage 2 (Paul field): V_PD | DeltaT ===")
    print(f"{'kperp':>8}", *[f"  m={m:<10d}" for m in MS])
    for i, kp in enumerate(KPERPS):
        cells = [f" {out[str(m)]['rows'][i]['V_PD']:.3f}|{out[str(m)]['rows'][i]['DeltaT']:.3f}" for m in MS]
        print(f"{kp:8.0e}", *cells)
    print("\nChecks (Paul sec 7):")
    for i, kp in enumerate(KPERPS):
        vpd = {m: out[str(m)]['rows'][i]['V_PD'] for m in MS}
        dT = {m: out[str(m)]['rows'][i]['DeltaT'] for m in MS}
        print(f" kperp={kp:.0e}: argmax V_PD=m{max(vpd,key=vpd.get)} (expect 4); "
              f"argmax DeltaT=m{max(dT,key=dT.get)} (expect 36); "
              f"argmin DeltaT=m{min(dT,key=dT.get)} (expect 4)")
    print("wrote results/stage2_paul.json + results/grid/*.npz")
