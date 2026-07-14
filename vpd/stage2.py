"""
Stage 2 -- V_PD and DeltaT for Paul's four fields (m = 4, 12, 20, 36) at
critical overlap (Chirikov s = 1), core rho in [0.25, 0.75], scanning
kperp = 1e-4, 1e-5, 1e-6 (Paul Fig 5, 6).

Predictions to check (Paul section 7):
  * V_PD largest for m=4.
  * DeltaT smallest for m=4 (leaks most) and large for m=36 (insulates even
    though its flux surfaces are the most broken).
  * As kperp -> 0 the four fields' V_PD approach one another.
"""
import json
import numpy as np
from vpd.field import paul_field, core_resonances
from vpd.solver3d import solve3d

MS = [4, 12, 20, 36]
KPERPS = [1e-4, 1e-5, 1e-6]
CORE = (0.25, 0.75)


def grid_for(m):
    """Resolution: Ntheta fixed (phi=m*theta reduction); Nzeta scales with the
    largest resonance n (~1.25 m); Nrho fine enough for the boundary layer."""
    n_max = max(core_resonances(m))
    Nzeta = int(2 ** np.ceil(np.log2(max(32, 2.6 * n_max))))
    return dict(Nrho=129, Ntheta=32, Nzeta=Nzeta)


def run():
    out = {}
    for m in MS:
        f = paul_field(m, chirikov=1.0)
        g = grid_for(m)
        x0 = None
        rows = []
        for kperp in KPERPS:      # large -> small, warm start
            r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / m, core=CORE, method="amg",
                        x0=x0, **g)
            x0 = r["T"]
            rows.append(dict(kperp=float(kperp), V_PD=r["V_PD"],
                             DeltaT=r["DeltaT"], residual=r["residual"],
                             n_iter=r["n_iter"], mp_min=r["mp_min"],
                             mp_max=r["mp_max"]))
            print(f"  m={m:2d} kperp={kperp:.0e} grid={g['Nzeta']:3d}z "
                  f"V_PD={r['V_PD']:.4f} DeltaT={r['DeltaT']:.4f} "
                  f"it={r['n_iter']} resid={r['residual']:.1e} "
                  f"mp=({r['mp_min']:.3f},{r['mp_max']:.3f})")
        out[str(m)] = dict(m=m, n_list=core_resonances(m), grid=g,
                           D_QL=f.quasilinear_D(), rows=rows)
    return out


if __name__ == "__main__":
    out = run()
    with open("results/stage2.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n=== Stage 2 summary (V_PD | DeltaT) ===")
    print(f"{'kperp':>8}", *[f" m={m:<12d}" for m in MS])
    for i, kp in enumerate(KPERPS):
        cells = []
        for m in MS:
            r = out[str(m)]["rows"][i]
            cells.append(f" {r['V_PD']:.3f}|{r['DeltaT']:.3f}   ")
        print(f"{kp:8.0e}", *cells)
    print("\nChecks:")
    for i, kp in enumerate(KPERPS):
        vpd = {m: out[str(m)]["rows"][i]["V_PD"] for m in MS}
        dT = {m: out[str(m)]["rows"][i]["DeltaT"] for m in MS}
        print(f" kperp={kp:.0e}: argmax V_PD = m={max(vpd, key=vpd.get)} "
              f"(expect 4); argmin DeltaT = m={min(dT, key=dT.get)} (expect 4)")
    # convergence of the four fields as kperp -> 0
    spread = [max(out[str(m)]['rows'][i]['V_PD'] for m in MS) -
              min(out[str(m)]['rows'][i]['V_PD'] for m in MS)
              for i in range(len(KPERPS))]
    print(f" V_PD spread across fields vs kperp {KPERPS}: "
          f"{[f'{s:.3f}' for s in spread]} (expect decreasing)")
    print("wrote results/stage2.json")
