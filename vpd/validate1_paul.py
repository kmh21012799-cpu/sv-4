"""Validation 1 on PAUL'S field (envelope) -- single island, reproduce Paul Fig 1.
Same gate as before: V_PD turns positive near eps_crit ~ sqrt(kperp), scales as
sqrt(eps). (m,n)=(2,1) -> resonance at rho=0.5; domain [0,1], envelope makes the
perturbation vanish at the boundaries.)"""
import json
import numpy as np
from vpd.paulsolve import paul_single_island
from vpd.solver2d import solve_single_island

M, N = 2, 1
RHO_MIN, RHO_MAX = 0.0, 1.0
CORE = (0.1, 0.9)
NRHO, NU = 257, 64


def run(kperp_list=(1e-2, 1e-4, 1e-6), n_eps=16):
    out = {}
    for kperp in kperp_list:
        ref = np.sqrt(kperp)            # order-of-magnitude anchor
        eps_grid = np.logspace(np.log10(ref) - 2.5, np.log10(ref) + 2.5, n_eps)
        rows = []
        for eps in eps_grid:
            f = paul_single_island(M, N, eps)
            r = solve_single_island(f, M, N, kpar=1.0, kperp=kperp,
                                    rho_min=RHO_MIN, rho_max=RHO_MAX,
                                    Nrho=NRHO, Nu=NU, core=CORE)
            rows.append(dict(eps=float(eps), V_PD=r["V_PD"], DeltaT=r["DeltaT"],
                             residual=r["residual"],
                             max_principle_ok=r["max_principle_ok"],
                             mp_min=r["mp_min"], mp_max=r["mp_max"]))
            print(f"  kperp={kperp:.0e} eps={eps:.3e} V_PD={r['V_PD']:.4f} "
                  f"resid={r['residual']:.1e} mp=({r['mp_min']:.3f},{r['mp_max']:.3f})")
        out[f"{kperp:.0e}"] = dict(eps_crit=float(np.sqrt(kperp) / 2), rows=rows)
    return out


def analyze(out):
    print("\n=== Validation 1 (Paul field) summary ===")
    for key, d in out.items():
        eps = np.array([r["eps"] for r in d["rows"]])
        vpd = np.array([r["V_PD"] for r in d["rows"]])
        ec = d["eps_crit"]
        idx = np.where(vpd > 0.02)[0]
        eps_tr = eps[idx[0]] if len(idx) else np.nan
        band = (vpd > 0.05) & (vpd < 0.6)
        slope = np.polyfit(np.log(eps[band]), np.log(vpd[band]), 1)[0] if band.sum() >= 3 else np.nan
        print(f" kperp={key}: eps_crit=sqrt(k)/2={ec:.2e} eps_trans={eps_tr:.2e} "
              f"ratio={eps_tr/ec:.2f} (within 1 decade: {0.1<eps_tr/ec<10}) "
              f"slope(~0.5)={slope:.2f}")


if __name__ == "__main__":
    out = run()
    analyze(out)
    json.dump(out, open("results/validation1_paul.json", "w"), indent=2)
    print("wrote results/validation1_paul.json")
