"""
Validation 1 -- single island chain (reproduce Paul et al. Fig 1).

Predictions to check:
  * V_PD turns positive as eps increases through eps_crit = sqrt(kperp)/2
    (agreement expected to within an order of magnitude, not exact).
  * At large eps, V_PD grows like sqrt(eps) (separatrix volume fraction),
    until the island saturates the measurement band.
"""
import json
import numpy as np
from vpd.field import single_island_field
from vpd.solver2d import solve_single_island

M, N = 2, 1                      # resonance iota = 1/2 at rho_* = 0
RHO_MIN, RHO_MAX = -0.6, 0.6
CORE = (-0.5, 0.5)
NRHO, NU = 257, 64


def run(kperp_list=(1e-2, 1e-4, 1e-6), n_eps=16):
    out = {}
    for kperp in kperp_list:
        eps_crit = np.sqrt(kperp) / 2.0
        # scan ~2 decades below to ~2 decades above eps_crit
        eps_grid = np.logspace(np.log10(eps_crit) - 2.2,
                               np.log10(eps_crit) + 2.2, n_eps)
        rows = []
        for eps in eps_grid:
            f = single_island_field(M, N, eps)
            r = solve_single_island(f, M, N, kpar=1.0, kperp=kperp,
                                    rho_min=RHO_MIN, rho_max=RHO_MAX,
                                    Nrho=NRHO, Nu=NU, core=CORE)
            rows.append(dict(eps=float(eps), V_PD=r["V_PD"], DeltaT=r["DeltaT"],
                             residual=r["residual"],
                             max_principle_ok=r["max_principle_ok"],
                             mp_min=r["mp_min"], mp_max=r["mp_max"]))
            print(f"  kperp={kperp:.0e} eps={eps:.3e} "
                  f"V_PD={r['V_PD']:.4f} resid={r['residual']:.1e} "
                  f"mp=({r['mp_min']:.4f},{r['mp_max']:.4f})")
        out[f"{kperp:.0e}"] = dict(eps_crit=float(eps_crit), rows=rows)
    return out


def analyze(out):
    print("\n=== Validation 1 summary ===")
    for key, d in out.items():
        eps_crit = d["eps_crit"]
        rows = d["rows"]
        eps = np.array([r["eps"] for r in rows])
        vpd = np.array([r["V_PD"] for r in rows])
        # transition eps: first eps where V_PD exceeds a small threshold
        thr = 0.02
        idx = np.where(vpd > thr)[0]
        eps_trans = eps[idx[0]] if len(idx) else np.nan
        ratio = eps_trans / eps_crit if np.isfinite(eps_trans) else np.nan
        # large-eps slope in log-log (fit where V_PD in [0.05,0.6])
        band = (vpd > 0.05) & (vpd < 0.6)
        slope = np.nan
        if band.sum() >= 3:
            slope = np.polyfit(np.log(eps[band]), np.log(vpd[band]), 1)[0]
        print(f" kperp={key}: eps_crit={eps_crit:.3e} eps_trans={eps_trans:.3e} "
              f"ratio={ratio:.2f} (within 1 order of mag: "
              f"{0.1 < ratio < 10 if np.isfinite(ratio) else False}) "
              f"loglog slope(~0.5 expected)={slope:.2f}")
    return out


if __name__ == "__main__":
    out = run()
    analyze(out)
    with open("results/validation1.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote results/validation1.json")
