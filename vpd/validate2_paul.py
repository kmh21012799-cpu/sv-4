"""Validation 2 on PAUL's field (envelope) -- strongly chaotic layer, Paul Fig 3.
m=12, n=2..10 (Paul's set), Chirikov S=2,2sqrt2,4. Transition kperp should track
the quasilinear D_QL (with the envelope amplitude at each resonance)."""
import json
import numpy as np
from vpd.paulsolve import paul_chaotic_layer
from vpd.solver3d import solve3d

MPOL = 12
NLIST = list(range(2, 11))            # Paul's n = 2..10
CHIRIKOVS = [2.0, 2.0 * np.sqrt(2.0), 4.0]
KPERPS = np.logspace(-2, -6, 9)
GRID = dict(Nrho=97, Ntheta=24, Nzeta=32)   # n<=10 -> Nzeta=32 ok


def _cross(kp, vp, level):
    for i in range(len(kp) - 1):
        if (vp[i] - level) * (vp[i + 1] - level) <= 0 and vp[i] != vp[i + 1]:
            f = (level - vp[i]) / (vp[i + 1] - vp[i])
            return float(np.exp(np.log(kp[i]) + f * (np.log(kp[i + 1]) - np.log(kp[i]))))
    return float("nan")


def run():
    out = {}
    for s in CHIRIKOVS:
        f = paul_chaotic_layer(MPOL, NLIST, chirikov=s)
        D_QL = f.quasilinear_D()
        rows = []; x0 = None
        for kperp in KPERPS:
            r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / MPOL, core=(0.25, 0.75),
                        method="amg", x0=x0, **GRID)
            x0 = r["T"]
            rows.append(dict(kperp=float(kperp), V_PD=r["V_PD"], DeltaT=r["DeltaT"],
                             residual=r["residual"], n_iter=r["n_iter"]))
            print(f"  S={s:.2f} kperp={kperp:.1e} V_PD={r['V_PD']:.4f} "
                  f"DeltaT={r['DeltaT']:.4f} it={r['n_iter']} resid={r['residual']:.1e}")
        kp = np.array([x["kperp"] for x in rows]); vp = np.array([x["V_PD"] for x in rows])
        kt = _cross(kp, vp, 0.5)
        out[f"{s:.3f}"] = dict(chirikov=s, D_QL=D_QL, kperp_trans=kt, rows=rows)
        print(f" => S={s:.2f} D_QL={D_QL:.3e} kperp_trans={kt:.2e} ratio={kt/D_QL:.3f}")
    return out


if __name__ == "__main__":
    out = run()
    json.dump(out, open("results/validation2_paul.json", "w"), indent=2)
    print("\n=== Validation 2 (Paul field) summary ===")
    for k, d in out.items():
        print(f" S={d['chirikov']:.2f}: D_QL={d['D_QL']:.3e} kperp_trans={d['kperp_trans']:.2e} "
              f"ratio={d['kperp_trans']/d['D_QL']:.3f} (within 1 decade: "
              f"{0.03 < d['kperp_trans']/d['D_QL'] < 30})")
    print("wrote results/validation2_paul.json")
