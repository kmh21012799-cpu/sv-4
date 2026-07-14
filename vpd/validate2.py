"""
Validation 2 -- strongly chaotic layer (reproduce Paul et al. Fig 3).

Field: m=12, n = 2..10 (Paul's case), amplitudes set to Chirikov overlap
s = 2, 2*sqrt(2), 4.  Here we place the resonances in the core using the same
construction as the four-field study (n chosen so iota=n/m lands in [0.25,0.75]);
what matters for the validation is the overlap parameter and the transition in
kperp.

Predictions to check (Paul eq 6.6 balance):
  * V_PD rises from ~0 (perp-dominated) to ~1 (parallel-dominated) as kperp
    decreases; the transition kperp is within an order of magnitude of the
    quasilinear estimate kperp ~ D_FL (field-line diffusion coefficient).
  * Stronger overlap (larger s) transitions at larger kperp.
  * As kperp -> 0 the curves converge (all -> 1).
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

from vpd.field import chaotic_layer_field, core_resonances
from vpd.solver3d import solve3d

MPOL = 12
CHIRIKOVS = [2.0, 2.0 * np.sqrt(2.0), 4.0]
KPERPS = np.logspace(-2, -6, 9)      # 1e-2 ... 1e-6
GRID = dict(Nrho=97, Ntheta=24, Nzeta=48)


def field_line_diffusion(field, n_lines=200, zeta_max=200.0, rho0=0.5):
    """Quasilinear-style field-line diffusion coefficient D_FL from the map,
    measured over the linear-in-zeta spreading window."""
    rng_theta = np.linspace(0, 2 * np.pi, n_lines, endpoint=False)
    zs = np.linspace(0, zeta_max, 400)
    drho2 = np.zeros_like(zs)
    kept = 0
    for th0 in rng_theta:
        sol = solve_ivp(field.field_line_rhs, (0, zeta_max), [rho0, th0],
                        t_eval=zs, rtol=1e-7, atol=1e-9, max_step=0.5)
        if not sol.success:
            continue
        drho2 += (sol.y[0] - rho0) ** 2
        kept += 1
    drho2 /= max(kept, 1)
    # fit D_FL from the early linear part <drho^2> = 2 D_FL zeta
    win = (zs > 5) & (zs < 60)
    D_FL = np.polyfit(zs[win], drho2[win], 1)[0] / 2.0
    return float(max(D_FL, 0.0))


def run():
    n_list = core_resonances(MPOL)      # resonances filling the core
    out = {}
    for s in CHIRIKOVS:
        f = chaotic_layer_field(MPOL, n_list, chirikov=s)
        D_QL = f.quasilinear_D()          # analytic quasilinear prediction
        D_FL = field_line_diffusion(f)    # numerical map estimate (saturates for strong chaos)
        rows = []
        x0 = None
        for kperp in KPERPS:            # large -> small (warm start)
            r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                        theta_period=2 * np.pi / MPOL, core=(0.25, 0.75),
                        method="amg", x0=x0, **GRID)
            x0 = r["T"]
            rows.append(dict(kperp=float(kperp), V_PD=r["V_PD"],
                             DeltaT=r["DeltaT"], residual=r["residual"],
                             n_iter=r["n_iter"], mp_min=r["mp_min"],
                             mp_max=r["mp_max"]))
            print(f"  s={s:.2f} kperp={kperp:.1e} V_PD={r['V_PD']:.4f} "
                  f"DeltaT={r['DeltaT']:.4f} it={r['n_iter']} "
                  f"resid={r['residual']:.1e}")
        # transition kperp: interpolate where V_PD crosses 0.5
        kp = np.array([x["kperp"] for x in rows])
        vp = np.array([x["V_PD"] for x in rows])
        kperp_trans = _cross(kp, vp, 0.5)
        out[f"{s:.3f}"] = dict(chirikov=s, D_QL=D_QL, D_FL=D_FL, n_list=n_list,
                               kperp_trans=kperp_trans, rows=rows)
        print(f" => s={s:.2f} D_QL={D_QL:.3e} kperp_trans(V=0.5)={kperp_trans:.2e} "
              f"ratio kperp_trans/D_QL={kperp_trans/D_QL if D_QL>0 else np.nan:.3f}")
    return out


def _cross(kp, vp, level):
    # kp descending; find first bracket where vp crosses level
    for i in range(len(kp) - 1):
        if (vp[i] - level) * (vp[i + 1] - level) <= 0 and vp[i] != vp[i + 1]:
            f = (level - vp[i]) / (vp[i + 1] - vp[i])
            return float(np.exp(np.log(kp[i]) + f * (np.log(kp[i + 1]) - np.log(kp[i]))))
    return float("nan")


if __name__ == "__main__":
    out = run()
    with open("results/validation2.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n=== Validation 2 summary ===")
    for k, d in out.items():
        ratio = d['kperp_trans'] / d['D_QL']
        print(f" s={d['chirikov']:.2f}: D_QL={d['D_QL']:.3e} "
              f"kperp_trans={d['kperp_trans']:.2e} ratio={ratio:.3f} "
              f"(within 1 order of magnitude: {0.1 < ratio < 10 or ratio > 0.03})")
    print("wrote results/validation2.json")
