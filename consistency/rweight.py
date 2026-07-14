"""
Secondary check: does a toroidal R-weighting of the volume element change V_PD?

V_PD = int Theta(...) w dV / int w dV, over the core.  Compare:
  * uniform  w = 1                         (slab / Paul's constant-Jacobian model)
  * toroidal w = 1 + eps_a cos(theta)      (large-aspect-ratio torus, R = R0(1+eps_a cos th))

Reuses the C3b solver's chi field; no PDE re-solve beyond one reference solve.
"""
import sys
import numpy as np
from vpd.field import paul_field
from vpd.solver3d import solve3d


def run(m=12, kperp=1e-6, eps_a=0.3):
    f = paul_field(m, chirikov=1.0)
    from vpd.stage2 import grid_for
    g = grid_for(m)
    r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                theta_period=2 * np.pi / m, core=(0.25, 0.75), method="amg", **g)
    chi = r["chi"]                      # (Nrho, Ntheta, Nzeta)
    rho = r["rho"]
    cmask = (rho >= 0.25) & (rho <= 0.75)
    chi_c = chi[cmask]                  # (Ncore, Ntheta, Nzeta)
    Ntheta = chi_c.shape[1]
    # theta is the REDUCED poloidal angle in [0, 2pi/m); the physical poloidal
    # angle for the R-weight is Theta = m * theta_reduced folded, but the torus
    # R-weight uses the PHYSICAL poloidal angle. In the reduced cell the physical
    # angle sweeps [0, 2pi) as theta_reduced sweeps [0, 2pi/m), so use phys angle.
    theta_phys = np.linspace(0, 2 * np.pi, Ntheta, endpoint=False)
    w_uni = np.ones(Ntheta)
    w_tor = 1.0 + eps_a * np.cos(theta_phys)

    def vpd(w):
        wf = w[None, :, None]
        return float((chi_c * wf).sum() / (np.ones_like(chi_c) * wf).sum())

    v_uni = vpd(w_uni)
    v_tor = vpd(w_tor)
    print(f"m={m} kperp={kperp:g} eps_a={eps_a}: "
          f"V_PD uniform={v_uni:.4f}  V_PD torus(1+{eps_a}cos)={v_tor:.4f}  "
          f"rel change={100*(v_tor-v_uni)/v_uni:+.2f}%")
    return dict(m=m, kperp=kperp, eps_a=eps_a, vpd_uniform=v_uni, vpd_torus=v_tor,
                rel_change_pct=100 * (v_tor - v_uni) / v_uni)


if __name__ == "__main__":
    import json
    res = [run(m=int(x)) for x in (sys.argv[1:] or [12])]
    json.dump(res, open("results/rweight.json", "w"), indent=2)
