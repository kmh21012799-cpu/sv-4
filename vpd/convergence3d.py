"""3D grid convergence: V_PD and DeltaT must converge under refinement, and the
max-principle overshoot must shrink (guards the Heaviside/conditioning traps)."""
import numpy as np
from vpd.field import chaotic_layer_field, core_resonances
from vpd.solver3d import solve3d

f = chaotic_layer_field(12, core_resonances(12), chirikov=2.0)
kperp = 1e-5
print(f"3D grid convergence: m=12 s=2 kperp={kperp:g}")
for (Nr, Nt, Nz) in [(65, 16, 32), (97, 24, 48), (129, 32, 64), (161, 40, 80)]:
    r = solve3d(f, kpar=1.0, kperp=kperp, rho_min=0.0, rho_max=1.0,
                theta_period=2 * np.pi / 12, Nrho=Nr, Ntheta=Nt, Nzeta=Nz,
                core=(0.25, 0.75), method="amg")
    over = max(-r["mp_min"], r["mp_max"] - 1)
    print(f"  ({Nr:3d},{Nt:2d},{Nz:2d}) N={Nr*Nt*Nz:7d}: V_PD={r['V_PD']:.4f} "
          f"DeltaT={r['DeltaT']:.4f} overshoot={over:.2e} "
          f"resid={r['residual']:.1e} it={r['n_iter']}")
