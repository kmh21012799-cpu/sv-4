"""Grid convergence for Validation 1: V_PD and the max-principle overshoot
must converge as the grid is refined (guards against the Heaviside/island-area
non-monotone error that burned C2)."""
import numpy as np
from vpd.field import single_island_field
from vpd.solver2d import solve_single_island

M, N = 2, 1
kperp = 1e-4
eps = 2.706e-2  # a point on the rising part of the curve
f = single_island_field(M, N, eps)
print(f"Grid convergence: (m,n)=({M},{N}) kperp={kperp:g} eps={eps:g}")
for (Nrho, Nu) in [(129, 32), (193, 48), (257, 64), (385, 96), (513, 128)]:
    r = solve_single_island(f, M, N, kpar=1.0, kperp=kperp,
                            rho_min=-0.6, rho_max=0.6, Nrho=Nrho, Nu=Nu,
                            core=(-0.5, 0.5))
    print(f"  Nrho={Nrho:4d} Nu={Nu:4d}: V_PD={r['V_PD']:.5f} "
          f"DeltaT={r['DeltaT']:.5f} overshoot={max(-r['mp_min'],r['mp_max']-1):.2e} "
          f"resid={r['residual']:.1e}")
