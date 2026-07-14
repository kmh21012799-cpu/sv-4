import numpy as np
from vpd.field import single_island_field, iota
from vpd.solver2d import solve_single_island

# resonance for (m,n)=(2,1): iota=n/m=0.5 -> rho_*=0. Center domain there.
m, n = 2, 1
rho_min, rho_max = -0.5, 0.5
core = (-0.35, 0.35)

# --- Test A: eps=0 must give linear T and V_PD=0 ---
f0 = single_island_field(m, n, 0.0)
r = solve_single_island(f0, m, n, kpar=1.0, kperp=1e-3, rho_min=rho_min,
                        rho_max=rho_max, Nrho=65, Nu=16, core=core)
T = r["T"]
rho = r["rho"]
Tlin = (rho - rho_min) / (rho_max - rho_min)
lin_err = np.max(np.abs(T - Tlin[:, None]))
print(f"[A eps=0] linear-T max err = {lin_err:.2e}  V_PD={r['V_PD']:.4f} "
      f"resid={r['residual']:.1e} maxprinc={r['max_principle_ok']} "
      f"(min={r['mp_min']:.3f},max={r['mp_max']:.3f})")

# --- Test B: moderate island, check finite V_PD, max principle, residual ---
eps = 0.05
f1 = single_island_field(m, n, eps)
r1 = solve_single_island(f1, m, n, kpar=1.0, kperp=1e-3, rho_min=rho_min,
                         rho_max=rho_max, Nrho=129, Nu=48, core=core)
print(f"[B eps={eps}] V_PD={r1['V_PD']:.4f} DeltaT={r1['DeltaT']:.4f} "
      f"resid={r1['residual']:.1e} maxprinc={r1['max_principle_ok']} "
      f"(min={r1['mp_min']:.3f},max={r1['mp_max']:.3f})")
print("island half-width =", f1.island_half_width(m, n, eps),
      " eps_crit=sqrt(kperp)/2 =", np.sqrt(1e-3) / 2)
