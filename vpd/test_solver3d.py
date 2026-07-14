import time
import numpy as np
from vpd.field import Field, chaotic_layer_field
from vpd.solver3d import solve3d

# --- Test A: unperturbed (eps=0) single m=12 mode => linear T, V_PD=0 ---
f0 = Field([(12, 12, 0.0)], label="unperturbed")
t0 = time.time()
r = solve3d(f0, kpar=1.0, kperp=1e-3, rho_min=0.0, rho_max=1.0,
            theta_period=2 * np.pi / 12, Nrho=33, Ntheta=8, Nzeta=16,
            core=(0.25, 0.75), return_fields=True)
rho = r["rho"]
Tlin = (rho - 0.0) / 1.0
lin_err = np.max(np.abs(r["T"] - Tlin[:, None, None]))
print(f"[A unpert] lin-T err={lin_err:.2e} V_PD={r['V_PD']:.4f} "
      f"DeltaT={r['DeltaT']:.4f} resid={r['residual']:.1e} "
      f"maxpr={r['max_principle_ok']} t={time.time()-t0:.1f}s")

# --- Test B: chaotic layer m=12, timing on a realistic grid ---
f1 = chaotic_layer_field(12, list(range(9, 16)), chirikov=2.0)
for (Nrho, Nth, Nze) in [(65, 16, 32), (97, 24, 48)]:
    t0 = time.time()
    r = solve3d(f1, kpar=1.0, kperp=1e-3, rho_min=0.0, rho_max=1.0,
                theta_period=2 * np.pi / 12, Nrho=Nrho, Ntheta=Nth, Nzeta=Nze,
                core=(0.25, 0.75))
    print(f"[B s=2] grid=({Nrho},{Nth},{Nze}) N={Nrho*Nth*Nze} "
          f"V_PD={r['V_PD']:.4f} DeltaT={r['DeltaT']:.4f} "
          f"resid={r['residual']:.1e} maxpr={r['max_principle_ok']} "
          f"mp=({r['mp_min']:.4f},{r['mp_max']:.4f}) t={time.time()-t0:.1f}s")
