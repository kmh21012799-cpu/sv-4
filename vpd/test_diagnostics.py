import time
import numpy as np
from vpd.field import Field, single_island_field, paul_field
from vpd.diagnostics import run_diagnostics

rho = np.linspace(0.28, 0.72, 40)
theta = np.linspace(0, 2 * np.pi / 4, 40, endpoint=False)  # m=4 poloidal period


def summarize(tag, f, th_period_m):
    th = np.linspace(0, 2 * np.pi / th_period_m, 40, endpoint=False)
    t0 = time.time()
    d = run_diagnostics(f, rho, th, n_periods=800, substeps=24, ck_periods=120)
    det = d["detected"]; dig = d["dig"]
    print(f"[{tag}] detected frac={det.mean():.2f}  "
          f"dig: median={np.median(dig):.2f} "
          f"frac(dig<5)={np.mean(dig < 5):.2f} "
          f"frac(dig>8)={np.mean(dig > 8):.2f}  "
          f"t_c median(detected)={np.nanmedian(d['t_c']):.1f}  "
          f"t={time.time()-t0:.1f}s")


# A: unperturbed -> nothing detected, dig high everywhere
summarize("unpert  ", Field([(4, 4, 0.0)]), 4)
# B: single small island (m=4,n=4 at rho=0.5). eps small -> thin island
summarize("sm-island", single_island_field(4, 4, 2e-4), 4)
# C: single large island
summarize("lg-island", single_island_field(4, 4, 5e-3), 4)
# D: Paul m=4 critical overlap (chaotic)
summarize("paul m=4 ", paul_field(4, chirikov=1.0), 4)
# E: Paul m=12
summarize("paul m=12", paul_field(12, chirikov=1.0), 12)
