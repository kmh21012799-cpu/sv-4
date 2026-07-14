"""
Poincare sections (section zeta = 0) for Paul's four critical-overlap fields,
m = 4, 12, 20, 36, on the CORRECTED envelope field (tools/field_paul.py).

Recompute (explicitly requested): the underlying orbit points were never saved.
Saves BOTH  paper/data/poincare_m{m}.npz  (committed, reproducibility)
      and   paper/data_matlab/poincare_m{m}.mat  (MATLAB).

Coordinates: rho = psi in [0,1], theta in [0, 2*pi).  kind: 0 = regular seed
(theta=0 line -> KAM curves + island O-points), 1 = chaos-filling orbit.
Seeding mirrors scripts/poincare_paul.py.
"""
import os
import sys
import numpy as np
from scipy.io import savemat

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36   # noqa: E402

FIELDS = {4: paul_m4, 12: paul_m12, 20: paul_m20, 36: paul_m36}

N_SEED = 60          # regular seeds along theta=0
N_CROSS = 220        # crossings per regular seed
N_CHAOS = 6000       # crossings per chaos orbit
N_ORBIT = 4          # chaos orbits (inter-resonance midpoints)


def make(m):
    f = FIELDS[m]()
    rho_r, th_r = [], []
    for p0 in np.linspace(0.01, 0.99, N_SEED):
        psis, ths = f.poincare(p0, 0.0, n_cross=N_CROSS, rtol=1e-8, atol=1e-10)
        rho_r.append(psis); th_r.append(np.mod(ths, 2 * np.pi))
    rho_r = np.concatenate(rho_r); th_r = np.concatenate(th_r)

    res_psi = sorted(n / mm for (mm, n, _) in f.modes)
    mids = [(a + b) / 2 for a, b in zip(res_psi[:-1], res_psi[1:])]
    pick = [mids[int(round(i))] for i in np.linspace(0, len(mids) - 1,
                                                     min(N_ORBIT, len(mids)))]
    rho_c, th_c = [], []
    for p0 in pick:
        psis, ths = f.poincare(p0, 1e-3, n_cross=N_CHAOS, rtol=1e-8, atol=1e-10)
        good = (psis > 0) & (psis < 1)
        rho_c.append(psis[good]); th_c.append(np.mod(ths[good], 2 * np.pi))
    rho_c = np.concatenate(rho_c); th_c = np.concatenate(th_c)

    rho = np.concatenate([rho_r, rho_c])
    theta = np.concatenate([th_r, th_c])
    kind = np.concatenate([np.zeros(rho_r.size, int), np.ones(rho_c.size, int)])
    res = np.array(res_psi)
    npz = f"paper/data/poincare_m{m}.npz"
    np.savez_compressed(npz, rho=rho, theta=theta, kind=kind, res_psi=res)
    savemat(f"paper/data_matlab/poincare_m{m}.mat",
            {"rho": rho.reshape(-1, 1), "theta": theta.reshape(-1, 1),
             "kind": kind.reshape(-1, 1).astype(float),
             "res_psi": res.reshape(-1, 1)}, do_compression=True)
    print(f"  m={m}: {rho.size} pts ({rho_r.size} regular + {rho_c.size} chaos), "
          f"{len(res)} resonances -> {npz}")


if __name__ == "__main__":
    os.makedirs("paper/data", exist_ok=True)
    os.makedirs("paper/data_matlab", exist_ok=True)
    for m in (4, 12, 20, 36):
        make(m)
    print("POINCARE DONE")
