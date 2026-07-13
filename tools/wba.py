"""
Weighted Birkhoff Average (WBA) convergence diagnostic 'dig'.

Reimplemented from the C3a spec (the external `wba` repo is not in this
session's scope). Reference: Das & Yorke; Das, Saiki, Sander, Yorke; and the
Duignan-Meiss application to field lines.

Weight (bump) function, w=1, normalised so int_0^1 g = 1:
    g(s) = C exp(-w / (s(1-s))),  s in (0,1),   0 otherwise,   C = 142.2503757771

Weighted Birkhoff average of an observable h over a window of length tau (in the
flow time zeta), computed as an extra ODE integrated WITH the trajectory:
    dW/dzeta = g(zeta/tau) h(phi_zeta(x0)),     WB = W(tau)/tau
(the substitution s=zeta/tau makes W(tau)/tau = int_0^1 g(s) h ds, a weighted
time-average).

Convergence judge (no true value needed): compare two consecutive windows
    A = WB over [0, tau]     (from x0)
    B = WB over [tau, 2 tau] (from x_tau)
    absdig = -log10 |A - B|
    reldig = -log10 ( |A - B| / max(|A|+|B|)/2, eps) )   (unusable if mean ~ 0)
    dig    = max(absdig, reldig)
Regular (quasiperiodic) orbits: A,B super-converge to the same value -> dig high
(limited by the integration noise floor). Chaotic / near-separatrix (very long
period): A != B -> dig low.
"""
import numpy as np
from scipy.integrate import solve_ivp

WBA_C = 142.2503757771     # C for w=1 so that int_0^1 g(s) ds = 1


def gbump(s, w=1.0):
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    m = (s > 0.0) & (s < 1.0)
    out[m] = WBA_C * np.exp(-w / (s[m] * (1.0 - s[m])))
    return out


def _digits(A, B, eps=1e-300):
    """absdig, reldig, dig for one scalar component."""
    d = abs(A - B)
    absdig = -np.log10(max(d, eps))
    scale = 0.5 * (abs(A) + abs(B))
    reldig = -np.log10(max(d / max(scale, eps), eps)) if scale > eps else absdig
    return absdig, reldig, max(absdig, reldig)


def wba_point(field, rho0, theta0, n_periods=1000, rtol=1e-11, atol=1e-13,
              max_step_frac=0.25):
    """WBA dig at (rho0, theta0) for observables psi=rho and position (x,y).

    Returns dict with dig_psi, dig_pos (=min over x,y), plus the raw A/B.
    tau = 2 pi n_periods (window length in zeta); integrate [0, 2 tau].
    """
    tau = 2.0 * np.pi * n_periods

    def rhs(zeta, y):
        rho, th = y[0], y[1]
        drho, dth = field.rhs(zeta, [rho, th])
        g1 = gbump(zeta / tau)
        g2 = gbump((zeta - tau) / tau)
        hpsi = rho
        hx = rho * np.cos(th)
        hy = rho * np.sin(th)
        return [drho, dth,
                g1 * hpsi, g1 * hx, g1 * hy,
                g2 * hpsi, g2 * hx, g2 * hy]

    y0 = [rho0, theta0, 0, 0, 0, 0, 0, 0]
    sol = solve_ivp(rhs, (0.0, 2.0 * tau), y0, method="DOP853",
                    rtol=rtol, atol=atol, max_step=2.0 * np.pi * max_step_frac)
    W = sol.y[:, -1]
    A = W[2:5] / tau      # [psi, x, y]
    B = W[5:8] / tau
    dpsi = _digits(A[0], B[0])
    dx = _digits(A[1], B[1])
    dy = _digits(A[2], B[2])
    return dict(dig_psi=dpsi[2], absdig_psi=dpsi[0], reldig_psi=dpsi[1],
                dig_pos=min(dx[2], dy[2]), dig_x=dx[2], dig_y=dy[2],
                A=A, B=B)


def dig_psi(field, rho0, theta0, n_periods=1000, **kw):
    return wba_point(field, rho0, theta0, n_periods=n_periods, **kw)["dig_psi"]


# ---- parallel grid map ---------------------------------------------------
def _row(args):
    field, th, rhos, n_periods, guard, which = args
    out = np.full(len(rhos), np.nan)
    for i, r in enumerate(rhos):
        if r <= guard or r >= 1.0 - guard:
            continue
        d = wba_point(field, float(r), float(th), n_periods=n_periods)
        out[i] = d[which]
    return out


def dig_map(field, rho_grid, theta_grid, n_periods=1000, which="dig_psi",
            nproc=None, guard=1e-4):
    """dig[i_rho, j_theta] over (rho, theta)."""
    import os
    import multiprocessing as mp
    if nproc is None:
        nproc = max(1, (os.cpu_count() or 1))
    tasks = [(field, th, rho_grid, n_periods, guard, which) for th in theta_grid]
    if nproc == 1:
        cols = [_row(t) for t in tasks]
    else:
        with mp.Pool(nproc) as pool:
            cols = pool.map(_row, tasks)
    return np.array(cols).T
