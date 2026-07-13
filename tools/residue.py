"""
Greene's residue for periodic field lines of the KMM / Paul flows.

A period-q orbit of the Poincare map P (advance the toroidal angle by 2pi)
is a fixed point of P^q.  Its stability is set by the 2x2 monodromy (tangent)
matrix M obtained by integrating the *variational* equations along the orbit:

        d(eta)/d(time) = J(s, time) . eta ,      J = D(RHS)/D(state)

Greene's residue:
        R = (2 - tr M) / 4
    R in (0,1)  : elliptic  (island O-point)
    R < 0 or >1 : hyperbolic (X-point)
    R = 0 or 1  : parabolic / period-doubling thresholds

Reference: J.M. Greene, J. Math. Phys. 20, 1183 (1979); MacKay's converse-KAM
lectures.  See tools/SOURCE.txt for provenance of this (re)implementation.
"""

import numpy as np
from scipy.integrate import solve_ivp


def _augmented_rhs(field, period_len):
    """Return f(t, Y) integrating state (2) + tangent matrix columns (4)."""
    def rhs(t, Y):
        psi, th = Y[0], Y[1]
        dpsi, dth = field.rhs(t, [psi, th])
        # Jacobian of (dpsi/dt, dth/dt) w.r.t (psi, th) by central differences.
        h = 1e-7
        f_pp = field.rhs(t, [psi + h, th]); f_pm = field.rhs(t, [psi - h, th])
        f_tp = field.rhs(t, [psi, th + h]); f_tm = field.rhs(t, [psi, th - h])
        J = np.array([
            [(f_pp[0] - f_pm[0]) / (2 * h), (f_tp[0] - f_tm[0]) / (2 * h)],
            [(f_pp[1] - f_pm[1]) / (2 * h), (f_tp[1] - f_tm[1]) / (2 * h)],
        ])
        Mcols = Y[2:6].reshape(2, 2)
        dM = J @ Mcols
        return [dpsi, dth, dM[0, 0], dM[0, 1], dM[1, 0], dM[1, 1]]
    return rhs


def monodromy(field, psi0, theta0, q, rtol=1e-11, atol=1e-13):
    """Tangent map M and end state after q toroidal turns from (psi0, theta0)."""
    Y0 = [psi0, theta0, 1.0, 0.0, 0.0, 1.0]
    T = 2.0 * np.pi * q
    sol = solve_ivp(_augmented_rhs(field, T), (0.0, T), Y0,
                    rtol=rtol, atol=atol, method="DOP853")
    end = sol.y[:, -1]
    M = end[2:6].reshape(2, 2)
    return M, np.array([end[0], end[1]])


def poincare_q(field, psi0, theta0, q, rtol=1e-11, atol=1e-13):
    """Advance q toroidal turns; return (psi, theta)."""
    y = [psi0, theta0]
    for _ in range(q):
        sol = solve_ivp(field.rhs, (0.0, 2.0 * np.pi), y,
                        rtol=rtol, atol=atol, method="DOP853")
        y = [sol.y[0, -1], sol.y[1, -1]]
    return np.array(y)


def refine_periodic(field, psi0, theta0, q, fix_theta=True, iters=40, tol=1e-13):
    """Newton-refine a period-q orbit.

    For orbits sitting on a symmetry line (theta = 0), the theta coordinate is
    pinned by symmetry, so we solve the 1-D equation psi_return(psi) - psi = 0
    (fix_theta=True).  Otherwise solve the full 2-D fixed-point of P^q.
    """
    if fix_theta:
        p = psi0
        for _ in range(iters):
            end = poincare_q(field, p, theta0, q)
            g = end[0] - p
            if abs(g) < tol:
                break
            h = 1e-8
            endp = poincare_q(field, p + h, theta0, q)
            dg = (endp[0] - (p + h) - g) / h
            p = p - g / dg
        return p, theta0
    # full 2-D Newton
    x = np.array([psi0, theta0])
    for _ in range(iters):
        end = poincare_q(field, x[0], x[1], q)
        g = end - x
        if np.linalg.norm(g) < tol:
            break
        h = 1e-8
        Jf = np.empty((2, 2))
        for k in range(2):
            xp = x.copy(); xp[k] += h
            ep = poincare_q(field, xp[0], xp[1], q)
            Jf[:, k] = (ep - xp - g) / h
        x = x - np.linalg.solve(Jf, g)
    return x[0], x[1]


def greene_residue(field, psi0, theta0, q):
    """Residue R = (2 - tr M)/4 of the period-q orbit through (psi0, theta0)."""
    M, _ = monodromy(field, psi0, theta0, q)
    return (2.0 - np.trace(M)) / 4.0
