"""
converse-KAM 3D  --  MacKay (2018) / Kallinikos-MacKay-Martinez-del-Rio (2023),
Theorem 3.1 (general form) and Theorem 3.2 (symmetry-line short-cut).

Idea: if an invariant torus transverse to the direction field xi passes through a
field line, it forbids an infinitesimal displacement vector from rotating from
one side of the torus to the other.  So if, along the integrated orbit, the
displacement vector rotates in the forbidden way, NO such torus passes through
that orbit.  This is a *non-existence proof*, not a diagnostic -- what Poincare
cannot give.

We use the KMM auxiliary field V (phi = time), so field lines obey
    d(psi,theta)/dphi = (V^psi, V^theta)
and the displacement (tangent) vector obeys the variational equation
    d(eta)/dphi = J . eta ,     J = d(V^psi,V^theta)/d(psi,theta).
Because V^phi = 1, a poloidal eta stays poloidal.

Ingredients (Sec. 3):
  xi   = d/dpsi                      (radial direction field; grad psi || d/dpsi
                                      in the adapted metric)
  beta = flux 2-form; on the poloidal section (V^phi=1)  beta(eta, xi) = -eta^theta
  metric  ds^2 = dpsi^2/(2 B0 psi) + (2 psi/B0) dtheta^2 + R0^2 dphi^2
  lambda(eta) = eta^psi - V^psi * g(V,eta)/g(V,V)      (R^2 factor cancels;
                lambda(V)=0, lambda(xi)>0)

Theorem 3.1 (general): non-existence through s0 once, at some phi,
   (i)  beta(eta, xi) changes sign   AND   (ii) lambda(eta) < 0.
Theorem 3.2 (symmetry line, s0 on theta=0/pi & phi=0/pi): only (i) is needed.

The WHOLE orbit through s0 is then in the non-existence region.
"""

import os
import numpy as np
from scipy.integrate import solve_ivp


def _lambda(field, psi, theta, phi, eta):
    """lambda(eta) = eta^psi - V^psi g(V,eta)/g(V,V) with eta poloidal, V^phi=1."""
    Vpsi, Vth, _ = field.V(psi, theta, phi)
    Vpsi = float(Vpsi); Vth = float(Vth)
    gpp, gtt, gff = field.metric(psi)
    gVeta = gpp * Vpsi * eta[0] + gtt * Vth * eta[1]
    gVV = gpp * Vpsi ** 2 + gtt * Vth ** 2 + gff
    return eta[0] - Vpsi * gVeta / gVV


def _aug_rhs(field):
    def rhs(phi, Y):
        psi, th = Y[0], Y[1]
        Vpsi, Vth, _ = field.V(psi, th, phi)
        J = field.jacobian(psi, th, phi)
        eta = Y[2:4]
        deta = J @ eta
        return [float(Vpsi), float(Vth), deta[0], deta[1]]
    return rhs


def detect_tc(field, psi0, theta0, phi0=0.0, t_f=200.0, symmetry=False,
              rtol=1e-7, atol=1e-9):
    """First converse-KAM detection 'time' phi (>=0) for the orbit through
    (psi0, theta0), or np.inf if none within t_f.

    beta(eta,xi) = -eta^theta ; xi = eta0 = d/dpsi = (1,0).  Detection = first
    zero-crossing of eta^theta (a beta sign change) at which, unless `symmetry`,
    lambda(eta) < 0.  Uses a terminal event so island orbits stop early.
    """
    rhs = _aug_rhs(field)

    def event(phi, Y):
        return Y[3]                    # eta^theta ; zero-crossing = beta sign flip
    event.terminal = True
    event.direction = 0

    # warm-up: advance a little so eta^theta leaves its exact zero start before
    # arming the zero-crossing event (avoids a spurious detection at phi~0).
    # The smallest real detection time is O(15), so a short warm-up is safe.
    phi_warm = 0.5
    sw = solve_ivp(rhs, (0.0, phi_warm), [psi0, theta0, 1.0, 0.0],
                   rtol=rtol, atol=atol, method="DOP853")
    Y = list(sw.y[:, -1])
    phi_start = phi_warm
    guard = 0
    while phi_start < t_f and guard < 200:
        guard += 1
        sol = solve_ivp(rhs, (phi_start, t_f), Y, events=event,
                        rtol=rtol, atol=atol, method="DOP853", max_step=2.0)
        if not sol.success:
            return np.inf
        if sol.t_events[0].size == 0:
            return np.inf              # eta^theta never returned to zero
        phc = float(sol.t_events[0][0])
        ye = sol.y_events[0][0]
        pp, tt, ep, et = ye
        if symmetry or _lambda(field, pp, tt, phc, [ep, et]) < 0.0:
            return phc
        # lambda>=0: step just past the crossing and keep looking
        phi_start = phc + 1e-6
        # advance state slightly past the zero so the event re-arms
        s2 = solve_ivp(rhs, (phc, phc + 1e-3), ye, rtol=rtol, atol=atol,
                       method="DOP853")
        Y = list(s2.y[:, -1])
        phi_start = float(s2.t[-1])
    return np.inf


def _row_worker(args):
    field, z, ygrid, t_f, symmetry = args
    from field_kmm import KMMField
    row = np.full(len(ygrid), np.inf)
    for iy, y in enumerate(ygrid):
        psi0, th0 = KMMField.from_symplectic(np.array(y), np.array(z), B0=field.B0)
        psi0 = float(psi0); th0 = float(th0)
        if psi0 < 1e-4:
            row[iy] = np.nan
            continue
        row[iy] = detect_tc(field, psi0, th0, t_f=t_f, symmetry=symmetry)
    return row


def detection_map(field, ygrid, zgrid, t_f=200.0, symmetry=False, nproc=None):
    """Compute t_c over a grid of symplectic coords (ytil, ztil), row-parallel.

    Returns tc[i,j] (phi at detection or inf) with i over z, j over y so that
    imshow(origin='lower', extent=[y0,y1,z0,z1]) is oriented naturally.
    """
    import multiprocessing as mp
    if nproc is None:
        nproc = max(1, (os.cpu_count() or 1))
    tasks = [(field, z, ygrid, t_f, symmetry) for z in zgrid]
    if nproc == 1:
        rows = [_row_worker(t) for t in tasks]
    else:
        with mp.Pool(nproc) as pool:
            rows = pool.map(_row_worker, tasks)
    return np.array(rows)


def nonexistence_area(tc, ygrid, zgrid, t_f):
    """Area of {t_c <= t_f} in symplectic coords = toroidal flux (paper eq. 12).

    S ~ (cell area) * count.  In (ytil, ztil), d ytil ^ d ztil = d psi ^ d theta,
    so this area IS the enclosed toroidal flux.
    """
    dy = (ygrid[-1] - ygrid[0]) / (len(ygrid) - 1)
    dz = (zgrid[-1] - zgrid[0]) / (len(zgrid) - 1)
    detected = np.isfinite(tc) & (tc <= t_f)
    return float(np.sum(detected) * dy * dz)
