# ★★ SUPERSEDED — DO NOT USE FOR NEW WORK ★★
#
# This file REIMPLEMENTS converse-KAM (C2) and WBA (C3a) from scratch.
# The consistency check (records/RECORD_A_consistency_check.md) found these
# reimplementations AGREE with the originals (converse-KAM detection 100%,
# WBA regular/chaotic classification 93–96%) — so the C3b bug was the magnetic
# field (vpd/field.py), NOT these diagnostics. They are still not a source of
# record.
#
# ★ The source of record is tools/converse_kam.py (C2) and tools/wba.py (C3a).
# ★ The corrected C3b v2 uses those originals via consistency/orig/ (pinned to
#   bf5ca06); see vpd/stage3_paul.py.
# ★ Retained only so that RECORD_C3b_vpd.md (SUPERSEDED) stays reproducible.

"""
Field-line diagnostics on the SAME field and (rho, theta) grid used for V_PD,
so the three axes can be correlated point by point.

  * WBA  (weighted Birkhoff average): dig = number of converged digits of the
    weighted time average of a smooth observable along the Poincare orbit.
    High dig  -> regular (KAM surface OR island: quasiperiodic).
    Low  dig  -> chaotic.  [dynamics]

  * converse-KAM: MacKay-style tangent cone-crossing test.  Integrate the
    field line together with a tangent vector initialised vertical (delta_rho=1,
    delta_theta=0).  Positive shear tips it to delta_theta>0; if the point lies
    on a rotational invariant surface (a graph rho=W(theta,zeta)) the tangent
    can never rotate back through delta_theta=0.  The first zeta at which
    delta_theta<0 proves non-existence: t_c.  Fires inside islands (the original
    rotational surface is destroyed) AND in chaos; silent on surviving KAM
    surfaces.  [topology]

NOTE: these are self-contained reimplementations of the C2 (converse-KAM) and
C3a (WBA) diagnostics -- the original project code is not in this repository.
They are calibrated to reproduce the qualitative behaviour, not the exact
numbers pre-filled in the comparison table.  See RECORD_C3b_vpd.md.
"""
import numpy as np
from .field import IOTA0, IOTAP


def _rhs(rho, theta, zeta, marr, narr, epsarr):
    ph = marr[:, None] * theta[None, :] - narr[:, None] * zeta
    s = np.sin(ph); c = np.cos(ph)
    drho = (epsarr[:, None] * marr[:, None] * s).sum(axis=0)
    dtheta = IOTA0 + IOTAP * rho
    # tangent coefficients
    a_rt = (epsarr[:, None] * marr[:, None] ** 2 * c).sum(axis=0)  # d(drho)/dtheta
    return drho, dtheta, a_rt


def _step(rho, theta, drho_t, dth_t, zeta, dz, marr, narr, epsarr):
    """One RK4 step of state + tangent."""
    def deriv(r, t, dr, dt, z):
        drho, dtheta, a_rt = _rhs(r, t, z, marr, narr, epsarr)
        # tangent: d(dr)/dz = a_rt*dt ; d(dt)/dz = IOTAP*dr
        return drho, dtheta, a_rt * dt, IOTAP * dr

    k1 = deriv(rho, theta, drho_t, dth_t, zeta)
    k2 = deriv(rho + 0.5 * dz * k1[0], theta + 0.5 * dz * k1[1],
               drho_t + 0.5 * dz * k1[2], dth_t + 0.5 * dz * k1[3], zeta + 0.5 * dz)
    k3 = deriv(rho + 0.5 * dz * k2[0], theta + 0.5 * dz * k2[1],
               drho_t + 0.5 * dz * k2[2], dth_t + 0.5 * dz * k2[3], zeta + 0.5 * dz)
    k4 = deriv(rho + dz * k3[0], theta + dz * k3[1],
               drho_t + dz * k3[2], dth_t + dz * k3[3], zeta + dz)
    rho = rho + dz / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    theta = theta + dz / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    drho_t = drho_t + dz / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
    dth_t = dth_t + dz / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
    return rho, theta, drho_t, dth_t


def _bump_weights(n):
    t = (np.arange(n) + 0.5) / n
    g = np.exp(-1.0 / (t * (1.0 - t)))
    return g / g.sum()


def run_diagnostics(field, rho_grid, theta_grid, n_periods=2000, substeps=32,
                    ck_periods=200):
    """Return dict of 2-D arrays over (rho, theta): dig, t_c, detected."""
    marr = np.array([m for m, n, e in field.modes], float)
    narr = np.array([n for m, n, e in field.modes], float)
    epsarr = np.array([e for m, n, e in field.modes], float)

    RR, TT = np.meshgrid(rho_grid, theta_grid, indexing="ij")
    shape = RR.shape
    rho = RR.ravel().copy(); theta = TT.ravel().copy()
    K = rho.size
    dr_t = np.ones(K); dt_t = np.zeros(K)     # vertical tangent
    tc = np.full(K, np.nan)
    detected = np.zeros(K, bool)
    went_pos = np.zeros(K, bool)

    dz = 2 * np.pi / substeps
    total = 2 * n_periods                      # need 2N for two WBA windows
    h1 = np.zeros((total, K)); h2 = np.zeros((total, K))
    zeta = 0.0
    # record period-0 sample
    for p in range(total):
        # stroboscopic observables at start of period p
        h1[p] = np.cos(theta); h2[p] = rho
        for s in range(substeps):
            rho, theta, dr_t, dt_t = _step(rho, theta, dr_t, dt_t, zeta, dz,
                                           marr, narr, epsarr)
            zeta += dz
        # converse-KAM monitoring (first ck_periods only, cheap)
        if p < ck_periods:
            went_pos |= (dt_t > 1e-12)
            flip = went_pos & (dt_t < 0) & (~detected)
            tc[flip] = zeta / (2 * np.pi)
            detected |= flip
        # always renormalise tangent to avoid overflow (sign preserved)
        nrm = np.hypot(dr_t, dt_t) + 1e-300
        dr_t /= nrm; dt_t /= nrm

    # weighted Birkhoff over two consecutive windows of length n_periods
    w = _bump_weights(n_periods)
    def dig_of(h):
        wb1 = (w[:, None] * h[:n_periods]).sum(axis=0)
        wb2 = (w[:, None] * h[n_periods:2 * n_periods]).sum(axis=0)
        err = np.abs(wb1 - wb2)
        return -np.log10(np.clip(err, 1e-16, None))
    dig = np.maximum(dig_of(h1), dig_of(h2))
    dig = np.clip(dig, 0, 16)

    return dict(dig=dig.reshape(shape), t_c=tc.reshape(shape),
                detected=detected.reshape(shape),
                rho=rho_grid, theta=theta_grid)
