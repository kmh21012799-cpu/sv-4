"""
Analytic / semi-analytic island geometry for the *integrable* KMM example 1
(single (2,1) resonance).  This provides the "answer key" island area S_I that
the converse-KAM non-existence region must converge to in C1 (paper Fig. 4).

For a single (m,n) resonance the field-line flow has the exact invariant
    Psi(psi, theta, phi) = -n psi - m A_phi(psi, theta, phi)
which, written with the resonant angle  alpha = m theta - n phi , is an
autonomous 1-DOF Hamiltonian in (psi, alpha).  Its separatrix bounds the
island; the enclosed area in (psi, theta) equals the toroidal-flux island area
because  d ytil ^ d ztil = d psi ^ d theta  in symplectic coordinates.
"""

import numpy as np
from scipy.optimize import brentq


def _fixed_points_21(field):
    """O- and X-point psi values for a single (2,1) KMM resonance (eps>0).

    O-point (elliptic) at alpha = 0, X-point (hyperbolic) at alpha = pi.
    Derived from d(alpha)/dphi = 0 with sin(alpha)=0 (see RECORD_C0).
    """
    (m, n, eps) = field.modes[0]
    assert (m, n) == (2, 1), "closed form here is for the (2,1) resonance"
    w1, w2 = field.w1, field.w2
    c = field._c  # = R0^2/B0 = 4
    # d(alpha)/dphi = 2 V^theta - 1
    #  = (2 w1 - 1) + 4 w2 psi + 2 eps (2 psi - c) cos(alpha) = 0
    # alpha = 0  (cos=+1):  psi_O = (1 - 2 w1 + 2 eps c) / (4 w2 + 4 eps)
    # alpha = pi (cos=-1):  psi_X = (1 - 2 w1 - 2 eps c) / (4 w2 - 4 eps)
    psiO = (1.0 - 2.0 * w1 + 2.0 * eps * c) / (4.0 * w2 + 4.0 * eps)
    psiX = (1.0 - 2.0 * w1 - 2.0 * eps * c) / (4.0 * w2 - 4.0 * eps)
    return psiO, psiX


def island_width_analytic(field, theta):
    """Island width Delta psi(theta) from the separatrix level set.

    Returns Delta psi = psi_high - psi_low at toroidal angle theta (phi=0),
    zero outside the island's angular extent.  At each theta we locate the local
    minimum of Psi(.,theta) (the island interior, since the O-line is a min of
    Psi) and bracket the two separatrix roots around it.
    """
    (m, n, eps) = field.modes[0]
    psiO, psiX = _fixed_points_21(field)

    def Psi(psi, th):
        return field.invariant_single(psi, th, phi=0.0)

    PsiX = Psi(psiX, np.pi / 2.0)  # X-point sits at alpha=pi -> 2*theta=pi -> theta=pi/2

    theta = np.atleast_1d(np.asarray(theta, float))
    width = np.zeros_like(theta)
    lo = max(1e-6, psiO - 6.0 * abs(psiO - psiX))
    hi = psiO + 6.0 * abs(psiO - psiX)
    scan = np.linspace(lo, hi, 800)
    for i, th in enumerate(theta):
        vals = Psi(scan, th)
        j = int(np.argmin(vals))
        if vals[j] >= PsiX:
            continue  # no island at this angle
        pcen = scan[j]
        g = lambda p: Psi(p, th) - PsiX
        try:
            pl = brentq(g, lo, pcen, xtol=1e-13, rtol=1e-13)
        except ValueError:
            pl = lo
        try:
            ph = brentq(g, pcen, hi, xtol=1e-13, rtol=1e-13)
        except ValueError:
            ph = hi
        width[i] = ph - pl
    return width


def island_area(field, n_theta=2000):
    """Total island toroidal-flux area  S_I = int_0^{2pi} Delta psi(theta) d theta.

    Covers both islands of the (2,1) chain.  Units: toroidal flux (= area in the
    symplectic (ytil, ztil) plane).
    """
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    w = island_width_analytic(field, th)
    dtheta = 2.0 * np.pi / n_theta
    return float(np.sum(w) * dtheta)


def island_area_gridcheck(field, n_theta=1200, n_psi=4000):
    """Independent cross-check of S_I by direct area counting of {Psi < Psi_X}."""
    (m, n, eps) = field.modes[0]
    psiO, psiX = _fixed_points_21(field)
    PsiX = field.invariant_single(psiX, np.pi / 2.0, 0.0)
    lo = max(1e-6, psiO - 6.0 * abs(psiO - psiX))
    hi = psiO + 6.0 * abs(psiO - psiX)
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    ps = np.linspace(lo, hi, n_psi)
    TH, PS = np.meshgrid(th, ps)
    Z = field.invariant_single(PS, TH, 0.0)
    inside = Z < PsiX
    cell = (2.0 * np.pi / n_theta) * ((hi - lo) / (n_psi - 1))
    return float(np.sum(inside) * cell)
