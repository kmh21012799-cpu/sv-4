"""
Magnetic field model for the V_PD (effective volume of parallel diffusion)
computation.

MODEL (reduced "slab"/screw-pinch field-line Hamiltonian)
---------------------------------------------------------
We use flat coordinates (rho, theta, zeta) with a Euclidean metric, theta and
zeta periodic with period 2*pi.  This is the geometry in which Paul, Hudson &
Helander's boundary-layer analysis (their sections 5-6) is derived, so it is
the natural setting in which to reproduce their two analytic validation cases.

The field-line Hamiltonian (poloidal flux function) is

    psi(rho, theta, zeta) = Psi(rho) + sum_{m,n} eps_{m,n} cos(m*theta - n*zeta)

with Psi'(rho) = iota(rho) the rotational transform.  The physical magnetic
field in Cartesian slab components (rho, theta, zeta) is

    B = ( -psi_theta ,  psi_rho ,  1 )
      = (  sum eps m sin(m th - n ze) ,  iota(rho) ,  1 )

which is divergence free (dB = 0).  The field-line flow (zeta as "time") is the
standard pendulum/standard-map system

    d rho / d zeta   = -psi_theta = sum eps m sin(m th - n ze)
    d theta / d zeta =  psi_rho   = iota(rho)

A single resonance (m,n) at iota(rho_*) = n/m produces an island of
half-width (in rho)

    W_half = 2 * sqrt(eps_{m,n} / iota')          (full width 4 sqrt(eps/iota'))

Adjacent same-m resonances n and n+1 sit a distance delta_rho = 1/(m iota')
apart, so the Chirikov overlap parameter is

    s = (W_n + W_{n+1}) / (2 delta_rho) = 4 m sqrt(eps iota')   (equal eps)

Critical overlap s = 1 therefore requires eps_crit(m) = 1 / (16 m^2 iota').

IMPORTANT: none of these amplitudes are tuned to reproduce Paul's numbers.
iota' = 1 and iota_0 = 0.5 are fixed once, and every field amplitude follows
from the Chirikov formula above.  See RECORD_C3b_vpd.md (LIMITATIONS).
"""

import numpy as np

# Fixed, untuned model constants -------------------------------------------
IOTA0 = 0.5   # iota(rho) = IOTA0 + IOTAP * rho
IOTAP = 1.0   # iota' (constant shear)


def iota(rho):
    return IOTA0 + IOTAP * rho


def iota_prime():
    return IOTAP


class Field:
    """A magnetic field = a list of resonant modes (m, n, eps)."""

    def __init__(self, modes, label=""):
        # modes: list of (m, n, eps)
        self.modes = [(int(m), int(n), float(eps)) for (m, n, eps) in modes]
        self.label = label

    # -- flux function and its derivatives (analytic) -----------------------
    def psi_theta(self, rho, theta, zeta):
        out = np.zeros(np.broadcast(rho, theta, zeta).shape)
        for m, n, eps in self.modes:
            out += eps * m * np.sin(m * theta - n * zeta)
        return out  # = -sum eps m sin ... with a sign; see below

    def B(self, rho, theta, zeta):
        """Return Cartesian slab components (B_rho, B_theta, B_zeta)."""
        shape = np.broadcast(rho, theta, zeta).shape
        b_rho = np.zeros(shape)
        for m, n, eps in self.modes:
            b_rho = b_rho + eps * m * np.sin(m * theta - n * zeta)
        b_theta = np.broadcast_to(iota(rho), shape).astype(float)
        b_zeta = np.ones(shape)
        return b_rho, b_theta, b_zeta

    def Bmag2(self, rho, theta, zeta):
        br, bt, bz = self.B(rho, theta, zeta)
        return br * br + bt * bt + bz * bz

    # -- field-line map (for island widths / Poincare, cross-checks) --------
    def field_line_rhs(self, zeta, y):
        rho, theta = y
        drho = 0.0
        for m, n, eps in self.modes:
            drho += eps * m * np.sin(m * theta - n * zeta)
        dtheta = iota(rho)
        return [drho, dtheta]

    # -- resonance bookkeeping ---------------------------------------------
    def island_half_width(self, m, n, eps):
        return 2.0 * np.sqrt(max(eps, 0.0) / iota_prime())

    def resonance_rho(self, m, n):
        # iota(rho) = n/m  ->  rho = (n/m - IOTA0)/IOTAP
        return (n / m - IOTA0) / IOTAP

    def quasilinear_D(self):
        """Analytic quasilinear field-line diffusion coefficient
        D_QL = sum_mn pi (eps_mn m)^2 / (m |iota'|).  Robust predictor of the
        chaotic-layer transport (the numerical map estimate saturates against
        the domain for strong chaos)."""
        return sum(np.pi * (eps * m) ** 2 / (m * abs(iota_prime()))
                   for (m, n, eps) in self.modes)


# --- factory helpers -------------------------------------------------------

def single_island_field(m, n, eps):
    """Validation case 1: one resonant island chain."""
    return Field([(m, n, eps)], label=f"island m={m} n={n} eps={eps:.3e}")


def chaotic_layer_field(m, n_list, chirikov):
    """Validation case 2 / Paul four-field construction.

    All modes share poloidal number m; n ranges over n_list; amplitudes set so
    the Chirikov overlap parameter equals `chirikov`.  s = 4 m sqrt(eps iota')
    => eps = (chirikov / (4 m))^2 / iota'.
    """
    eps = (chirikov / (4.0 * m)) ** 2 / iota_prime()
    modes = [(m, n, eps) for n in n_list]
    return Field(modes, label=f"m={m} |n|={len(n_list)} s={chirikov:g}")


def core_resonances(m, rho_lo=0.25, rho_hi=0.75):
    """Integer n such that the (m,n) resonance iota=n/m lies in [rho_lo,rho_hi]."""
    ilo, ihi = iota(rho_lo), iota(rho_hi)
    n_lo = int(np.ceil(m * ilo))
    n_hi = int(np.floor(m * ihi))
    return list(range(n_lo, n_hi + 1))


def paul_field(m, chirikov=1.0, rho_lo=0.25, rho_hi=0.75):
    """One of Paul's four fields: poloidal number m at critical overlap,
    resonances filling the core [rho_lo, rho_hi]."""
    n_list = core_resonances(m, rho_lo, rho_hi)
    return chaotic_layer_field(m, n_list, chirikov)
