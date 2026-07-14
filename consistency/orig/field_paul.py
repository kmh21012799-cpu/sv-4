"""
Paul model magnetic field  --  Paul, Hudson, Helander (2022)
                               "Heat conduction in an irregular magnetic field.
                                Part 2 ..." J. Plasma Phys. / arXiv:2108.06328

Field-line Hamiltonian (paper eq. 4.1), with chi playing the role of the
Hamiltonian, psi the momentum, theta the coordinate, zeta the time:

    chi(psi, theta, zeta) = (iota'/2) psi^2
                            + sum_mn eps_mn psi (psi - psibar) cos(m theta - n zeta)

    B = grad psi x grad theta - grad chi x grad zeta                      (eq. 4.2)

Field-line ODEs (zeta as time):
    dpsi/dzeta   = -d chi / d theta =  sum_mn eps_mn m psi(psi-psibar) sin(m th - n z)
    dtheta/dzeta =  d chi / d psi   =  iota' psi
                                       + sum_mn eps_mn (2 psi - psibar) cos(m th - n z)

Coordinates:  psi = psibar rho / rhobar,  rho in [0, rhobar],  theta,zeta in [0,2pi).
Normalisation (paper Sec. 7):  psibar = rhobar = L_zeta = iota' = 1   => psi = rho.

Key property: the perturbation carries the factor psi(psi - psibar), which
*vanishes at both boundaries* rho=0 and rho=rhobar -- no parallel flux leaves
the domain.

Resonances:  iota(psi) = iota' psi = n/m  =>  psi = n/(m iota').
Island half-width (paper eq. 5.2):
    W_mn = (2 rhobar/psibar) sqrt( eps_tilde_mn(psi) / iota'(psi) )   at psi = n/(m iota'),
    with eps_tilde_mn(psi) = eps_mn psi (psi - psibar).
(|.| taken inside the sqrt; the sign of eps_mn only sets the island phase.)
"""

import numpy as np
from scipy.integrate import solve_ivp


class PaulField:
    """Paul (2022) model field-line flow with zeta as time."""

    def __init__(self, modes, iota_p=1.0, psibar=1.0, rhobar=1.0):
        self.modes = [(int(m), int(n), float(eps)) for (m, n, eps) in modes]
        self.iota_p = iota_p
        self.psibar = psibar
        self.rhobar = rhobar
        # cached arrays for a vectorised RHS
        self._m = np.array([m for (m, n, e) in self.modes], float)
        self._n = np.array([n for (m, n, e) in self.modes], float)
        self._e = np.array([e for (m, n, e) in self.modes], float)

    # ---- Hamiltonian ---------------------------------------------------
    def chi(self, psi, theta, zeta):
        psi = np.asarray(psi, float)
        val = 0.5 * self.iota_p * psi ** 2
        for (m, n, eps) in self.modes:
            val = val + eps * psi * (psi - self.psibar) * np.cos(m * theta - n * zeta)
        return val

    # ---- field-line RHS (zeta = time) ----------------------------------
    def rhs(self, zeta, y):
        psi, theta = y
        ang = self._m * theta - self._n * zeta
        s = np.sin(ang)
        c = np.cos(ang)
        pfac = psi * (psi - self.psibar)
        dpsi = np.dot(self._e * self._m, s) * pfac
        dth = self.iota_p * psi + (2.0 * psi - self.psibar) * np.dot(self._e, c)
        return [dpsi, dth]

    def iota(self, psi):
        return self.iota_p * np.asarray(psi, float)

    # ---- converse-KAM ingredients (zeta = time) ------------------------
    # State (psi=rho, theta); the "time" component V^zeta = 1, so a poloidal
    # tangent stays poloidal -- exactly as for the KMM V-flow.
    def V(self, psi, theta, zeta):
        """(V^psi, V^theta, V^zeta) = (dpsi/dzeta, dtheta/dzeta, 1)."""
        dpsi, dth = self.rhs(zeta, [psi, theta])
        return dpsi, dth, 1.0

    def jacobian(self, psi, theta, zeta):
        """Analytic d(dpsi/dzeta, dtheta/dzeta)/d(psi, theta) at a scalar point."""
        psi = float(psi)
        ang = self._m * theta - self._n * zeta
        s = np.sin(ang); c = np.cos(ang)
        pb = self.psibar
        # V^psi = sum e m psi(psi-pb) sin ;  d/dpsi psi(psi-pb) = 2psi-pb
        dVpsi_dpsi = np.dot(self._e * self._m, s) * (2.0 * psi - pb)
        dVpsi_dth = np.dot(self._e * self._m * self._m, c) * (psi * (psi - pb))
        # V^theta = iota' psi + sum e (2psi-pb) cos
        dVth_dpsi = self.iota_p + np.dot(self._e, c) * 2.0
        dVth_dth = -np.dot(self._e * self._m, s) * (2.0 * psi - pb)
        return np.array([[dVpsi_dpsi, dVpsi_dth], [dVth_dpsi, dVth_dth]])

    def metric(self, psi):
        """Euclidean adapted metric: (rho,theta,zeta) orthogonal with |grad rho|=1
        (paper: rho_hat x theta_hat . zeta_hat = 1).  grad psi || d/dpsi holds."""
        return (1.0, 1.0, 1.0)

    def invariant_single(self, psi, theta, zeta=0.0):
        """K = chi - (n/m) psi : exact field-line invariant for ONE resonance.

        Derived in RECORD_C2; conserved because dK/dzeta = 0 for a single mode.
        Its level sets give the analytic island (the Paul-coordinate answer key).
        """
        if len(self.modes) != 1:
            raise ValueError("single-resonance invariant only")
        m, n, eps = self.modes[0]
        psi = np.asarray(psi, float)
        return (0.5 * self.iota_p * psi ** 2 - (n / m) * psi
                + eps * psi * (psi - self.psibar) * np.cos(m * theta - n * zeta))

    # ---- island half-width (eq. 5.2) -----------------------------------
    def island_half_width(self, m, n, eps):
        psi_res = n / (m * self.iota_p)
        eps_tilde = eps * psi_res * (psi_res - self.psibar)
        return (2.0 * self.rhobar / self.psibar) * np.sqrt(abs(eps_tilde) / self.iota_p)

    # ---- Poincare map (section zeta = 0 mod 2pi) -----------------------
    def poincare(self, psi0, theta0, n_cross=400, rtol=1e-9, atol=1e-11):
        """Single long integration sampled at zeta = 2 pi k (fast).

        For chaotic orbits we only need a faithful sampling of the invariant set,
        not a bit-exact trajectory, so one integration with t_eval beats
        restarting the solver every crossing.
        """
        t_eval = 2.0 * np.pi * np.arange(n_cross + 1)
        sol = solve_ivp(self.rhs, (0.0, t_eval[-1]), [psi0, theta0],
                        rtol=rtol, atol=atol, method="DOP853", t_eval=t_eval)
        return sol.y[0], sol.y[1]


# ---------------------------------------------------------------------------
# Critical-overlap chaotic fields (paper Sec. 7)
# ---------------------------------------------------------------------------
def critical_overlap_amplitude(m, n):
    """eps_mn making the island half-width equal to half the resonance spacing.

    Spacing between neighbouring n/m resonances is Delta psi = 1/m.  Setting the
    half-width W = Delta psi/2 = 1/(2m):
        W = 2 sqrt(|eps| |psi(psi-1)|),  psi = n/m,  |psi(psi-1)| = n(m-n)/m^2
        => |eps| = 1 / (16 n (m-n)).
    (psibar = rhobar = iota' = 1.)
    """
    return 1.0 / (16.0 * n * (m - n))


def critical_field(m, ns):
    """Build a PaulField with a chain of n/m resonances at critical overlap.

    ns : iterable of toroidal mode numbers n (so resonances at psi = n/m).
    """
    modes = [(m, int(n), critical_overlap_amplitude(m, int(n))) for n in ns]
    return PaulField(modes)


# The four fields of paper Sec. 7 (Fig. 6).
def paul_m4():
    return critical_field(4, [1, 2, 3])


def paul_m12():
    return critical_field(12, range(2, 11))       # n = 2..10


def paul_m20():
    return critical_field(20, range(3, 18))       # n = 3..17


def paul_m36():
    return critical_field(36, range(5, 32))       # n = 5..31
