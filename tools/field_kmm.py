"""
KMM magnetic field  --  Kallinikos, MacKay, Martinez-del-Rio, Tigkas, Turnbull (2023)
                        "The converse KAM method ..." arXiv:2304.09613

We implement the *auxiliary* vector field V, with  B = (B0 R0 / R^2) V.
Because V^phi = 1 identically, phi is a legitimate time variable and the
field-line flow reduces to a 1.5-DOF (periodically forced) system

        dpsi/dphi = V^psi(psi, theta, phi)
        dtheta/dphi = V^theta(psi, theta, phi)

which has the *same* invariant tori as B.  This is exactly the simplification
the paper exploits (Sec. 2) and is what we use for both Poincare and
converse-KAM.

Vector potential (covariant):
    A_psi   = 0
    A_theta = psi
    A_phi   = -[ w1 psi + w2 psi^2
                 + sum_mn eps_mn psi^{m/2} f_mn(psi) cos(m theta - n phi + zeta_mn) ]

Auxiliary field (contravariant):
    V^psi   =  sum_mn m eps_mn psi^{m/2} f_mn(psi) sin(m theta - n phi + zeta_mn)
    V^theta =  w1 + 2 w2 psi
               + sum_mn eps_mn psi^{m/2-1} [ (m/2) f_mn + psi f'_mn ] cos(...)
    V^phi   =  1

Standard parameters (paper eq. 13):
    w1 = 1/4,  w2 = 1,  B0 = 1,  R0 = 2,  zeta_mn = 0
    f(psi) = psi - R0^2/B0 = psi - 4          (same f for every mode)

Symplectic (area-preserving) plotting coordinates:
    ytil = sqrt(2 psi / B0) cos theta
    ztil = sqrt(2 psi / B0) sin theta
In these coordinates  d ytil ^ d ztil = d psi ^ d theta, i.e. *area = toroidal flux*.
"""

import numpy as np
from scipy.integrate import solve_ivp


class KMMField:
    """KMM auxiliary vector field V and its field-line flow (phi as time).

    modes : list of (m, n, eps) tuples.  f_mn(psi) = psi - R0^2/B0 for every mode
            (as in every example of the paper), so f'_mn = 1.
    """

    def __init__(self, modes, w1=0.25, w2=1.0, B0=1.0, R0=2.0, zeta=0.0):
        self.modes = [(int(m), int(n), float(eps)) for (m, n, eps) in modes]
        self.w1 = w1
        self.w2 = w2
        self.B0 = B0
        self.R0 = R0
        self.zeta = zeta
        # f(psi) = psi - R0^2/B0 ,  f'(psi) = 1
        self._c = R0 ** 2 / B0

    # ---- building blocks ------------------------------------------------
    def f(self, psi):
        return psi - self._c

    def fp(self, psi):
        return np.ones_like(np.asarray(psi, dtype=float))

    # ---- contravariant components of V ---------------------------------
    def V(self, psi, theta, phi):
        """Return (V^psi, V^theta, V^phi)."""
        psi = np.asarray(psi, dtype=float)
        theta = np.asarray(theta, dtype=float)
        f = self.f(psi)
        fp = self.fp(psi)

        Vpsi = np.zeros_like(psi)
        Vth = self.w1 + 2.0 * self.w2 * psi
        for (m, n, eps) in self.modes:
            ang = m * theta - n * phi + self.zeta
            # psi^{m/2}, guard against tiny negative psi from integrator noise
            ph = np.power(np.clip(psi, 0.0, None), m / 2.0)
            ph1 = np.power(np.clip(psi, 1e-300, None), m / 2.0 - 1.0)
            Vpsi = Vpsi + m * eps * ph * f * np.sin(ang)
            Vth = Vth + eps * ph1 * ((m / 2.0) * f + psi * fp) * np.cos(ang)
        Vphi = np.ones_like(psi)
        return Vpsi, Vth, Vphi

    def A_phi(self, psi, theta, phi):
        """Covariant A_phi (used for invariant / island geometry)."""
        psi = np.asarray(psi, dtype=float)
        val = self.w1 * psi + self.w2 * psi ** 2
        for (m, n, eps) in self.modes:
            ang = m * theta - n * phi + self.zeta
            ph = np.power(np.clip(psi, 0.0, None), m / 2.0)
            val = val + eps * ph * self.f(psi) * np.cos(ang)
        return -val

    # ---- field-line RHS (phi = time) -----------------------------------
    def rhs(self, phi, y):
        """dy/dphi with y = [psi, theta].  Since V^phi = 1, this is V^psi, V^theta."""
        psi, theta = y
        Vpsi, Vth, _ = self.V(psi, theta, phi)
        return [float(Vpsi), float(Vth)]

    # ---- Poincare map (section phi = 0 mod 2pi) ------------------------
    def poincare(self, psi0, theta0, n_cross=400, rtol=1e-10, atol=1e-12):
        """Iterate the field line, sampling at phi = 2 pi k.  Returns psi[], theta[]."""
        psis = np.empty(n_cross + 1)
        thetas = np.empty(n_cross + 1)
        psis[0], thetas[0] = psi0, theta0
        y = [psi0, theta0]
        for k in range(n_cross):
            sol = solve_ivp(self.rhs, (0.0, 2.0 * np.pi), y,
                            rtol=rtol, atol=atol, method="DOP853", dense_output=False)
            y = [sol.y[0, -1], sol.y[1, -1]]
            psis[k + 1] = y[0]
            thetas[k + 1] = y[1]
        return psis, thetas

    # ---- invariant for a *single* (m,n) resonance ----------------------
    def invariant_single(self, psi, theta, phi=0.0):
        """Psi = -n psi - m A_phi  (conserved iff exactly one resonance).

        For example 1 this is an exact field-line invariant; its level sets
        coincide with the Poincare section at phi = 0.
        """
        if len(self.modes) != 1:
            raise ValueError("invariant defined only for a single-resonance field")
        m, n, _ = self.modes[0]
        return -n * np.asarray(psi, float) - m * self.A_phi(psi, theta, phi)

    # ---- symplectic coordinates ----------------------------------------
    def to_symplectic(self, psi, theta):
        r = np.sqrt(2.0 * np.asarray(psi, float) / self.B0)
        return r * np.cos(theta), r * np.sin(theta)

    @staticmethod
    def from_symplectic(y, z, B0=1.0):
        psi = 0.5 * B0 * (y ** 2 + z ** 2)
        theta = np.arctan2(z, y)
        return psi, theta


# ---------------------------------------------------------------------------
# Factory helpers for the three worked examples (paper Sec. 4 / eqs. 14-16)
# ---------------------------------------------------------------------------
def example1(eps=0.004):
    """Integrable: single (2,1) resonance.  A_phi = -[psi/4 + psi^2 + eps psi(psi-4)cos(2th-phi)]."""
    return KMMField([(2, 1, eps)])


def example2(eps=0.003):
    """2/1 + 3/2 resonances (same eps for both)."""
    return KMMField([(2, 1, eps), (3, 2, eps)])


def example3(eps21=0.001, eps54=0.01):
    """2/1 + 5/4 resonances."""
    return KMMField([(2, 1, eps21), (5, 4, eps54)])
