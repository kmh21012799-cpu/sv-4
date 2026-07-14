"""
Adapter so the C3b V_PD solver (vpd/solver3d.py, solver2d.py) runs on Paul's
ACTUAL field -- the ORIGINAL PaulField (consistency/orig/field_paul.py, pinned
to bf5ca06), which carries the psi(psi-psibar) envelope of paper eq. 4.1.

The field is used verbatim from the original module; this wrapper only exposes
the Cartesian slab components B = (dpsi/dzeta, dtheta/dzeta, 1) that the solver
needs to build the conductivity tensor.  Field-line dynamics come straight from
PaulField.rhs, so nothing about the field is reimplemented.

Paul field (psi=rho, psibar=1, iota'=1):
    B^rho   = sum eps_mn m rho(rho-1) sin(m th - n ze)      (envelope!)
    B^theta = iota' rho + (2 rho - 1) sum eps_mn cos(m th - n ze)
    B^zeta  = 1
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "consistency", "orig"))
import field_paul as _fp   # ORIGINAL Paul field module (bf5ca06)


class PaulSolveField:
    def __init__(self, paulfield, label=""):
        self.pf = paulfield
        self.modes = paulfield.modes
        self.iota_p = paulfield.iota_p
        self.psibar = paulfield.psibar
        self.label = label or getattr(paulfield, "label", "")

    def B(self, rho, theta, zeta):
        shape = np.broadcast(rho, theta, zeta).shape
        rho_b = np.broadcast_to(rho, shape).astype(float)
        brho = np.zeros(shape)
        bth_pert = np.zeros(shape)
        for (m, n, eps) in self.modes:
            ang = m * theta - n * zeta
            brho = brho + eps * m * rho_b * (rho_b - self.psibar) * np.sin(ang)
            bth_pert = bth_pert + eps * np.cos(ang)
        bth = self.iota_p * rho_b + (2.0 * rho_b - self.psibar) * bth_pert
        bzeta = np.ones(shape)
        return brho, bth, bzeta

    def field_line_rhs(self, zeta, y):
        return self.pf.rhs(zeta, y)

    def island_half_width(self, m, n, eps):
        return self.pf.island_half_width(m, n, eps)

    def resonance_rho(self, m, n):
        return n / (m * self.iota_p)

    def quasilinear_D(self):
        """QL field-line diffusion with the envelope: radial-kick amplitude at
        each resonance a_mn = eps m rho_res(rho_res-1); D = sum pi a^2/(m iota')."""
        D = 0.0
        for (m, n, eps) in self.modes:
            rr = n / (m * self.iota_p)
            a = eps * m * rr * (rr - self.psibar)
            D += np.pi * a * a / (m * abs(self.iota_p))
        return D


# --- constructors (original four fields + parametric) ----------------------
def paul_critical(m):
    f = {4: _fp.paul_m4, 12: _fp.paul_m12, 20: _fp.paul_m20, 36: _fp.paul_m36}[m]()
    return PaulSolveField(f, label=f"Paul m={m} (S=1)")


def paul_chaotic_layer(m, n_list, chirikov):
    """eps set so island half-width = chirikov x (spacing/2): eps = S^2/(16 n(m-n))."""
    modes = [(m, int(n), chirikov ** 2 / (16.0 * n * (m - n))) for n in n_list]
    return PaulSolveField(_fp.PaulField(modes), label=f"Paul m={m} S={chirikov:g}")


def paul_single_island(m, n, eps):
    return PaulSolveField(_fp.PaulField([(m, n, eps)]), label=f"Paul island {m},{n}")
