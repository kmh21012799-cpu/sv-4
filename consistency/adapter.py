"""
Adapter exposing the C3b (vpd) field through the interface the ORIGINAL
C2/C3a tools (consistency/orig/converse_kam.py, wba.py) expect.

This lets us run the *original* converse-KAM and WBA algorithms on the *same
C3b field* used for the C3b correlations, so the only thing that differs from
the C3b diagnostics is the diagnostic algorithm itself.

C3b field dynamics (from vpd/field.py, vpd/diagnostics.py):
    drho/dzeta   = sum eps_mn m sin(m theta - n zeta)
    dtheta/dzeta = iota(rho) = IOTA0 + IOTAP*rho          (no perturbation term)
Euclidean slab metric (1,1,1) -- identical convention to PaulField.metric.
"""
import numpy as np
from vpd.field import IOTA0, IOTAP


class C3bAdapter:
    def __init__(self, field):
        self._m = np.array([m for (m, n, e) in field.modes], float)
        self._n = np.array([n for (m, n, e) in field.modes], float)
        self._e = np.array([e for (m, n, e) in field.modes], float)
        # attributes some original routines read directly
        self.iota_p = IOTAP
        self.psibar = 1.0
        self.rhobar = 1.0
        self.B0 = 1.0

    # ---- field-line RHS (zeta = time), matches vpd/diagnostics._rhs -------
    def rhs(self, zeta, y):
        rho, theta = y[0], y[1]
        ang = self._m * theta - self._n * zeta
        drho = np.dot(self._e * self._m, np.sin(ang))
        dth = IOTA0 + IOTAP * rho
        return [drho, dth]

    # ---- converse-KAM ingredients ---------------------------------------
    def V(self, psi, theta, zeta):
        drho, dth = self.rhs(zeta, [psi, theta])
        return drho, dth, 1.0

    def jacobian(self, psi, theta, zeta):
        ang = self._m * theta - self._n * zeta
        c = np.cos(ang)
        dVpsi_dpsi = 0.0
        dVpsi_dth = np.dot(self._e * self._m * self._m, c)
        dVth_dpsi = IOTAP
        dVth_dth = 0.0
        return np.array([[dVpsi_dpsi, dVpsi_dth], [dVth_dpsi, dVth_dth]])

    def metric(self, psi):
        return (1.0, 1.0, 1.0)
