"""
2-D anisotropic heat-diffusion solver for the single-helicity (single island
chain) validation case, using the helical symmetry T = T(rho, u), u = m th - n ze.

Solves   div( kappa . grad T ) = 0,   kappa = kperp I + (kpar-kperp) bb/|B|^2
with Dirichlet T=0 at rho_min, T=1 at rho_max, periodic in u.

Bilinear (Q1) finite elements on a structured grid => symmetric positive
definite stiffness matrix (guarantees the discrete maximum principle in the
diagonal-dominant regime and is checked a posteriori).
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .field import iota

# 2x2 Gauss quadrature on the reference element [0,1]^2
_g = 1.0 / np.sqrt(3.0)
_GP = np.array([[0.5 - 0.5 * _g, 0.5 - 0.5 * _g],
                [0.5 + 0.5 * _g, 0.5 - 0.5 * _g],
                [0.5 + 0.5 * _g, 0.5 + 0.5 * _g],
                [0.5 - 0.5 * _g, 0.5 + 0.5 * _g]])
_GW = np.full(4, 0.25)  # weights (area of ref element = 1)


def _shape(xi, eta):
    """Bilinear shape functions N and grads dN/dxi, dN/deta at (xi,eta)."""
    N = np.array([(1 - xi) * (1 - eta), xi * (1 - eta), xi * eta, (1 - xi) * eta])
    dNxi = np.array([-(1 - eta), (1 - eta), eta, -eta])
    dNeta = np.array([-(1 - xi), -xi, xi, (1 - xi)])
    return N, dNxi, dNeta


def _kappa_tensor_3d(field, m, n, rho, u, kpar, kperp):
    """Full 3x3 conductivity tensor at (rho,u) via the field's own B.

    Single-helicity: any (theta,zeta) with m*theta-n*zeta=u gives the same B, so
    evaluate at theta=u/m, zeta=0. Works for any field exposing .B (old constant
    amplitude OR Paul-envelope), so validation runs identically on both."""
    br, bt, bz = field.B(rho, u / m, 0.0)
    B = np.array([float(br), float(bt), float(bz)])
    B2 = B @ B
    K = kperp * np.eye(3) + (kpar - kperp) * np.outer(B, B) / B2
    return K, B, B2


def solve_single_island(field, m, n, kpar, kperp, rho_min, rho_max,
                        Nrho=257, Nu=64, core=(0.25, 0.75)):
    """Assemble and solve; return dict with V_PD, DeltaT, residual, maxprinciple."""
    hu = 2.0 * np.pi / Nu
    rho = np.linspace(rho_min, rho_max, Nrho)
    hrho = rho[1] - rho[0]

    # helical map matrix Mmap: grad3 = Mmap @ (T_rho, T_u)
    Mmap = np.array([[1.0, 0.0], [0.0, m], [0.0, -n]])

    def node(i, j):
        return i * Nu + (j % Nu)

    N = Nrho * Nu
    rows, cols, vals = [], [], []

    # precompute per (element in rho) since tensor depends on rho and u only
    for i in range(Nrho - 1):
        for j in range(Nu):
            # element corners (i,j),(i+1,j),(i+1,j+1),(i,j+1)
            enodes = [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]
            xi_rho = [rho[i], rho[i + 1], rho[i + 1], rho[i]]
            xi_u = [j * hu, j * hu, (j + 1) * hu, (j + 1) * hu]
            Ke = np.zeros((4, 4))
            for q in range(4):
                xi, eta = _GP[q]
                Nsh, dNxi, dNeta = _shape(xi, eta)
                # physical coords at gauss point
                rho_q = sum(Nsh[a] * xi_rho[a] for a in range(4))
                u_q = sum(Nsh[a] * xi_u[a] for a in range(4))
                # Jacobian (axis-aligned rectangle): drho/dxi=hrho, du/deta=hu
                J = np.array([[hrho, 0.0], [0.0, hu]])
                detJ = hrho * hu
                Jinv = np.array([[1.0 / hrho, 0.0], [0.0, 1.0 / hu]])
                # grad of shape wrt (rho,u)
                dN = np.vstack([dNxi, dNeta])          # 2x4 wrt (xi,eta)
                gradN = Jinv.T @ dN                    # 2x4 wrt (rho,u)
                K3, _, _ = _kappa_tensor_3d(field, m, n, rho_q, u_q, kpar, kperp)
                K2 = Mmap.T @ K3 @ Mmap                # effective 2x2 tensor
                # element stiffness contribution
                Ke += _GW[q] * detJ * (gradN.T @ K2 @ gradN)
            for a in range(4):
                for b in range(4):
                    rows.append(enodes[a]); cols.append(enodes[b]); vals.append(Ke[a, b])

    K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))

    # Dirichlet: i=0 -> T=0, i=Nrho-1 -> T=1
    Tb = np.zeros(N)
    is_bnd = np.zeros(N, dtype=bool)
    for j in range(Nu):
        is_bnd[node(0, j)] = True
        Tb[node(0, j)] = 0.0
        is_bnd[node(Nrho - 1, j)] = True
        Tb[node(Nrho - 1, j)] = 1.0
    free = ~is_bnd
    Kff = K[free][:, free]
    Kfb = K[free][:, is_bnd]
    rhs = -Kfb @ Tb[is_bnd]
    Tf = spla.spsolve(Kff.tocsc(), rhs)
    T = Tb.copy()
    T[free] = Tf

    # residual (interior equation residual norm relative to rhs)
    res = np.linalg.norm(Kff @ Tf - rhs) / (np.linalg.norm(rhs) + 1e-300)

    Tgrid = T.reshape(Nrho, Nu)

    # maximum principle: interior extrema must lie within [0,1]
    interior = Tgrid[1:-1, :]
    mp_min, mp_max = interior.min(), interior.max()
    max_principle_ok = (mp_min > -1e-6) and (mp_max < 1.0 + 1e-6)

    # ---- local V_PD indicator and V_PD on the core ----
    # gradients via spectral in u, central/FD in rho
    Tu = _du(Tgrid, hu)
    Trho = _drho(Tgrid, hrho)
    RHO = rho[:, None] * np.ones((1, Nu))
    U = (np.arange(Nu) * hu)[None, :] * np.ones((Nrho, 1))
    b_rho, b_theta, b_zeta = field.B(RHO, U / m, 0.0 * U)
    B2 = b_rho**2 + b_theta**2 + b_zeta**2
    # grad3 T = (Trho, m Tu, -n Tu)
    g_rho, g_th, g_ze = Trho, m * Tu, -n * Tu
    Bdot = b_rho * g_rho + b_theta * g_th + b_zeta * g_ze
    bdotgrad = Bdot / np.sqrt(B2)
    grad2 = g_rho**2 + g_th**2 + g_ze**2
    par2 = bdotgrad**2
    perp2 = grad2 - par2
    chi = (kpar * par2 - kperp * perp2 > 0.0).astype(float)

    core_mask = (rho >= core[0]) & (rho <= core[1])
    # volume element in reduced problem is uniform in (rho,u); ratio cancels
    V_PD = chi[core_mask, :].mean()

    # DeltaT: angle-averaged T at rho=core_hi and rho=core_lo
    def avg_at(rho_target):
        idx = np.argmin(np.abs(rho - rho_target))
        return Tgrid[idx, :].mean()
    DeltaT = avg_at(core[1]) - avg_at(core[0])

    return dict(T=Tgrid, rho=rho, V_PD=float(V_PD), DeltaT=float(DeltaT),
                residual=float(res), max_principle_ok=bool(max_principle_ok),
                mp_min=float(mp_min), mp_max=float(mp_max), chi=chi,
                Nrho=Nrho, Nu=Nu)


def _du(T, hu):
    """Spectral derivative in periodic u (axis=1)."""
    Nu = T.shape[1]
    k = np.fft.fftfreq(Nu, d=hu / (2 * np.pi))  # integer wavenumbers
    Th = np.fft.fft(T, axis=1)
    return np.real(np.fft.ifft(1j * k[None, :] * Th, axis=1))


def _drho(T, hrho):
    """2nd-order central difference in rho (axis=0), one-sided at ends."""
    dT = np.zeros_like(T)
    dT[1:-1, :] = (T[2:, :] - T[:-2, :]) / (2 * hrho)
    dT[0, :] = (T[1, :] - T[0, :]) / hrho
    dT[-1, :] = (T[-1, :] - T[-2, :]) / hrho
    return dT
