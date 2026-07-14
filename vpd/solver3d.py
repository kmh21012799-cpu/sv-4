"""
3-D anisotropic heat-diffusion solver on a structured (rho, theta, zeta) grid.

  div( kappa . grad T ) = 0,   kappa = kperp I + (kpar-kperp) B B / |B|^2
  Dirichlet T=0 at rho_min, T=1 at rho_max; periodic in theta and zeta.

Trilinear (Q1) finite elements => symmetric positive-definite stiffness matrix.
Assembly is fully vectorised over elements (the reference-element gradients are
identical for every cell on the uniform grid, only the tensor varies).

Poloidal-periodicity reduction: when every mode of the field shares the same
poloidal number m, T is periodic in theta with period 2*pi/m, so we solve on
theta in [0, 2*pi/m).  Pass theta_period accordingly.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .field import iota

_g = 1.0 / np.sqrt(3.0)
_QP1 = np.array([0.5 - 0.5 * _g, 0.5 + 0.5 * _g])   # 1-D Gauss nodes on [0,1]


def _ref_grads():
    """Return (gauss offsets (8,3), reference grads gradref (8pts,3,8nodes),
    weights (8,))."""
    # node ordering l = a_rho*4 + a_th*2 + a_ze, a in {0,1}
    nodes = [(ar, at, az) for ar in (0, 1) for at in (0, 1) for az in (0, 1)]
    gps = [(x, y, z) for x in _QP1 for y in _QP1 for z in _QP1]
    gradref = np.zeros((8, 3, 8))
    offs = np.zeros((8, 3))
    for q, (gx, gy, gz) in enumerate(gps):
        offs[q] = (gx, gy, gz)
        for l, (ar, at, az) in enumerate(nodes):
            # shape N_l = fx*fy*fz with fx = gx if ar==1 else 1-gx
            fx = gx if ar else 1 - gx
            fy = gy if at else 1 - gy
            fz = gz if az else 1 - gz
            dfx = 1.0 if ar else -1.0
            dfy = 1.0 if at else -1.0
            dfz = 1.0 if az else -1.0
            gradref[q, 0, l] = dfx * fy * fz
            gradref[q, 1, l] = fx * dfy * fz
            gradref[q, 2, l] = fx * fy * dfz
    w = np.full(8, 0.125)
    return offs, gradref, w, nodes


def _amg_solve(Kff, rhs, x0=None, tol=1e-8):
    """CG preconditioned by smoothed-aggregation AMG (SPD anisotropic)."""
    import pyamg
    ml = pyamg.smoothed_aggregation_solver(Kff.tocsr(), max_coarse=2000)
    M = ml.aspreconditioner(cycle="V")
    res = []
    x = spla.cg(Kff, rhs, rtol=tol, maxiter=2000, M=M, x0=x0,
                callback=lambda xk: res.append(0))
    xf, info = x
    return xf, len(res), info


def solve3d(field, kpar, kperp, rho_min, rho_max, theta_period,
            Nrho, Ntheta, Nzeta, core=(0.25, 0.75), return_fields=False,
            method="amg", x0=None):
    rho = np.linspace(rho_min, rho_max, Nrho)
    hrho = rho[1] - rho[0]
    htheta = theta_period / Ntheta
    hzeta = 2.0 * np.pi / Nzeta
    theta = np.arange(Ntheta) * htheta
    zeta = np.arange(Nzeta) * hzeta

    offs, gradref, gw, _ = _ref_grads()
    # physical-space gradients (constant per element on the uniform grid)
    scale = np.array([1.0 / hrho, 1.0 / htheta, 1.0 / hzeta])
    gradphys = gradref * scale[None, :, None]        # (8pts,3,8nodes)
    detJ = hrho * htheta * hzeta

    Nre, Nte, Nze = Nrho - 1, Ntheta, Nzeta
    Nelem = Nre * Nte * Nze
    ii, jj, kk = np.meshgrid(np.arange(Nre), np.arange(Nte), np.arange(Nze),
                             indexing="ij")
    ii = ii.ravel(); jj = jj.ravel(); kk = kk.ravel()

    def gidx(i, j, k):
        return i * (Ntheta * Nzeta) + (j % Ntheta) * Nzeta + (k % Nzeta)

    # global node index per element per local node (Nelem, 8)
    nodes = [(ar, at, az) for ar in (0, 1) for at in (0, 1) for az in (0, 1)]
    gnodes = np.zeros((Nelem, 8), dtype=np.int64)
    for l, (ar, at, az) in enumerate(nodes):
        gnodes[:, l] = gidx(ii + ar, jj + at, kk + az)

    N = Nrho * Ntheta * Nzeta
    Ke = np.zeros((Nelem, 8, 8))
    rho0 = rho[ii]; th0 = theta[jj]; ze0 = zeta[kk]
    for q in range(8):
        rq = rho0 + offs[q, 0] * hrho
        tq = th0 + offs[q, 1] * htheta
        zq = ze0 + offs[q, 2] * hzeta
        br, bt, bz = field.B(rq, tq, zq)
        B = np.stack([br, bt, bz], axis=1)            # (Nelem,3)
        B2 = (B * B).sum(axis=1)
        outer = B[:, :, None] * B[:, None, :]         # (Nelem,3,3)
        Kq = (kperp * np.eye(3)[None] +
              (kpar - kperp) * outer / B2[:, None, None])
        G = gradphys[q]                               # (3,8)
        # Ke[e,a,b] += w*detJ * G[i,a] Kq[e,i,j] G[j,b]
        GK = np.einsum('ia,eij->eaj', G, Kq)          # (Nelem,8,3)
        Ke += gw[q] * detJ * np.einsum('eaj,jb->eab', GK, G)

    I = np.broadcast_to(gnodes[:, :, None], (Nelem, 8, 8))
    J = np.broadcast_to(gnodes[:, None, :], (Nelem, 8, 8))
    K = sp.coo_matrix((Ke.ravel(), (I.ravel(), J.ravel())), shape=(N, N)).tocsr()

    # Dirichlet
    Tb = np.zeros(N)
    is_bnd = np.zeros(N, dtype=bool)
    for j in range(Ntheta):
        for k in range(Nzeta):
            is_bnd[gidx(0, j, k)] = True
            is_bnd[gidx(Nrho - 1, j, k)] = True
            Tb[gidx(Nrho - 1, j, k)] = 1.0
    free = ~is_bnd
    Kff = K[free][:, free].tocsc()
    rhs = -(K[free][:, is_bnd] @ Tb[is_bnd])
    if method == "direct":
        Tf = spla.spsolve(Kff, rhs)
        n_it = 0
    else:
        x0f = np.asarray(x0).ravel()[free] if x0 is not None else None
        Tf, n_it, info = _amg_solve(Kff, rhs, x0=x0f)
    res = np.linalg.norm(Kff @ Tf - rhs) / (np.linalg.norm(rhs) + 1e-300)
    T = Tb.copy(); T[free] = Tf
    Tg = T.reshape(Nrho, Ntheta, Nzeta)

    interior = Tg[1:-1]
    mp_min, mp_max = float(interior.min()), float(interior.max())
    mp_ok = (mp_min > -1e-6) and (mp_max < 1 + 1e-6)

    diag = _diagnostics(field, Tg, rho, theta, zeta, hrho, htheta, hzeta,
                        kpar, kperp, core)
    diag.update(residual=float(res), max_principle_ok=bool(mp_ok),
                mp_min=mp_min, mp_max=mp_max, Nrho=Nrho, Ntheta=Ntheta,
                Nzeta=Nzeta, n_iter=int(n_it))
    diag["T"] = Tg
    diag["rho"] = rho
    return diag


def _diagnostics(field, Tg, rho, theta, zeta, hrho, htheta, hzeta,
                 kpar, kperp, core):
    Nrho, Ntheta, Nzeta = Tg.shape
    # spectral derivatives in theta, zeta; central FD in rho
    kth = np.fft.fftfreq(Ntheta, d=htheta / (2 * np.pi))
    kze = np.fft.fftfreq(Nzeta, d=hzeta / (2 * np.pi))
    Tth = np.real(np.fft.ifft(1j * kth[None, :, None] *
                              np.fft.fft(Tg, axis=1), axis=1))
    Tze = np.real(np.fft.ifft(1j * kze[None, None, :] *
                              np.fft.fft(Tg, axis=2), axis=2))
    Trho = np.zeros_like(Tg)
    Trho[1:-1] = (Tg[2:] - Tg[:-2]) / (2 * hrho)
    Trho[0] = (Tg[1] - Tg[0]) / hrho
    Trho[-1] = (Tg[-1] - Tg[-2]) / hrho

    R = rho[:, None, None]; TH = theta[None, :, None]; ZE = zeta[None, None, :]
    br, bt, bz = field.B(R, TH, ZE)
    B2 = br * br + bt * bt + bz * bz
    Bdot = br * Trho + bt * Tth + bz * Tze
    bdotgrad = Bdot / np.sqrt(B2)
    grad2 = Trho**2 + Tth**2 + Tze**2
    par2 = bdotgrad**2
    perp2 = grad2 - par2
    chi = (kpar * par2 - kperp * perp2 > 0.0).astype(float)

    cmask = (rho >= core[0]) & (rho <= core[1])
    V_PD = float(chi[cmask].mean())

    def avg_at(rt):
        idx = int(np.argmin(np.abs(rho - rt)))
        return float(Tg[idx].mean())
    DeltaT = avg_at(core[1]) - avg_at(core[0])
    return dict(V_PD=V_PD, DeltaT=DeltaT, chi=chi, rho=rho,
                core_mask=cmask)
