"""
Analytic island geometry for a SINGLE-resonance (integrable) Paul field.
The Paul-coordinate answer key for C2-B: converse-KAM run on a single-resonance
field must reproduce this island area.

Invariant  K = chi - (n/m) psi = (1/2) iota' psi^2 - (n/m) psi
                                  + eps psi(psi-1) cos(m theta - n zeta).
O-point at alpha=0 (elliptic, min of K), X-point at alpha=pi (saddle).
Island area measured in (psi=rho, theta) = toroidal flux (dpsi ^ dtheta).
"""
import numpy as np


def fixed_points(field):
    (m, n, eps) = field.modes[0]
    ip = field.iota_p
    # d K/d psi = 0 with sin(alpha)=0:
    #   alpha=0 : psi_O = (n/m + eps) / (iota' + 2 eps)
    #   alpha=pi: psi_X = (n/m - eps) / (iota' - 2 eps)
    psiO = (n / m + eps) / (ip + 2 * eps)
    psiX = (n / m - eps) / (ip - 2 * eps)
    return psiO, psiX


def island_half_width_pendulum(field):
    """W = 2 sqrt(|eps_tilde|/iota'), eps_tilde = eps psi_s(psi_s-1) (eq. 5.2)."""
    (m, n, eps) = field.modes[0]
    return field.island_half_width(m, n, eps)


def island_area_pendulum(field):
    """Leading pendulum estimate: total (m-island) chain area = 8 W in (psi,theta)."""
    return 8.0 * island_half_width_pendulum(field)


def island_area_gridcount(field, n_theta=1600, n_psi=4000):
    """Exact answer key: cell-count of {K < K_X} over (psi,theta) in [0,1]x[0,2pi]."""
    (m, n, eps) = field.modes[0]
    psiO, psiX = fixed_points(field)
    KX = field.invariant_single(psiX, np.pi / m, 0.0)   # X at alpha=pi -> m*theta=pi
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    ps = np.linspace(1e-4, 1.0, n_psi)
    TH, PS = np.meshgrid(th, ps)
    inside = field.invariant_single(PS, TH, 0.0) < KX
    cell = (2.0 * np.pi / n_theta) * ((1.0 - 1e-4) / (n_psi - 1))
    return float(inside.sum()) * cell
