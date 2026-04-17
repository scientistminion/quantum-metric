"""
Optical conductivity integrals from the imaginary dielectric function.

Given eps_2(omega) from VASP, we compute:
  sigma(omega)         = (omega / 4 pi) * eps_2(omega)   [VASP / Gaussian units]
  omega_p^2            = (2/pi) * integral[ omega * eps_2(omega) dw ]   (f-sum rule)
  I                    = integral[ sigma(omega) / omega dw ]
  sigma_int            = integral[ sigma(omega) dw ]
  wsigma               = integral[ omega * sigma(omega) dw ]
  sigma_over_wsquare   = integral[ sigma(omega) / omega^2 dw ]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from quantum_metric.io import DielectricData


@dataclass
class OpticalIntegrals:
    """Integrals of the optical conductivity along a single direction.

    All quantities are in the units that follow from eV-based energies and the
    VASP convention sigma = (omega / 4pi) * eps_2. The user is responsible for
    unit-converting if needed downstream.
    """

    omega_p_squared: float   # (2/pi) * int omega * eps_2 domega  (eV^2)
    I: float                 # int sigma/omega domega
    sigma_int: float         # int sigma domega
    wsigma: float            # int omega * sigma domega
    sigma_over_wsquare: float  # int sigma / omega^2 domega


@dataclass
class OpticalIntegralsAllDirections:
    """Optical integrals for xx, yy, zz."""

    xx: OpticalIntegrals
    yy: Optional[OpticalIntegrals] = None
    zz: Optional[OpticalIntegrals] = None


def compute_optical_integrals(
    dielectric: DielectricData,
    *,
    e_min: float = 0.0,
    e_max: Optional[float] = None,
) -> OpticalIntegralsAllDirections:
    """Compute all optical conductivity integrals from eps_2(omega).

    Parameters
    ----------
    dielectric : DielectricData
        Imaginary part of the dielectric function as a function of energy.
    e_min, e_max : float
        Optional integration window in eV. `e_min` is exclusive (we always exclude
        omega == 0 to avoid division-by-zero in the sigma/omega integrand).
    """
    xx = _integrate_one_direction(dielectric.energy, dielectric.eps_imag_xx, e_min, e_max)

    yy = None
    zz = None
    if dielectric.has_anisotropic:
        yy = _integrate_one_direction(dielectric.energy, dielectric.eps_imag_yy, e_min, e_max)
        zz = _integrate_one_direction(dielectric.energy, dielectric.eps_imag_zz, e_min, e_max)

    return OpticalIntegralsAllDirections(xx=xx, yy=yy, zz=zz)


def _integrate_one_direction(
    energy: np.ndarray,
    eps_imag: np.ndarray,
    e_min: float,
    e_max: Optional[float],
) -> OpticalIntegrals:
    mask = energy > max(e_min, 0.0)
    if e_max is not None:
        mask &= energy <= e_max

    if not np.any(mask):
        raise ValueError("No energy points in integration window")

    e = energy[mask]
    eps = eps_imag[mask]

    # sigma(omega) = (omega / 4 pi) * eps_2(omega)
    sigma = (e / (4.0 * np.pi)) * eps

    omega_p_squared = (2.0 / np.pi) * np.trapezoid(e * eps, e)
    I = np.trapezoid(sigma / e, e)
    sigma_int = np.trapezoid(sigma, e)
    wsigma = np.trapezoid(sigma * e, e)
    sigma_over_wsquare = np.trapezoid(sigma / (e**2), e)

    return OpticalIntegrals(
        omega_p_squared=float(omega_p_squared),
        I=float(I),
        sigma_int=float(sigma_int),
        wsigma=float(wsigma),
        sigma_over_wsquare=float(sigma_over_wsquare),
    )
