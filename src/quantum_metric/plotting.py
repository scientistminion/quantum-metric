"""
Plotting utilities for the dielectric function and optical conductivity.

Uses matplotlib only (lazy-imported so non-plotting users don't pay the cost).

Unit conventions
----------------
VASP gives ε₂(ω) dimensionless and ω in eV. The library internally uses

    σ_code(ω) = (ω / 4π) · ε₂(ω)        [ω in eV, σ_code in eV]

which is the Gaussian-CGS expression. For *plotting* we convert to SI units
so the y-axis has a physical meaning users can compare against experiment:

    σ_SI(ω) = ω · ε₀ · ε₂(ω)            [ω in rad/s, σ_SI in S/m]

To go from σ_code (with ω in eV) to σ_SI (with ω in rad/s):

    σ_SI = (4π · ε₀ · e / ℏ) · σ_code   when σ_code is the numerical value
                                         from (ω/4π)·ε₂ with ω in eV.

Equivalently, σ_SI = ε₀ · (e/ℏ) · ω_eV · ε₂  — i.e. just compute it directly.
That's what we do here. Result is in S/m; we plot in MS/m for readability.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from quantum_metric.io import DielectricData

# Physical constants (SI, CODATA)
EPSILON0 = 8.854_187_8128e-12   # F/m
E_CHARGE = 1.602_176_634e-19    # C
HBAR     = 1.054_571_817e-34    # J·s

# eV (energy) → rad/s (angular frequency):  ω_SI = ω_eV · e / ℏ
EV_TO_RAD_PER_S = E_CHARGE / HBAR    # ≈ 1.519e15 rad/s per eV


def _sigma_SI(energy_eV: np.ndarray, eps_imag: np.ndarray) -> np.ndarray:
    """Optical conductivity in SI units (S/m) given ω in eV and ε₂ dimensionless.

    σ_SI(ω) = ω · ε₀ · ε₂(ω),  with ω in rad/s.
    """
    omega_rad_per_s = energy_eV * EV_TO_RAD_PER_S
    return omega_rad_per_s * EPSILON0 * eps_imag   # units: S/m


def plot_dielectric(
    dielectric: DielectricData,
    *,
    output: Optional[str | Path] = None,
    title: Optional[str] = None,
    show: bool = False,
    e_max: Optional[float] = None,
):
    """Plot the imaginary dielectric function eps_2(omega) along available directions."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))

    e = dielectric.energy
    mask = np.ones_like(e, dtype=bool)
    if e_max is not None:
        mask = e <= e_max

    ax.plot(e[mask], dielectric.eps_imag_xx[mask], label=r"$\varepsilon_2^{xx}$")
    if dielectric.has_anisotropic:
        ax.plot(e[mask], dielectric.eps_imag_yy[mask], label=r"$\varepsilon_2^{yy}$")
        ax.plot(e[mask], dielectric.eps_imag_zz[mask], label=r"$\varepsilon_2^{zz}$")

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"$\varepsilon_2(\omega)$  (dimensionless)")
    if title:
        ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=200)
    if show:
        plt.show()
    return fig, ax


def plot_optical_conductivity(
    dielectric: DielectricData,
    *,
    output: Optional[str | Path] = None,
    title: Optional[str] = None,
    show: bool = False,
    e_max: Optional[float] = None,
    units: str = "MS/m",
):
    """Plot σ(ω) along available directions in proper SI units.

    Computes σ(ω) = ω · ε₀ · ε₂(ω) in S/m, then rescales for the chosen y-axis
    units. The default "MS/m" (megasiemens per meter) gives values of order
    ~0.01–10 for typical solids, which is the same scale VASP itself reports
    in the OUTCAR's σ tensor block.

    Parameters
    ----------
    units : str
        One of "S/m", "kS/m", "MS/m" (default).
    """
    import matplotlib.pyplot as plt

    rescale = {"S/m": 1.0, "kS/m": 1.0e-3, "MS/m": 1.0e-6}
    if units not in rescale:
        raise ValueError(f"Unknown units {units!r}. Use one of {list(rescale)}.")
    scale = rescale[units]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    e = dielectric.energy
    mask = e > 0
    if e_max is not None:
        mask &= e <= e_max
    e_m = e[mask]

    ax.plot(e_m, _sigma_SI(e_m, dielectric.eps_imag_xx[mask]) * scale,
            label=r"$\sigma_{xx}$")
    if dielectric.has_anisotropic:
        ax.plot(e_m, _sigma_SI(e_m, dielectric.eps_imag_yy[mask]) * scale,
                label=r"$\sigma_{yy}$")
        ax.plot(e_m, _sigma_SI(e_m, dielectric.eps_imag_zz[mask]) * scale,
                label=r"$\sigma_{zz}$")

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(rf"$\sigma(\omega)$  ({units})")
    if title:
        ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=200)
    if show:
        plt.show()
    return fig, ax
