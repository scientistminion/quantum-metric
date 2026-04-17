"""
Plotting helpers for the optical conductivity and dielectric function.
Uses matplotlib only (lazy-imported so non-plotting users don't pay the cost).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from quantum_metric.io import DielectricData


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
    ax.set_ylabel(r"$\varepsilon_2(\omega)$")
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
):
    """Plot sigma(omega) = (omega / 4pi) * eps_2(omega) along available directions."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))

    e = dielectric.energy
    mask = e > 0
    if e_max is not None:
        mask &= e <= e_max

    def _sigma(eps_imag):
        return (e[mask] / (4.0 * np.pi)) * eps_imag[mask]

    ax.plot(e[mask], _sigma(dielectric.eps_imag_xx), label=r"$\sigma_{xx}$")
    if dielectric.has_anisotropic:
        ax.plot(e[mask], _sigma(dielectric.eps_imag_yy), label=r"$\sigma_{yy}$")
        ax.plot(e[mask], _sigma(dielectric.eps_imag_zz), label=r"$\sigma_{zz}$")

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"$\sigma(\omega)$  (arb. units)")
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
