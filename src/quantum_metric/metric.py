"""
Quantum metric calculation.

Formula (applied for both 'kai' and 'fsum' electron-counting methods):

    sqrtG = sqrt( prefactor * I / n_bound^(1/3) )

where:
  - prefactor = 0.0694 Angstrom^-1 eV^-1 (unit-conversion constant)
  - I         = integral[ sigma(omega) / omega dw ]   (from optics.py)
  - n_bound   = bound electron density (1/Angstrom^3)

The method (kai vs fsum) only affects how n_bound is computed, not the metric
formula itself. See `electrons.py` for the two electron-counting methods.

Note: historical column name is `sqrtG_over_A_bound`, but no division by a_len
is actually performed — this name is kept for backwards compatibility with
existing TSV outputs from the original pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

# Prefactor from the original pipeline (units: Angstrom^-1 eV^-1).
DEFAULT_PREFACTOR = 0.0694


@dataclass
class QuantumMetricResult:
    """Quantum metric along each direction."""

    sqrtG_xx: float
    sqrtG_yy: Optional[float] = None
    sqrtG_zz: Optional[float] = None
    prefactor: float = DEFAULT_PREFACTOR

    # Backwards-compatible aliases — the old workflow named these 'sqrtG_over_A_*'
    @property
    def sqrtG_over_A_xx(self) -> float:
        return self.sqrtG_xx

    @property
    def sqrtG_over_A_yy(self) -> Optional[float]:
        return self.sqrtG_yy

    @property
    def sqrtG_over_A_zz(self) -> Optional[float]:
        return self.sqrtG_zz


def compute_quantum_metric(
    I_xx: float,
    bound_electron_density: float,
    *,
    I_yy: Optional[float] = None,
    I_zz: Optional[float] = None,
    prefactor: float = DEFAULT_PREFACTOR,
) -> QuantumMetricResult:
    """Compute sqrtG = sqrt(prefactor * I / n^(1/3)) for each direction.

    Parameters
    ----------
    I_xx, I_yy, I_zz : float
        Optical conductivity integrals int(sigma/omega domega) along each direction.
    bound_electron_density : float
        N_bound / V in units of 1/Angstrom^3.
    prefactor : float
        Unit-conversion constant (default 0.0694 A^-1 eV^-1).
    """
    if bound_electron_density == 0:
        raise ValueError(
            "Bound electron density is exactly 0 → N_bound = 0. "
            "Check OUTCAR parsing or try the other --method."
        )
    if bound_electron_density < 0:
        raise ValueError(
            f"Bound electron density is negative ({bound_electron_density:.3g} /Å³). "
            "With --method fsum this means the intraband plasma freq implies more "
            "itinerant electrons than NELECT. Try --method kai."
        )

    den_third = bound_electron_density ** (1.0 / 3.0)

    def _one(I):
        arg = prefactor * I / den_third
        if arg < 0:
            raise ValueError(f"sqrtG would be complex: prefactor*I/n^(1/3) = {arg:.3g} < 0")
        return float(np.sqrt(arg))

    return QuantumMetricResult(
        sqrtG_xx=_one(I_xx),
        sqrtG_yy=_one(I_yy) if I_yy is not None else None,
        sqrtG_zz=_one(I_zz) if I_zz is not None else None,
        prefactor=prefactor,
    )
