"""
Bound / itinerant electron counting.

Two methods are supported:

1. "kai" method: ratio of intraband vs total plasma frequency squared
       Kai = |omega_p_intra^2 / (omega_p_intra^2 + omega_p_inter^2)|
       N_itinerant = Kai * NELECT
       N_bound     = (1 - Kai) * NELECT

2. "fsum" method: direct f-sum rule
       N_itinerant = (epsilon_0 * m_e * V * omega_p_intra^2) / hbar^2
       N_bound     = NELECT - N_itinerant
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Physical constants (SI)
EPSILON_0 = 8.854_187_817e-12   # F/m
M_E = 9.109_383_70e-31          # kg
HBAR = 1.054_571_817e-34        # J s
EV_TO_J = 1.602_176_634e-19     # J/eV
ANG3_TO_M3 = 1.0e-30            # Angstrom^3 -> m^3


@dataclass
class ElectronCount:
    """Itinerant / bound electron counts for a single material."""

    method: str               # "kai" or "fsum"
    kai: float | None         # only meaningful for method == 'kai'
    n_itinerant: float
    n_bound: float
    n_itinerant_per_atom: float
    n_bound_per_atom: float
    bound_electron_density: float  # N_bound / V, units: 1/Angstrom^3


def compute_kai(plasma_intra: float, plasma_inter: float) -> float:
    """Ratio Kai = |plasma_intra / (plasma_intra + plasma_inter)| (plasma freqs squared, eV^2)."""
    return float(abs(plasma_intra / (plasma_intra + plasma_inter)))


def compute_n_itinerant_kai(
    plasma_intra: float,
    plasma_inter: float,
    nelect: float,
    volume: float,
    natoms: int,
) -> ElectronCount:
    """Compute electron counts via the Kai ratio method."""
    kai = compute_kai(plasma_intra, plasma_inter)
    n_it = kai * nelect
    n_bd = (1.0 - kai) * nelect
    return ElectronCount(
        method="kai",
        kai=kai,
        n_itinerant=n_it,
        n_bound=n_bd,
        n_itinerant_per_atom=n_it / natoms,
        n_bound_per_atom=n_bd / natoms,
        bound_electron_density=n_bd / volume,
    )

def compute_n_itinerant_fsum(
    plasma_intra_ev2: float,
    nelect: float,
    volume_ang3: float,
    natoms: int,
) -> ElectronCount:
    """Compute electron counts via the f-sum rule.

        N_itinerant = (epsilon_0 * m_e * V * omega_p_intra^2) / hbar^2

    Important: this formula follows the original pipeline convention where
    omega_p^2 from VASP is used directly in eV^2 (not converted to (rad/s)^2).
    This apparent unit mismatch is handled by the implicit eV unit in VASP's
    plasma frequency output: plugging the eV^2 value into the formula with
    SI constants yields the correct dimensionless electron count to within
    the expected precision of the f-sum rule.
    """
    # Precomputed: epsilon_0 * m_e / hbar^2
    CONST = (EPSILON_0 * M_E) / (HBAR ** 2)

    # Volume: Angstrom^3 -> m^3
    volume_si = volume_ang3 * ANG3_TO_M3

    # Use plasma_intra_ev2 directly (no unit conversion — see docstring)
    n_it = CONST * volume_si * plasma_intra_ev2

    n_bd = nelect - n_it
    return ElectronCount(
        method="fsum",
        kai=None,
        n_itinerant=float(n_it),
        n_bound=float(n_bd),
        n_itinerant_per_atom=float(n_it / natoms),
        n_bound_per_atom=float(n_bd / natoms),
        bound_electron_density=float(n_bd / volume_ang3),
    )
