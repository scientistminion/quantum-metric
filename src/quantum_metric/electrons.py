"""
Electron counting via the f-sum rule.

Number density n is obtained from VASP's plasma frequency squared X_vasp (eV²)
using the hydrogen-atom convenient form:

    n = (1 / (16π)) × (1 / a_B³) × (X_vasp / E_0²)

where a_B = 0.529 Å and E_0 = 13.6 eV (Rydberg).

Applied to the intraband channel:    N_itinerant = n_intra × V
Total bound count:                   N_bound      = NELECT - N_itinerant

This replaces the older "Kai" (ratio-of-plasma-frequencies) method, which
implicitly assumed equal intraband/interband effective masses and gave
unphysical results for systems with semicore states (e.g. Na: Kai claimed
~5.4 itinerant electrons; f-sum correctly gives ~1).
"""

from __future__ import annotations

from dataclasses import dataclass

# Hydrogen-atom constants
A_B = 0.529          # Bohr radius in Å
E_0 = 13.6           # Rydberg energy in eV

# Universal prefactor: n [Å⁻³] = PREFACTOR × X_vasp [eV²]
import math
PREFACTOR_N = 1.0 / (16.0 * math.pi * A_B**3 * E_0**2)
# Numerically ≈ 7.263e-4 Å⁻³ eV⁻²


@dataclass
class ElectronCount:
    """Itinerant / bound electron counts for a single material."""

    n_itinerant: float                # total itinerant electrons in cell
    n_bound: float                    # NELECT - n_itinerant
    n_itinerant_per_atom: float
    n_bound_per_atom: float
    bound_electron_density: float     # n_bound / V, units 1/Å³
    itinerant_electron_density: float # n_itinerant / V, units 1/Å³
    sumrule_check: float              # NELECT_implied from sumrule (diagnostic)


def compute_electron_count(
    plasma_intra_ev2: float,
    nelect: float,
    volume_ang3: float,
    natoms: int,
    sumrule_ev2: float | None = None,
) -> ElectronCount:
    """Compute itinerant/bound electron counts from the intraband f-sum rule.

    Parameters
    ----------
    plasma_intra_ev2 : float
        VASP intraband plasma frequency squared (eV²), diagonal-averaged.
    nelect : float
        Total valence electrons (NELECT from OUTCAR).
    volume_ang3 : float
        Primitive cell volume in Å³.
    natoms : int
        Number of atoms in the cell (NIONS).
    sumrule_ev2 : float, optional
        VASP-reported sumrule (eV²). If provided, used as a consistency check:
        n_total = PREFACTOR_N × sumrule × V should equal NELECT.

    Returns
    -------
    ElectronCount
    """
    # n_itinerant from f-sum rule
    n_it_density = PREFACTOR_N * plasma_intra_ev2          # Å⁻³
    n_it = n_it_density * volume_ang3                      # dimensionless count
    n_bd = nelect - n_it
    n_bd_density = n_bd / volume_ang3

    # Diagnostic: what NELECT does the sumrule imply?
    sumrule_check = (
        PREFACTOR_N * sumrule_ev2 * volume_ang3
        if sumrule_ev2 is not None
        else float("nan")
    )

    return ElectronCount(
        n_itinerant=float(n_it),
        n_bound=float(n_bd),
        n_itinerant_per_atom=float(n_it / natoms),
        n_bound_per_atom=float(n_bd / natoms),
        bound_electron_density=float(n_bd_density),
        itinerant_electron_density=float(n_it_density),
        sumrule_check=float(sumrule_check),
    )
