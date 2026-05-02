"""
Quantum metric calculation from the Souza-Wilkens-Martin sum rule.

Derivation
----------
The SWM sum rule (in SI units) reads

    ∫₀^∞ dω  Re[σ_µν(ω)] / ω  =  (π e² / ℏ) · (1/V) · Q_µν                    (1)

where Q_µν is the full localization tensor (units of L²) and V is the system
volume. In a periodic system V = N · V_uc and Q_µν = N · N_bound · g_µν, so

    g_µν  =  (ℏ / π e²) · (1 / n_bound) · I_µν                                 (2)

with the optical-conductivity integral

    I_µν  =  ∫₀^∞ dω  Re[σ_µν(ω)] / ω

To compare across materials we form the dimensionless ratio

    κ_µ  =  n_bound^{-(1/2 - 1/d)} · √g_µµ                                     (3)

where d is the spatial dimension (3 for bulk crystals).

Unit conversions (carried out explicitly below, no magic numbers)
----------------------------------------------------------------
VASP gives ε₂(ω) dimensionless and ω in eV. The library computes

    σ_code(ω) = (ω / 4π) ε₂(ω)        [ω in eV, σ_code in eV]

This is the Gaussian-CGS convention. The SI conductivity is

    σ_SI(ω) = ω · ε₀ · ε₂(ω)          [ω in rad/s, σ_SI in S/m]

Working through the conversion (see derivation below in code):

    I_SI [in SI units]  =  (4π ε₀ e / ℏ) · I_code [dimensionless number]

Plugging into (2):

    g [m²]  =  (ℏ / π e²) · (1/n_bound_SI) · I_SI
            =  (ℏ / π e²) · (1/n_bound_SI) · (4π ε₀ e / ℏ) · I_code
            =  (4 ε₀ / e) · (1 / n_bound_SI) · I_code

All factors of ℏ and one factor of e cancel — only ε₀ and e survive.
Converting volume to Å^-3 input gives the final form used in the code.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Fundamental constants (CODATA, SI units)
# ---------------------------------------------------------------------------
HBAR     = 1.054_571_817e-34   # J·s         (reduced Planck constant)
E_CHARGE = 1.602_176_634e-19   # C           (elementary charge)
M_E      = 9.109_383_7015e-31  # kg          (electron mass)
EPSILON0 = 8.854_187_8128e-12  # F/m         (vacuum permittivity)

# Unit conversions
ANG_TO_M = 1.0e-10             # 1 Å = 1e-10 m
EV_TO_J  = E_CHARGE            # 1 eV = e [in C] joules


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class QuantumMetricResult:
    """Quantum metric (g_µµ in Å²) and dimensionless ratio κ_µ along each axis."""

    # Dimensionful per-electron metric, units Å²
    g_xx: float
    g_yy: Optional[float] = None
    g_zz: Optional[float] = None

    # Dimensionless geometric ratio κ = n_bound^{-(1/2 - 1/d)} · √g
    kappa_xx: float = 0.0
    kappa_yy: Optional[float] = None
    kappa_zz: Optional[float] = None

    # Spatial dimension used in the κ normalization
    dim: int = 3

    # Backwards-compatible aliases — the original pipeline named these sqrtG_over_A_*
    @property
    def sqrtG_over_A_xx(self) -> float:
        return self.kappa_xx

    @property
    def sqrtG_over_A_yy(self) -> Optional[float]:
        return self.kappa_yy

    @property
    def sqrtG_over_A_zz(self) -> Optional[float]:
        return self.kappa_zz


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def compute_quantum_metric(
    I_xx: float,
    bound_electron_density: float,
    *,
    I_yy: Optional[float] = None,
    I_zz: Optional[float] = None,
    dim: int = 3,
) -> QuantumMetricResult:
    """Compute g_µµ [Å²] and dimensionless κ_µ along each direction.

    Parameters
    ----------
    I_xx, I_yy, I_zz : float
        Optical-conductivity integrals  ∫ σ_code(ω)/ω dω  along each direction.
        These are the dimensionless numbers produced by `optics.compute_optical_integrals`,
        with ω in eV and σ_code = (ω/4π)·ε₂.
    bound_electron_density : float
        n_bound in units of 1/Å³.
    dim : int
        Spatial dimension d (default 3 for bulk crystals).
        Used in the normalization exponent (1/2 - 1/d).

    Returns
    -------
    QuantumMetricResult with g_µµ in Å² and dimensionless κ_µ.

    Notes
    -----
    Step-by-step (in the body of the function):

      1. Convert n_bound from Å⁻³ to m⁻³ to work in SI.
      2. Convert I_code (dimensionless, ω-in-eV convention) to I_SI (SI units of σ/ω):
             I_SI = (4π ε₀ e / ℏ) · I_code
         The factor (4π ε₀) restores the Gaussian-to-SI conductivity convention,
         and (e/ℏ) restores the eV-to-rad/s frequency conversion.
      3. Compute g_SI [m²] from Eq. (2) of the module docstring:
             g_SI = (ℏ / π e²) · (1/n_bound_SI) · I_SI
                  = (4 ε₀ / e) · (1/n_bound_SI) · I_code        (after cancellation)
      4. Convert g back to Å²: g_Å² = g_SI · 1e20.
      5. Form κ = n_bound^{-(1/2 - 1/d)} · √g (n_bound in Å⁻³, g in Å²
         → κ has units Å · Å^{-(1/2 - 1/d)·3} = dimensionless when d = 3).
    """

    # 1. Bound density to SI (m⁻³)
    n_bound_SI = bound_electron_density / ANG_TO_M**3   # = n_bound_Å · 1e30

    def _g_in_Ang_squared(I_code: float) -> float:
        """Convert one I_code value to g [Å²] via the explicit derivation."""
        # 2. I_code  →  I_SI
        I_SI = (4.0 * np.pi * EPSILON0 * E_CHARGE / HBAR) * I_code
        # 3. SWM:  g = (ℏ / π e²) · (1/n_bound) · I_SI
        g_SI = (HBAR / (np.pi * E_CHARGE**2)) * I_SI / n_bound_SI
        # (Algebraically equivalent to:  g_SI = (4 ε₀ / e) / n_bound_SI · I_code )
        # 4. m² → Å²
        return float(g_SI / ANG_TO_M**2)

    def _kappa(g_AngSq: float) -> float:
        """Dimensionless κ from g [Å²] and n_bound [Å⁻³]."""
        # κ = n_bound^{-(1/2 - 1/d)} · √g
        exponent = -(0.5 - 1.0 / dim)
        return float(bound_electron_density**exponent * np.sqrt(g_AngSq))

    # xx (always present)
    g_xx = _g_in_Ang_squared(I_xx)
    kxx  = _kappa(g_xx)

    # yy, zz (optional)
    g_yy = _g_in_Ang_squared(I_yy) if I_yy is not None else None
    g_zz = _g_in_Ang_squared(I_zz) if I_zz is not None else None
    kyy  = _kappa(g_yy) if g_yy is not None else None
    kzz  = _kappa(g_zz) if g_zz is not None else None

    return QuantumMetricResult(
        g_xx=g_xx, g_yy=g_yy, g_zz=g_zz,
        kappa_xx=kxx, kappa_yy=kyy, kappa_zz=kzz,
        dim=dim,
    )
