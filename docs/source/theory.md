# Theory

This page explains the physical quantities computed by `quantum-metric`. The library implements the Souza–Wilkens–Martin sum rule directly, with all fundamental constants in SI — no fitted prefactors.

## Pipeline overview

Given a VASP optical calculation, the pipeline is:

1. Parse plasma frequency tensors, sumrule, volume, NELECT, NIONS from `OUTCAR`
2. Read the imaginary dielectric function $\varepsilon_2(\omega)$ from `vasprun.xml`
3. Compute the optical conductivity $\sigma(\omega) = \tfrac{\omega}{4\pi}\,\varepsilon_2(\omega)$
4. Integrate $\sigma(\omega)$ to get various moments
5. Compute bound / itinerant electron counts (f-sum rule)
6. Compute the per-electron metric tensor $g_{\mu\nu}$ and the dimensionless ratio $\kappa$

## The optical conductivity integrals

From $\varepsilon_2(\omega)$ (imaginary part of the dielectric function), we compute

$$
\sigma(\omega) = \frac{\omega}{4\pi}\,\varepsilon_2(\omega)
$$

and the moments

$$
\begin{aligned}
\omega_p^2 &= \tfrac{2}{\pi}\!\int \omega\,\varepsilon_2(\omega)\,d\omega \quad\text{(f-sum rule)} \\
I_{\mu\nu} &= \int_{0^+}^\infty \frac{\sigma_{\mu\nu}(\omega)}{\omega}\,d\omega
\end{aligned}
$$

All integrals are evaluated numerically on the VASP energy grid using the trapezoidal rule. Values for $\omega = 0$ are excluded to avoid the $1/\omega$ divergence in the $I$-integrand.

## Bound vs itinerant electrons (f-sum rule)

VASP's `OUTCAR` reports two plasma-frequency-squared tensors:

- **Intraband** ($\omega^2_{p,\mathrm{intra}}$): Drude-like contribution from electrons at the Fermi surface
- **Interband** ($\omega^2_{p,\mathrm{inter}}$): contribution from optical transitions to higher bands

The library uses the f-sum rule, expressed via hydrogen-atom relations, to split `NELECT` into itinerant and bound parts. VASP provides $X_{\rm vasp} = (\hbar\omega_p)^2$ in eV², and the goal is to compute the number density

$$
n = \frac{\varepsilon_0 m_e}{e^2 \hbar^2}\,X_{\rm vasp}
$$

Importing the standard hydrogen-atom relations $E_0 = e^2/(8\pi\varepsilon_0 a_B) = 13.6$ eV and $a_B = 4\pi\varepsilon_0\hbar^2/(e^2 m_e) = 0.529$ Å gives the clean form

$$
\boxed{\; n = \frac{1}{16\pi}\,\cdot\,\frac{1}{a_B^3}\,\cdot\,\frac{X_{\rm vasp}}{E_0^2}\;}
$$

Numerically $1/(16\pi a_B^3 E_0^2) \approx 7.263\times 10^{-4}$ Å⁻³ eV⁻². Applied to the intraband channel, $N_{\rm itinerant} = n_{\rm intra}\cdot V$, and we obtain

$$
N_{\rm bound} = N_{\rm e} - N_{\rm itinerant}, \qquad n_{\rm bound} = \frac{N_{\rm bound}}{V}
$$

where $N_{\rm e}$ is `NELECT`. A free diagnostic comes from applying the same formula to VASP's reported total sumrule: $V/(16\pi a_B^3 E_0^2) \cdot \text{Sumrule}$ should recover `NELECT` if the f-sum integration is well-converged. The library reports this as `sumrule_check_NELECT`.

## The quantum metric from the SWM sum rule

The Souza–Wilkens–Martin sum rule (in SI units) reads

$$
\int_0^\infty d\omega\,\frac{\mathrm{Re}[\sigma_{\mu\nu}(\omega)]}{\omega} = \frac{\pi e^2}{\hbar}\,\frac{1}{V}\,\mathcal{Q}_{\mu\nu}
$$

where $\mathcal{Q}_{\mu\nu}$ is the full localization tensor (units of $L^2$, independent of dimension) and $V$ is the system volume.

In a periodic system we split $V = N\cdot V_{\rm uc}$ and, in the spirit of equipartition, $\mathcal{Q}_{\mu\nu} = N\cdot N_{\rm bound}\cdot g_{\mu\nu}$. The sum rule becomes

$$
\boxed{\;\frac{1}{n_{\rm bound}}\int_{0^+}^\infty d\omega\,\frac{\mathrm{Re}[\sigma_{\mu\nu}(\omega)]}{\omega} = \frac{\pi e^2}{\hbar}\,g_{\mu\nu}\;}
$$

so the per-electron metric tensor is

$$
g_{\mu\nu} = \frac{\hbar}{\pi e^2 n_{\rm bound}}\int_{0^+}^\infty d\omega\,\frac{\mathrm{Re}[\sigma_{\mu\nu}(\omega)]}{\omega}
$$

with units of length squared, regardless of dimension. The library computes this directly using SI values of $\hbar$, $e$, $\varepsilon_0$, with explicit unit conversions in `metric.py` — there is no hidden prefactor.

## Cross-material comparison: the dimensionless ratio $\kappa$

To compare different materials, $g_{\mu\nu}$ (units of $L^2$) is normalized by the effective inter-particle spacing $n_{\rm bound}^{-1/d}$, where $d$ is the spatial dimension. The dimensionless ratio is

$$
\boxed{\;\kappa_\mu = \frac{1}{n_{\rm bound}^{1/2 - 1/d}}\,\sqrt{g_{\mu\mu}}\;}
$$

For 3D bulk crystals ($d = 3$) the exponent is $1/6$. The library exposes `--dim` (CLI) / `dim=` (Python) so 2D systems can be analysed correctly.

For anisotropic materials, $g_{\mu\mu}$ and $\kappa_\mu$ are computed separately along xx, yy, zz using the corresponding direction-resolved $I_{\mu\mu}$.

## Worked example: Ag fcc

For silver (fcc, $V = 16.39\;\text{Å}^3$, `NELECT` = 11 for Ag with $4d^{10}5s^1$):

| Quantity                             | Value     |
|--------------------------------------|-----------|
| $X_{\rm vasp}^{\rm intra}$           | 101.969 eV² |
| Sumrule                              | 925.581 eV² |
| $I_{xx}$                             | 3.6795    |
| $N_{\rm itinerant}$                  | 1.214     |
| $N_{\rm bound}$                      | 9.786     |
| $n_{\rm bound}$                      | 0.597 Å⁻³ |
| $g_{xx}$                             | **0.1362 Å²** |
| $\sqrt{g_{xx}}$                      | 0.369 Å   |
| $\kappa_{xx}$ ($d=3$)                | 0.402     |

The geometric length $\sqrt{g_{xx}} \approx 0.37$ Å is a physically reasonable scale for the spread of bound-electron Wannier orbitals in a transition metal.

## References

- Souza, Wilkens, and Martin, *Phys. Rev. B* **62**, 1666 (2000) — the foundational SWM sum rule
- VASP optics documentation: <https://www.vasp.at/wiki/index.php/LOPTICS>
- M. Gajdoš et al., *Phys. Rev. B* **73**, 045112 (2006) — VASP linear optical properties
