# Theory

This page explains the physical quantities computed by `quantum-metric`. You're encouraged to read the original papers referenced below for full derivations.

## Pipeline overview

Given a VASP optical calculation, the pipeline is:

1. Parse plasma frequency tensors, sumrule, volume, NELECT, NIONS from `OUTCAR`
2. Read the imaginary dielectric function $\varepsilon_2(\omega)$ from `vasprun.xml`
3. Compute the optical conductivity $\sigma(\omega) = \tfrac{\omega}{4\pi}\,\varepsilon_2(\omega)$
4. Integrate $\sigma(\omega)$ to get various moments
5. Compute bound / itinerant electron counts (f-sum rule)
6. Compute the quantum metric $\sqrt{G}$

## The optical conductivity integrals

From $\varepsilon_2(\omega)$ (imaginary part of the dielectric function), we compute

$$
\sigma(\omega) = \frac{\omega}{4\pi}\,\varepsilon_2(\omega)
$$

and then the moments

$$
\begin{aligned}
\omega_p^2 &= \tfrac{2}{\pi}\!\int \omega\,\varepsilon_2(\omega)\,d\omega \quad\text{(f-sum rule)} \\
I &= \int \sigma(\omega)/\omega\,d\omega \\
\int \sigma(\omega)\,d\omega,\quad &\int \omega\sigma(\omega)\,d\omega,\quad \int \sigma(\omega)/\omega^2\,d\omega
\end{aligned}
$$

All integrals are evaluated numerically on the VASP energy grid using the trapezoidal rule. Values for $\omega=0$ are excluded to avoid the $1/\omega$ divergence in the $I$-integrand.

The quantity $I = \int \sigma(\omega)/\omega\,d\omega$ is proportional to the **trace of the quantum metric tensor** integrated over the Brillouin zone — this is the connection to the geometric quantity we care about.

## Bound vs itinerant electrons: the f-sum rule

VASP's `OUTCAR` reports two plasma-frequency-squared tensors:

- **Intraband** ($\omega^2_{p,\mathrm{intra}}$): Drude-like contribution from electrons at the Fermi surface
- **Interband** ($\omega^2_{p,\mathrm{inter}}$): contribution from optical transitions to higher bands

It just so happens that, in perturbation theory, organizing factors of $n e^2/(m_e \varepsilon_0)$ in terms of an energy scale has several advantages. One can therefore use the sum rule to define "plasma frequencies", but it should be kept in mind that this is merely a relabeling. We note this because VASP uses this convention.

### Conversion factor via hydrogen-atom relations

VASP provides the square of the plasma frequency in $\text{eV}^2$, which we denote $X_{\rm vasp}$. The goal is to use the relation

$$
n = \frac{\varepsilon_0 m_e}{e^2}\,\omega_p^2
  = \frac{\varepsilon_0 m_e}{e^2 \hbar^2}\,(\hbar\omega_p)^2
  = \left(\frac{\varepsilon_0 m_e}{e^2 \hbar^2}\right) X_{\rm vasp}
$$

to find the number density $n$. While one can use standard values of the universal constants to simplify the prefactors, it is easier to import the standard hydrogen-atom relations

$$
E_0 = \frac{e^2}{8\pi\varepsilon_0 a_B} = 13.6\;\text{eV},\qquad
a_B = \frac{4\pi\varepsilon_0\hbar^2}{e^2 m_e} = 0.529\;\text{Å}
$$

to write the expression as

$$
\boxed{\; n = \frac{1}{16\pi}\,\cdot\,\frac{1}{a_B^3}\,\cdot\,\frac{X_{\rm vasp}}{E_0^2}\;}
$$

Numerically, the prefactor evaluates to

$$
\frac{1}{16\pi\,a_B^3\,E_0^2} \approx 7.263\times 10^{-4}\;\text{Å}^{-3}\,\text{eV}^{-2}
$$

### Applying the formula

The library uses the **intraband** channel for the itinerant count, and obtains the bound count by subtraction:

$$
N_{\rm itinerant} = n_{\rm intra}\cdot V,\qquad n_{\rm intra} = \frac{1}{16\pi\,a_B^3\,E_0^2}\,X_{\rm vasp}^{\rm intra}
$$

$$
N_{\rm bound} = N_{\rm e} - N_{\rm itinerant}
$$

where $N_{\rm e}$ is `NELECT` from `OUTCAR`.

### Sumrule consistency check

As a free diagnostic, the same prefactor applied to VASP's reported total sumrule should recover `NELECT`:

$$
N_{\rm e}^{\rm implied} = \frac{V}{16\pi\,a_B^3\,E_0^2}\cdot\text{Sumrule}
$$

Deviation from the true `NELECT` flags numerical issues such as truncation of the optical integral at finite energy. The library reports `sumrule_check_NELECT` for every calculation.

## The bound electron density

$$
n_{\rm bound} = \frac{N_{\rm bound}}{V}\qquad\text{(units: Å}^{-3}\text{)}
$$

where $V$ is the cell volume.

## The quantum metric

The library computes

$$
\sqrt{G} = \sqrt{\;\alpha\,\frac{I}{n_{\rm bound}^{1/3}}\;}
$$

with the default prefactor $\alpha = 0.0694\;\text{Å}^{-1}\,\text{eV}^{-1}$ (a unit-conversion constant; overridable with `--prefactor`).

For anisotropic materials, $\sqrt{G}$ is computed separately along xx, yy, zz using the corresponding direction-resolved $I_{ii}$.

:::{note}
In the original workflow that this library replaces, this quantity appeared in column names as `sqrtG_over_A_bound` — a naming accident. No division by the lattice length $|a|$ is actually performed. The column names are kept for backwards compatibility.
:::

## Worked example: Na bcc

For sodium (bcc, $V = 39.28\;\text{Å}^3$, `NELECT` $= 7$ for the Na_pv pseudopotential including the $2p^6$ semicore):

| Quantity                             | Value     |
|--------------------------------------|-----------|
| $X_{\rm vasp}^{\rm intra}$           | 37.394 eV² |
| Sumrule                              | 245.701 eV² |
| $n_{\rm intra}$                      | 0.0272 Å⁻³ |
| $N_{\rm itinerant}$                  | **1.067** |
| $N_{\rm bound}$                      | 5.933     |
| $n_{\rm bound}$                      | 0.151 Å⁻³ |
| $N_{\rm e}^{\rm implied}$ (sumrule check) | 7.011 ✓ |

The itinerant count comes out to ~1 per atom, as expected for the lone $3s^1$ conduction electron. The sumrule check lands almost exactly on `NELECT`, validating both the prefactor and the f-sum integration.

## References

- VASP optics documentation: <https://www.vasp.at/wiki/index.php/LOPTICS>
- Original VASP paper on linear optical properties: M. Gajdoš et al., *Phys. Rev. B* **73**, 045112 (2006)
- For background on the quantum metric in solids, see e.g. reviews on Berry curvature and band geometry
