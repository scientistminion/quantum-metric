# Theory

This page explains the physical quantities computed by `quantum-metric`. You're encouraged to read the original papers referenced below for full derivations.

## Pipeline overview

Given a VASP optical calculation, the pipeline is:

1. Parse plasma frequency tensors, sumrule, volume, NELECT, NIONS from `OUTCAR`
2. Read the imaginary dielectric function $\varepsilon_2(\omega)$ from `vasprun.xml`
3. Compute the optical conductivity $\sigma(\omega) = \tfrac{\omega}{4\pi}\,\varepsilon_2(\omega)$
4. Integrate $\sigma(\omega)$ to get various moments
5. Compute bound / itinerant electron counts
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

## Bound vs itinerant electrons: two methods

VASP's `OUTCAR` reports two plasma-frequency-squared tensors:

- **Intraband** ($\omega^2_{p,\mathrm{intra}}$): Drude-like contribution from electrons at the Fermi surface
- **Interband** ($\omega^2_{p,\mathrm{inter}}$): contribution from optical transitions to higher bands

`quantum-metric` offers two methods for splitting the total electron count $N_{\rm e}$ (= `NELECT` in `OUTCAR`) into itinerant and bound parts.

### Method 1: Kai ratio (`--method kai`, default)

$$
K = \left|\frac{\omega^2_{p,\mathrm{intra}}}{\omega^2_{p,\mathrm{intra}} + \omega^2_{p,\mathrm{inter}}}\right|
$$

$$
N_{\rm itinerant} = K \cdot N_{\rm e},\qquad N_{\rm bound} = (1-K)\,N_{\rm e}
$$

This is simple and transparent — the itinerant fraction is just the Drude weight as a fraction of the total spectral weight.

### Method 2: f-sum rule (`--method fsum`)

Both itinerant and bound counts come directly from their respective plasma frequencies:

$$
N_{\rm itinerant} = \frac{\varepsilon_0 m_e V}{\hbar^2}\,\omega^2_{p,\mathrm{intra}}
$$

$$
N_{\rm bound} = \frac{\varepsilon_0 m_e V}{\hbar^2}\,\omega^2_{p,\mathrm{inter}}
$$

If the f-sum rule is exactly satisfied, $N_{\rm itinerant} + N_{\rm bound} = N_{\rm e}$.
The library reports the residual $N_{\rm e} - (N_{\rm itinerant} + N_{\rm bound})$ as a
diagnostic of sum-rule satisfaction.

This uses the absolute magnitude of the intraband plasma frequency to compute a Drude-based electron count, rather than the relative ratio.

The two methods can give noticeably different answers for the same material, especially when the interband contribution is large. Which one is "right" depends on what you're trying to capture.

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

## References

- VASP optics documentation: <https://www.vasp.at/wiki/index.php/LOPTICS>
- Original VASP paper on linear optical properties: M. Gajdoš et al., *Phys. Rev. B* **73**, 045112 (2006)
- For background on the quantum metric in solids, see e.g. reviews on Berry curvature and band geometry
