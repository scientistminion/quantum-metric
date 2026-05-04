# quantum-metric
Compute the **quantum metric** and related optical quantities from [VASP](https://www.vasp.at/) output — in one command.
Given a VASP optical calculation (`LOPTICS = .TRUE.`), this tool extracts:
- Plasma frequencies (intraband, interband) and the f-sum rule, directly from `OUTCAR`
- The frequency-dependent dielectric function ε₂(ω) from `vasprun.xml` (no `sumo` required) or from pre-computed `*_eps_imag.dat`
- Optical conductivity integrals: ∫σ(ω)/ω dω, ∫σ(ω) dω, ∫ωσ(ω) dω, …
- Bound / itinerant electron counts via the f-sum rule with hydrogen-atom relations:

$$n = \frac{1}{16\pi}\,\cdot\,\frac{1}{a_B^3}\,\cdot\,\frac{X_{\rm vasp}}{E_0^2}, \qquad N_{\rm itinerant} = n_{\rm intra}\cdot V, \qquad N_{\rm bound} = N_{\rm e} - N_{\rm itinerant}$$

  with a built-in sumrule consistency check.
- The per-electron quantum metric tensor `g_µµ` (Å²) and dimensionless ratio `κ_µ`, computed directly from the Souza–Wilkens–Martin sum rule:

$$g_{\mu\nu} = \frac{\hbar}{\pi e^2 n_{\rm bound}}\int_{0^+}^\infty\frac{\mathrm{Re}[\sigma_{\mu\nu}(\omega)]}{\omega}d\omega, \qquad \kappa_\mu = \frac{1}{n_{\rm bound}^{1/2 - 1/d}}\sqrt{g_{\mu\mu}}$$

  with all fundamental constants ($\hbar$, $e$, $\varepsilon_0$) in SI — no fitted prefactors.
## Installation
```bash
pip install quantum-metric
```
## Quick start
Point it at a VASP directory and go:
```bash
cd my_vasp_run/
quantum-metric
```
or
```bash
quantum-metric --dir /path/to/my_vasp_run/
```
You'll get a pretty table with all the quantities computed.
### Overriding individual files
```bash
quantum-metric --outcar ../other_run/OUTCAR \
               --poscar ./POSCAR \
               --dielectric ./vasprun.xml
```
### 2D systems
```bash
quantum-metric --dim 2
```
### Different output formats
```bash
quantum-metric --format table    # pretty table (default)
quantum-metric --format json --output results.json
quantum-metric --format tsv  --output results.tsv
quantum-metric --format csv  --output results.csv
```
### Restrict the integration window
```bash
quantum-metric --e-max 15
```
### Plot optical conductivity
```bash
quantum-metric plot --kind optics --output sigma.png
quantum-metric plot --kind epsilon --output eps2.png
```
### Inspect what files will be used
```bash
quantum-metric info
```
## Python API
```python
from quantum_metric import QMetricCalculator
# Auto-discover files in a directory
calc = QMetricCalculator.from_directory("./MoS2")
result = calc.compute()
print(result.metric.g_xx)              # Å²
print(result.metric.kappa_xx)          # dimensionless (d=3)
print(result.electrons.n_bound)
print(result.electrons.sumrule_check)  # diagnostic — should ≈ NELECT
print(result.to_dict())                # flat dict, ready for pandas
```
Or build it from explicit file paths:
```python
calc = QMetricCalculator(
    outcar="./OUTCAR",
    poscar="./POSCAR",
    dielectric="./vasprun.xml",
)
result = calc.compute(dim=3)
```
## Requirements
- Python ≥ 3.9
- `numpy`, `pandas`, `matplotlib`, `typer`, `rich`
- A VASP calculation run with `LOPTICS = .TRUE.` that produced `OUTCAR`, `POSCAR`, and either `vasprun.xml` or an `*_eps_imag.dat` file (e.g. from `sumo-optplot eps_imag --anisotropic`).
## How it works
The computation follows this pipeline:
1. Parse `OUTCAR` for intraband / interband plasma frequency tensors, sumrule, volume, NELECT, NIONS.
2. Parse `POSCAR` for the lattice length `|a|`.
3. Load ε₂(ω) from `vasprun.xml` (preferred) or an `eps_imag.dat` file.
4. Compute σ(ω) = (ω / 4π) ε₂(ω) and the integrals `I_µµ = ∫σ_µµ/ω dω`.
5. Compute the itinerant electron count via the f-sum rule applied to the intraband plasma frequency, then `N_bound = NELECT − N_itinerant`.
6. Compute `g_µµ` (in Å²) directly from the SWM sum rule with all fundamental constants in SI, then `κ_µ = n_bound^{−(1/2 − 1/d)} √g_µµ` for cross-material comparison.
## License
MIT
