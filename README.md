# quantum-metric

Compute the **quantum metric** and related optical quantities from [VASP](https://www.vasp.at/) output — in one command.

Given a VASP optical calculation (`LOPTICS = .TRUE.`), this tool extracts:

- Plasma frequencies (intraband, interband) and the f-sum rule, directly from `OUTCAR`
- The frequency-dependent dielectric function ε₂(ω) from `vasprun.xml` (no `sumo` required) or from pre-computed `*_eps_imag.dat`
- Optical conductivity integrals: ∫σ(ω)/ω dω, ∫σ(ω) dω, ∫ωσ(ω) dω, ...
- Bound / itinerant electron counts via two methods:
  - **Kai ratio:** Kᵢ = ω²_p,intra / (ω²_p,intra + ω²_p,inter), N_itinerant = Kᵢ·NELECT
  - **f-sum rule:** N_itinerant = (ε₀ mₑ V ω²_p,intra) / ħ²
- The quantum metric **√G / A** along xx, yy, zz directions.

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

### Choose the electron-counting method

```bash
quantum-metric --method fsum     # f-sum rule
quantum-metric --method kai      # ratio method (default)
```

### Different output formats

```bash
quantum-metric --format table    # pretty table (default)
quantum-metric --format json --output results.json
quantum-metric --format tsv  --output results.tsv
quantum-metric --format csv  --output results.csv
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
result = calc.compute(method="kai")

print(result.metric.sqrtG_over_A_xx)
print(result.electrons.n_bound)
print(result.to_dict())           # flat dict, ready for pandas
```

Or build it from explicit file paths:

```python
calc = QMetricCalculator(
    outcar="./OUTCAR",
    poscar="./POSCAR",
    dielectric="./vasprun.xml",
)
result = calc.compute(method="fsum")
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
4. Compute σ(ω) = (ω / 4π) ε₂(ω) and the integrals `I_xx = ∫σ/ω dω`, etc.
5. Compute bound electron count `N_bound` via the chosen method.
6. Compute the metric

   √G / A = √( prefactor · I_xx / n_bound^(1/3) ) / |a|

   where *prefactor* = 0.0694 Å⁻¹·eV⁻¹ (unit-conversion constant; overridable with `--prefactor`).

## License

MIT
