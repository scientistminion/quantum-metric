# Python API
`quantum-metric` can also be used as a library from Python scripts or Jupyter notebooks. Everything the CLI does is built on top of a small public API.
## High-level: `QMetricCalculator`
The simplest entry point. It auto-discovers files in a directory, runs the full pipeline, and returns a result object you can query or serialize.
```python
from quantum_metric import QMetricCalculator
calc = QMetricCalculator.from_directory("./Ag_fcc")
result = calc.compute()
print(f"Material:        {result.material}")
print(f"N_itinerant:     {result.electrons.n_itinerant:.4f}")
print(f"N_bound:         {result.electrons.n_bound:.4f}")
print(f"Sumrule check:   {result.electrons.sumrule_check:.4f}  (should ≈ NELECT)")
print(f"√G (xx):         {result.metric.sqrtG_xx:.6f}")
```
### Customising the integration window
```python
result = calc.compute(e_min=0.1, e_max=20.0)
```
### Overriding the prefactor
```python
result = calc.compute(prefactor=0.0694)  # default value, shown for reference
```
### Serializing results
Every result has a `.to_dict()` that returns a flat dictionary, perfect for writing to TSV/CSV or feeding into pandas:
```python
import pandas as pd
rows = []
for material in ["Ag_fcc", "Cu_fcc", "Au_fcc"]:
    calc = QMetricCalculator.from_directory(f"./{material}")
    rows.append(calc.compute().to_dict())
df = pd.DataFrame(rows)
df.to_csv("metals.tsv", sep="\t", index=False)
```
## Building a calculator from explicit file paths
If files are spread across multiple directories:
```python
calc = QMetricCalculator(
    outcar="./run_a/OUTCAR",
    poscar="./structures/POSCAR",
    dielectric="./run_b/vasprun.xml",
    material="my_material",
)
result = calc.compute()
```
## Low-level: individual steps
If you want finer control, call each step manually.
```python
from quantum_metric import (
    read_outcar, read_poscar, read_dielectric,
    compute_optical_integrals,
    compute_electron_count,
    compute_quantum_metric,
)
# 1. Parse files
outcar = read_outcar("./OUTCAR")
poscar = read_poscar("./POSCAR")
eps = read_dielectric("./vasprun.xml")
# 2. Integrate the optical conductivity
optical = compute_optical_integrals(eps, e_min=0.0, e_max=20.0)
# 3. Count electrons via the f-sum rule
electrons = compute_electron_count(
    plasma_intra_ev2=outcar.plasma_intra_xx,
    nelect=outcar.nelect,
    volume_ang3=outcar.volume,
    natoms=outcar.natoms,
    sumrule_ev2=outcar.sumrule,
)
# 4. Compute the metric
metric = compute_quantum_metric(
    I_xx=optical.xx.I,
    I_yy=optical.yy.I,
    I_zz=optical.zz.I,
    bound_electron_density=electrons.bound_electron_density,
)
print(metric.sqrtG_xx)
```
## The prefactor as a constant
The library exposes the f-sum prefactor $1/(16\pi a_B^3 E_0^2)$ as a module-level constant for direct use:
```python
from quantum_metric.electrons import PREFACTOR_N
# n [Å⁻³] = PREFACTOR_N × X_vasp [eV²]
# numerically ≈ 7.263e-4 Å⁻³ eV⁻²
n_intra = PREFACTOR_N * 37.394   # Na bcc intraband example
print(f"n_intra = {n_intra:.4e} Å⁻³")
```
See [Theory](theory.md) for the derivation.
## Plotting from Python
```python
from quantum_metric import read_dielectric
from quantum_metric.plotting import plot_optical_conductivity
eps = read_dielectric("./vasprun.xml")
fig, ax = plot_optical_conductivity(eps, output="sigma.png", e_max=20)
```
Both `plot_dielectric` and `plot_optical_conductivity` return a `(fig, ax)` tuple so you can further customise the plot:
```python
fig, ax = plot_optical_conductivity(eps, e_max=20)
ax.axvline(2.5, color="red", linestyle="--", label="Some feature")
ax.legend()
fig.savefig("annotated.png", dpi=300)
```
## Next steps
- {doc}`api/index` for full API reference (every function and class)
- {doc}`theory` for the physics
