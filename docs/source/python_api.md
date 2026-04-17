# Python API

`quantum-metric` can also be used as a library from Python scripts or Jupyter notebooks. Everything the CLI does is built on top of a small public API.

## High-level: `QMetricCalculator`

The simplest entry point. It auto-discovers files in a directory, runs the full pipeline, and returns a result object you can query or serialize.

```python
from quantum_metric import QMetricCalculator

calc = QMetricCalculator.from_directory("./Ag_fcc")
result = calc.compute(method="kai")

print(f"Material:     {result.material}")
print(f"Kai ratio:    {result.electrons.kai:.4f}")
print(f"N_bound:      {result.electrons.n_bound:.4f}")
print(f"√G (xx):      {result.metric.sqrtG_xx:.6f}")
```

### Switching methods

```python
result_kai  = calc.compute(method="kai")
result_fsum = calc.compute(method="fsum")
```

### Customising the integration window

```python
result = calc.compute(method="kai", e_min=0.1, e_max=20.0)
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
    compute_n_itinerant_kai,
    compute_quantum_metric,
)

# 1. Parse files
outcar = read_outcar("./OUTCAR")
poscar = read_poscar("./POSCAR")
eps = read_dielectric("./vasprun.xml")

# 2. Integrate the optical conductivity
optical = compute_optical_integrals(eps, e_min=0.0, e_max=20.0)

# 3. Count electrons
electrons = compute_n_itinerant_kai(
    plasma_intra=outcar.plasma_intra_xx,
    plasma_inter=outcar.plasma_inter_xx,
    nelect=outcar.nelect,
    volume=outcar.volume,
    natoms=outcar.natoms,
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
