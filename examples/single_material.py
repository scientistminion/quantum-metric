"""
Example: compute the quantum metric for a single VASP calculation.

Assumes you are inside a directory that contains:
  - OUTCAR (from an LOPTICS=.TRUE. run)
  - POSCAR
  - vasprun.xml (or a *_eps_imag.dat file from sumo-optplot)
"""

from quantum_metric import QMetricCalculator

# Auto-discover files in the current directory
calc = QMetricCalculator.from_directory(".")

# Compute using the Kai-ratio method (default)
result = calc.compute(method="kai")

print(f"Material: {result.material}")
print(f"Volume:   {result.volume:.2f} Å³")
print(f"Kai:      {result.electrons.kai:.4f}")
print(f"N_bound:  {result.electrons.n_bound:.4f}")
print(f"√G / A (xx): {result.metric.sqrtG_over_A_xx:.6f}")

if result.metric.sqrtG_over_A_yy is not None:
    print(f"√G / A (yy): {result.metric.sqrtG_over_A_yy:.6f}")
    print(f"√G / A (zz): {result.metric.sqrtG_over_A_zz:.6f}")

# Switch to f-sum rule method
result_fsum = calc.compute(method="fsum")
print(f"\nf-sum: N_itinerant = {result_fsum.electrons.n_itinerant:.4f}")

# Dump to a flat dict (e.g. for pandas)
import pandas as pd
df = pd.DataFrame([result.to_dict()])
df.to_csv("results.tsv", sep="\t", index=False)
print("\nSaved results.tsv")
