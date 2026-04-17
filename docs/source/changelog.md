# Changelog

## 0.1.0 — initial release

- Single-directory workflow: `quantum-metric` runs on any VASP directory with `LOPTICS=.TRUE.` output
- CLI commands: default compute, `plot`, `info`
- Two electron-counting methods: `kai` (default) and `fsum`
- Reads dielectric function from `vasprun.xml` (no `sumo` dependency) with fallback to `*_eps_imag.dat`
- Output formats: rich table, JSON, TSV, CSV
- Python API via `QMetricCalculator`
- Plotting helpers for σ(ω) and ε₂(ω)
- Validated against reference data for fcc Ag
