# quantum-metric

A Python library and command-line tool for computing the **quantum metric** and related optical quantities from [VASP](https://www.vasp.at/) output.

[![PyPI](https://img.shields.io/pypi/v/quantum-metric)](https://pypi.org/project/quantum-metric/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/docs-latest-blue)](https://yourusername.github.io/quantum-metric/)

## What it does

Given a VASP optical calculation run with `LOPTICS = .TRUE.`, `quantum-metric`:

- Parses `OUTCAR` for intraband / interband plasma frequencies, sumrule, volume, NELECT, NIONS
- Reads the frequency-dependent dielectric function ε₂(ω) directly from `vasprun.xml` (no `sumo` needed) or from pre-computed `*_eps_imag.dat`
- Computes optical conductivity integrals: ∫σ(ω)/ω dω, ∫σ(ω) dω, ∫ωσ(ω) dω, …
- Computes bound / itinerant electron counts via two methods (Kai ratio or f-sum rule)
- Computes the **quantum metric** $\sqrt{G}$ along xx, yy, zz directions

## Quick start

```bash
pip install quantum-metric
```

Then, in any directory containing VASP output:

```bash
cd my_vasp_run/
quantum-metric
```

You'll get a pretty table with all the quantities computed.

## Documentation

```{toctree}
:maxdepth: 2
:caption: Getting Started

installation
quickstart
```

```{toctree}
:maxdepth: 2
:caption: User Guide

cli
python_api
theory
```

```{toctree}
:maxdepth: 2
:caption: Examples

examples/index
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
changelog
```

## Citation

If you use `quantum-metric` in your research, please cite:

```
Your Name. (2026). quantum-metric: A Python package for the quantum metric from VASP.
GitHub: https://github.com/yourusername/quantum-metric
```

## License

MIT License — see [LICENSE](https://github.com/yourusername/quantum-metric/blob/main/LICENSE).
