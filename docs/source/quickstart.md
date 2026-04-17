# Quickstart

This page walks you through a complete computation on one material, start to finish.

## Prerequisites

You need a VASP calculation directory that contains:

- `OUTCAR` — from an optics run with `LOPTICS = .TRUE.`
- `POSCAR` — the structure
- `vasprun.xml` — (strongly preferred) or a `*_eps_imag.dat` file produced by [sumo](https://smtg-bham.github.io/sumo/)

See the [VASP optics tutorial](https://www.vasp.at/tutorials/latest/response/part1/) if you're unsure how to run one.

## One-liner

```bash
cd my_vasp_run/
quantum-metric
```

This uses sensible defaults: the Kai-ratio method, the whole energy range, and a pretty-printed table.

## Example output

For fcc Ag:

```
         Quantum Metric Results: Ag_fcc
┏━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┓
┃ Quantity                ┃    Value ┃ Units    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━┩
│ Method                  │      kai │          │
│ Volume                  │    16.39 │ Å³       │
│ NELECT                  │       11 │          │
│ NAtoms (NIONS)          │        1 │          │
│ a_len                   │  2.85105 │ Å        │
├─────────────────────────┼──────────┼──────────┤
│ ω²_p (intraband, xx)    │  101.969 │ eV²      │
│ ω²_p (interband, xx)    │  348.767 │ eV²      │
│ Sumrule                 │  925.581 │ eV²      │
├─────────────────────────┼──────────┼──────────┤
│ I_xx = ∫σ/ω dω          │  3.67951 │          │
│ Kai ratio               │ 0.226228 │          │
│ N_itinerant             │  2.48851 │          │
│ N_bound                 │  8.51149 │          │
│ Bound electron density  │  0.51931 │ 1/Å³     │
├─────────────────────────┼──────────┼──────────┤
│ √G  (xx)                │ 0.563642 │          │
└─────────────────────────┴──────────┴──────────┘
```

## Common workflows

### Switch the electron-counting method

```bash
quantum-metric --method fsum
```

See {doc}`theory` for the difference between `kai` and `fsum`.

### Save results to a file

```bash
quantum-metric --format tsv --output ag.tsv
quantum-metric --format json --output ag.json
```

### Override individual files

```bash
quantum-metric --outcar ../other_run/OUTCAR \
               --poscar ./POSCAR \
               --dielectric ./vasprun.xml
```

### Plot the optical conductivity

```bash
quantum-metric plot --kind optics --output sigma.png
quantum-metric plot --kind epsilon --output eps2.png
```

### Check what files will be used

```bash
quantum-metric info
```

Helpful when your directory contains multiple `*_eps_imag.dat` files.

## Next steps

- {doc}`cli` — full CLI reference
- {doc}`python_api` — call from Python scripts and notebooks
- {doc}`theory` — what the numbers mean
- {doc}`examples/index` — worked examples on real materials
