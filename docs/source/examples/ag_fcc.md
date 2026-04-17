# fcc Silver (Ag)

Ag is a good test case: a simple cubic metal with a well-known Drude response and a clean interband onset around 4 eV.

## Setting up the VASP calculation

Brief INCAR flags for the optics run:

```text
LOPTICS = .TRUE.
LPEAD   = .TRUE.
NEDOS   = 2000
SIGMA   = 0.01
```

Run VASP. Afterward, the directory should contain `OUTCAR`, `POSCAR`, `vasprun.xml`, and (optionally) a `*_eps_imag.dat` file from `sumo-optplot`.

## Running `quantum-metric`

```bash
cd Ag_fcc/
quantum-metric
```

## Results

### Kai method (default)

```text
         Quantum Metric Results: Ag_fcc
┏━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┓
┃ Quantity                ┃    Value ┃ Units    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━┩
│ Method                  │      kai │          │
│ Volume                  │    16.39 │ Å³       │
│ NELECT                  │       11 │          │
│ NAtoms (NIONS)          │        1 │          │
│ a_len                   │  2.85105 │ Å        │
│ ω²_p (intraband, xx)    │  101.969 │ eV²      │
│ ω²_p (interband, xx)    │  348.767 │ eV²      │
│ I_xx = ∫σ/ω dω          │  3.67951 │          │
│ Kai ratio               │ 0.226228 │          │
│ N_itinerant             │  2.48851 │          │
│ N_bound                 │  8.51149 │          │
│ Bound electron density  │  0.51931 │ 1/Å³     │
│ √G (xx)                 │ 0.563642 │          │
└─────────────────────────┴──────────┴──────────┘
```

### f-sum method

```bash
quantum-metric --method fsum
```

```text
│ N_itinerant             │  1.21208 │          │
│ N_bound                 │  9.78792 │          │
│ Bound electron density  │ 0.597188 │ 1/Å³     │
│ √G (xx)                 │ 0.550667 │          │
```

## Interpretation

- fcc Ag has cubic symmetry, so $I_{xx} = I_{yy} = I_{zz}$ and $\sqrt{G}$ is identical along all three directions — a good sanity check that the parser is reading the diagonal tensor components correctly.
- Kai ≈ 0.23 reflects a substantial **interband** contribution to the total plasma weight, as expected for Ag where the d→s transitions kick in around 4 eV.
- The f-sum method gives a smaller $N_{\rm itinerant}$ (~1.2 vs ~2.5) because it uses the absolute intraband plasma frequency rather than the Kai ratio.

## Plotting

```bash
quantum-metric plot --kind optics --output sigma_Ag.png --no-show
quantum-metric plot --kind epsilon --output eps2_Ag.png --no-show
```

The σ(ω) plot shows the expected Drude peak below ~1 eV and the strong interband absorption feature near 4 eV.
