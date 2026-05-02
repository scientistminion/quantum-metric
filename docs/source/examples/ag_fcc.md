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
```text
            Quantum Metric Results: Ag_fcc
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┓
┃ Quantity                   ┃    Value ┃ Units    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━┩
│ Volume                     │    16.39 │ Å³       │
│ NELECT                     │       11 │          │
│ NAtoms (NIONS)             │        1 │          │
│ a_len                      │  2.85105 │ Å        │
├────────────────────────────┼──────────┼──────────┤
│ ω²_p (intraband, xx)       │  101.969 │ eV²      │
│ ω²_p (interband, xx)       │  348.767 │ eV²      │
│ Sumrule                    │  925.581 │ eV²      │
├────────────────────────────┼──────────┼──────────┤
│ I_xx = ∫σ/ω dω             │  3.67951 │          │
├────────────────────────────┼──────────┼──────────┤
│ N_itinerant                │  1.21421 │          │
│ N_bound                    │  9.78579 │          │
│ N_itinerant / atom         │  1.21421 │          │
│ N_bound / atom             │  9.78579 │          │
│ Itinerant electron density │ 0.074083 │ 1/Å³     │
│ Bound electron density     │ 0.597059 │ 1/Å³     │
│ Sumrule check (≈ NELECT)   │  11.0241 │          │
├────────────────────────────┼──────────┼──────────┤
│ √G  (xx)                   │ 0.550697 │          │
│ prefactor used             │   0.0694 │ Å⁻¹ eV⁻¹ │
└────────────────────────────┴──────────┴──────────┘
```
## Interpretation
- fcc Ag has cubic symmetry, so $I_{xx} = I_{yy} = I_{zz}$ and $\sqrt{G}$ is identical along all three directions — a good sanity check that the parser is reading the diagonal tensor components correctly.
- $N_{\rm itinerant} \approx 1.2$ from the f-sum rule corresponds to roughly one conduction electron per atom — the $5s^1$ valence electron of silver, with a small renormalization that reflects the deviation of the intraband effective mass from $m_e$.
- $N_{\rm bound} \approx 9.8$ accounts for the filled $4d^{10}$ shell that contributes to interband transitions but not to the Drude weight.
- The `Sumrule check` value (≈ 11.0) recovers `NELECT` = 11 to within ~0.2%, confirming that the optical integration is well-converged and the prefactor is correctly applied.
## Plotting
```bash
quantum-metric plot --kind optics --output sigma_Ag.png --no-show
quantum-metric plot --kind epsilon --output eps2_Ag.png --no-show
```
The σ(ω) plot shows the expected Drude peak below ~1 eV and the strong interband absorption feature near 4 eV.
