# Command-line interface
## `quantum-metric`
The top-level command computes the quantum metric for a single VASP calculation directory.
```bash
quantum-metric [OPTIONS] [COMMAND]
```
By default it processes the current directory. Subcommands are available for plotting and introspection.
### Options
| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--dir`, `-d` | path | cwd | VASP calculation directory |
| `--outcar` | path | auto | Path to OUTCAR |
| `--poscar` | path | auto | Path to POSCAR |
| `--dielectric`, `--eps` | path | auto | Path to vasprun.xml or `*_eps_imag.dat` |
| `--dim` | int | `3` | Spatial dimension d (3 = bulk, 2 = monolayer, …) for the dimensionless κ |
| `--e-min` | float | `0.0` | Lower integration bound (eV) |
| `--e-max` | float | ∞ | Upper integration bound (eV) |
| `--format`, `-f` | `table` \| `json` \| `tsv` \| `csv` | `table` | Output format |
| `--output`, `-o` | path | stdout | Write output to file instead of stdout |
| `--version`, `-V` | — | — | Show version and exit |
| `--help`, `-h` | — | — | Show help and exit |

The quantum metric is computed directly from the Souza–Wilkens–Martin sum rule, with all fundamental constants (ℏ, e, ε₀) in SI. The output is the per-electron metric tensor `g_µµ` in Å² and the dimensionless ratio `κ_µ = n_bound^{-(1/2 − 1/d)} √g_µµ`. Electron counting uses the f-sum rule applied to the intraband plasma frequency, with `N_bound = NELECT − N_itinerant`. See [Theory](theory.md) for the full derivation.

### Examples
Run on the current directory with defaults:
```bash
quantum-metric
```
Save results as TSV:
```bash
quantum-metric --format tsv --output results.tsv
```
Limit integration to the 0–15 eV window:
```bash
quantum-metric --e-max 15
```
Run on a 2D monolayer:
```bash
quantum-metric --dim 2
```
Point at a different directory and override a specific file:
```bash
quantum-metric --dir /path/to/run --outcar /other/path/OUTCAR
```
## `quantum-metric plot`
Plot the optical conductivity σ(ω) or the imaginary dielectric function ε₂(ω).
```bash
quantum-metric plot [OPTIONS]
```
### Options
| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--kind`, `-k` | `optics` \| `epsilon` | `optics` | What to plot |
| `--dir`, `-d` | path | cwd | VASP calculation directory |
| `--dielectric`, `--eps` | path | auto | Path to a specific dielectric file |
| `--output`, `-o` | path | — | Save plot to file instead of showing |
| `--e-max` | float | — | Upper x-axis bound (eV) |
| `--show` / `--no-show` | flag | `--show` | Open interactive window (use `--no-show` on headless) |
### Examples
Save σ(ω) as a PNG on a headless compute node:
```bash
quantum-metric plot --kind optics --output sigma.png --no-show
```
Plot ε₂(ω) zoomed into the low-energy range:
```bash
quantum-metric plot --kind epsilon --e-max 10
```
## `quantum-metric info`
Inspect a directory and report which VASP files are auto-discovered.
```bash
quantum-metric info [OPTIONS]
```
### Example
```bash
quantum-metric info --dir /path/to/run
```
