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
| `--method`, `-m` | `kai` \| `fsum` | `kai` | Electron-counting method |
| `--prefactor` | float | `0.0694` | Unit-conversion constant (Å⁻¹ eV⁻¹) |
| `--e-min` | float | `0.0` | Lower integration bound (eV) |
| `--e-max` | float | ∞ | Upper integration bound (eV) |
| `--format`, `-f` | `table` \| `json` \| `tsv` \| `csv` | `table` | Output format |
| `--output`, `-o` | path | stdout | Write output to file instead of stdout |
| `--version`, `-V` | — | — | Show version and exit |
| `--help`, `-h` | — | — | Show help and exit |

### Examples

Run on the current directory with defaults:

```bash
quantum-metric
```

Use the f-sum rule method and save as TSV:

```bash
quantum-metric --method fsum --format tsv --output results.tsv
```

Limit integration to the 0–15 eV window:

```bash
quantum-metric --e-max 15
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

Inspect a directory and report which VASP files are auto-discovered. Handy when multiple `*_eps_imag.dat` files exist or when you want to verify which `OUTCAR` / `POSCAR` / `vasprun.xml` will be used.

```bash
quantum-metric info [OPTIONS]
```

### Example

```bash
quantum-metric info --dir /path/to/run
```
