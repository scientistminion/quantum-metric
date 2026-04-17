# Installation

## From PyPI (recommended)

```bash
pip install quantum-metric
```

## From source

```bash
git clone https://github.com/yourusername/quantum-metric.git
cd quantum-metric
pip install -e ".[dev]"
```

The `[dev]` extra installs testing and linting tools. Leave it off for a minimal install.

## Requirements

- Python ≥ 3.9
- `numpy`, `pandas`, `matplotlib`, `typer`, `rich`

All dependencies are installed automatically by `pip`.

## Verifying the installation

```bash
quantum-metric --version
quantum-metric --help
```

If you see the version number and a help message with subcommands, you're good to go.

## Running the tests

If you installed from source with the `[dev]` extra:

```bash
pytest -v
```

All tests should pass. This is a good way to confirm your environment is set up correctly before running on real VASP output.
