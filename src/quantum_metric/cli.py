"""
Command-line interface for quantum-metric.

Usage:
    quantum-metric                                    # current directory
    quantum-metric --dir path/to/vasp/run
    quantum-metric --outcar OUTCAR --poscar POSCAR --dielectric vasprun.xml
    quantum-metric --format json
    quantum-metric plot --kind optics
"""

from __future__ import annotations
import os
import json
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from quantum_metric._version import __version__
from quantum_metric.calculator import QMetricCalculator, QMetricResult
from quantum_metric.kresolved_plot import plot_band_with_metric
 
 
 

app = typer.Typer(
    help="Compute the quantum metric and optical quantities from VASP output.",
    add_completion=False,
    no_args_is_help=False,
    context_settings={"help_option_names": ["-h", "--help"]}
    )
console = Console()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _version_callback(value: bool):
    if value:
        console.print(f"quantum-metric {__version__}")
        raise typer.Exit()


def _build_calculator(
    directory: Optional[Path],
    outcar: Optional[Path],
    poscar: Optional[Path],
    dielectric: Optional[Path],
) -> QMetricCalculator:
    """Decide between from_directory and explicit-files modes.

    If any individual file is passed, they override the directory-based discovery.
    """
    base_dir = directory or Path.cwd()

    if outcar is None and poscar is None and dielectric is None:
        # Pure directory mode
        return QMetricCalculator.from_directory(base_dir)

    # Mixed mode: start from directory defaults, override individual files
    calc_dir = QMetricCalculator.from_directory(base_dir)
    return QMetricCalculator(
        outcar=outcar or calc_dir.outcar_path,
        poscar=poscar or calc_dir.poscar_path,
        dielectric=dielectric or calc_dir.dielectric_path,
        dielectric_source="auto",
        material=base_dir.name,
    )


def _print_table(result: QMetricResult):
    """Pretty-print a single-material result as a rich table."""
    table = Table(
        title=f"[bold cyan]Quantum Metric Results: {result.material}[/bold cyan]",
        show_header=True,
        header_style="bold magenta",
    )
    table.add_column("Quantity", style="cyan", no_wrap=True)
    table.add_column("Value", justify="right")
    table.add_column("Units", style="dim")

    def row(name, value, units=""):
        if value is None:
            table.add_row(name, "[dim]—[/dim]", units)
        elif isinstance(value, float):
            table.add_row(name, f"{value:.6g}", units)
        else:
            table.add_row(name, str(value), units)

    # --- inputs
    table.add_section()
    row("Volume", result.volume, "Å³")
    row("NELECT", result.nelect)
    row("NAtoms (NIONS)", result.natoms)
    row("a_len", result.a_len, "Å")
    row("dimension d", result.metric.dim)

    # --- plasma & sumrule
    table.add_section()
    row("ω²_p (intraband, xx)", result.plasma_intra_xx, "eV²")
    row("ω²_p (interband, xx)", result.plasma_inter_xx, "eV²")
    row("Sumrule", result.sumrule, "eV²")

    # --- optical integrals
    table.add_section()
    row("ω²_p from integral (xx)", result.optical.xx.omega_p_squared, "eV²")
    row("I_xx = ∫σ/ω dω", result.optical.xx.I)
    row("∫σ_xx dω", result.optical.xx.sigma_int)
    row("∫ω·σ_xx dω", result.optical.xx.wsigma)

    if result.optical.yy is not None:
        row("I_yy = ∫σ/ω dω", result.optical.yy.I)
    if result.optical.zz is not None:
        row("I_zz = ∫σ/ω dω", result.optical.zz.I)

    # --- electrons (f-sum rule via hydrogen relations)
    table.add_section()
    row("N_itinerant", result.electrons.n_itinerant)
    row("N_bound", result.electrons.n_bound)
    row("N_itinerant / atom", result.electrons.n_itinerant_per_atom)
    row("N_bound / atom", result.electrons.n_bound_per_atom)
    row("Itinerant electron density", result.electrons.itinerant_electron_density, "1/Å³")
    row("Bound electron density", result.electrons.bound_electron_density, "1/Å³")
    row("Sumrule check (≈ NELECT)", result.electrons.sumrule_check)

    # --- quantum metric (g in Å² and dimensionless κ)
    table.add_section()
    row("[bold green]g_xx[/bold green]", result.metric.g_xx, "Å²")
    if result.metric.g_yy is not None:
        row("[bold green]g_yy[/bold green]", result.metric.g_yy, "Å²")
    if result.metric.g_zz is not None:
        row("[bold green]g_zz[/bold green]", result.metric.g_zz, "Å²")
    row("[bold green]κ_xx[/bold green]", result.metric.kappa_xx)
    if result.metric.kappa_yy is not None:
        row("[bold green]κ_yy[/bold green]", result.metric.kappa_yy)
    if result.metric.kappa_zz is not None:
        row("[bold green]κ_zz[/bold green]", result.metric.kappa_zz)

    console.print(table)


def _emit_result(result: QMetricResult, fmt: str, output: Optional[Path]):
    """Emit the result in the requested format."""
    if fmt == "table":
        _print_table(result)
        return

    d = result.to_dict()

    if fmt == "json":
        text = json.dumps(d, indent=2, default=float)
    elif fmt == "tsv":
        keys = list(d.keys())
        vals = [str(d[k]) if d[k] is not None else "" for k in keys]
        text = "\t".join(keys) + "\n" + "\t".join(vals) + "\n"
    elif fmt == "csv":
        keys = list(d.keys())
        vals = [str(d[k]) if d[k] is not None else "" for k in keys]
        text = ",".join(keys) + "\n" + ",".join(vals) + "\n"
    else:
        raise typer.BadParameter(f"Unknown format: {fmt}")

    if output:
        Path(output).write_text(text)
        console.print(f"[green]✓[/green] Wrote {output}")
    else:
        console.print(text)


# ---------------------------------------------------------------------------
# Main "compute" command — default if no subcommand given
# ---------------------------------------------------------------------------
@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    directory: Optional[Path] = typer.Option(
        None, "--dir", "-d", help="VASP calculation directory (default: current dir)."
    ),
    outcar: Optional[Path] = typer.Option(
        None, "--outcar", help="Path to OUTCAR (overrides directory default)."
    ),
    poscar: Optional[Path] = typer.Option(
        None, "--poscar", help="Path to POSCAR (overrides directory default)."
    ),
    dielectric: Optional[Path] = typer.Option(
        None,
        "--dielectric",
        "--eps",
        help="Path to vasprun.xml or *_eps_imag.dat (overrides auto-discovery).",
    ),
    dim: int = typer.Option(
        3,
        "--dim",
        help="Spatial dimension d for the dimensionless κ normalization (default 3).",
    ),
    e_min: float = typer.Option(0.0, "--e-min", help="Lower integration bound (eV)."),
    e_max: Optional[float] = typer.Option(
        None, "--e-max", help="Upper integration bound (eV)."
    ),
    fmt: str = typer.Option(
        "table", "--format", "-f", help="Output format: table | json | tsv | csv."
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Write output to file instead of stdout."
    ),
    version: Optional[bool] = typer.Option(
        None, "--version", "-V", callback=_version_callback, is_eager=True, help="Show version."
    ),
):
    """
    Compute the quantum metric and optical quantities from a VASP calculation.

    By default, runs on the current directory. Override the directory with --dir,
    or override individual files with --outcar / --poscar / --dielectric.

    The quantum metric is computed directly from the Souza-Wilkens-Martin sum rule:

        g_µν = (ℏ / π e²) · (1 / n_bound) · ∫ σ_µν(ω)/ω dω        (units of Å²)
        κ_µ  = n_bound^{-(1/2 - 1/d)} · √g_µµ                      (dimensionless)

    with all fundamental constants in SI — no fitted prefactors. Electron
    counting uses the f-sum rule applied to the intraband plasma frequency,
    with N_bound = NELECT − N_itinerant.

    \b
    Examples:
      quantum-metric                              # current dir, pretty table
      quantum-metric --dir /path/to/run           # different directory
      quantum-metric --dim 2                      # 2D system (e.g. monolayer)
      quantum-metric --format tsv -o results.tsv  # dump TSV to file
      quantum-metric --format json                # print JSON to stdout
      quantum-metric --e-max 15                   # integrate ε₂ only up to 15 eV
      quantum-metric --outcar ../run2/OUTCAR      # use OUTCAR from elsewhere

    \b
    Subcommands:
      plot   Plot the optical conductivity or dielectric function
      info   Inspect which VASP files are auto-discovered in a directory
    """

    if ctx.invoked_subcommand is not None:
        return

    try:
        calc = _build_calculator(directory, outcar, poscar, dielectric)
        result = calc.compute(dim=dim, e_min=e_min, e_max=e_max)
    except FileNotFoundError as e:
        console.print(f"[red]✗[/red] {e}")
        raise typer.Exit(code=1)
    except ValueError as e:
        console.print(f"[red]✗[/red] {e}")
        raise typer.Exit(code=1)

    _emit_result(result, fmt, output)


# ---------------------------------------------------------------------------
# `quantum-metric plot` subcommand
# ---------------------------------------------------------------------------
@app.command("plot")
def plot_cmd(
    kind: str = typer.Option(
        "optics", "--kind", "-k", help="What to plot: 'optics' (sigma) or 'epsilon' (eps_2)."
    ),
    directory: Optional[Path] = typer.Option(None, "--dir", "-d"),
    dielectric: Optional[Path] = typer.Option(None, "--dielectric", "--eps"),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Save plot to file instead of showing."
    ),
    e_max: Optional[float] = typer.Option(None, "--e-max", help="Upper x-axis bound (eV)."),
    show: bool = typer.Option(True, "--show/--no-show", help="Open in window."),
):
    """
    Plot the optical conductivity σ(ω) or imaginary dielectric function ε₂(ω).

    \b
    Examples:
      quantum-metric plot                               # σ(ω) from current dir, show window
      quantum-metric plot --kind epsilon                # ε₂(ω) instead of σ(ω)
      quantum-metric plot --output sigma.png --no-show  # save to PNG (good for headless)
      quantum-metric plot --e-max 10                    # zoom into 0–10 eV
      quantum-metric plot --dir /path/to/run --output sigma.png
      quantum-metric plot --dielectric ./my_eps.dat     # use a specific eps_imag.dat file
    """
    from quantum_metric.io import find_vasp_files, read_dielectric
    from quantum_metric.plotting import plot_dielectric, plot_optical_conductivity

    base_dir = directory or Path.cwd()

    if dielectric is None:
        files = find_vasp_files(base_dir)
        dielectric = files["vasprun"] or files["eps_dat"]
        if dielectric is None:
            console.print(f"[red]✗[/red] No dielectric data found in {base_dir}")
            raise typer.Exit(code=1)

    data = read_dielectric(dielectric)

    title = base_dir.name
    if kind == "optics":
        plot_optical_conductivity(data, output=output, title=title, show=show, e_max=e_max)
    elif kind == "epsilon":
        plot_dielectric(data, output=output, title=title, show=show, e_max=e_max)
    else:
        console.print(f"[red]✗[/red] Unknown plot kind: {kind}. Use 'optics' or 'epsilon'.")
        raise typer.Exit(code=1)

    if output:
        console.print(f"[green]✓[/green] Plot saved to {output}")


# ---------------------------------------------------------------------------
# `quantum-metric info` — inspect what files are found
# ---------------------------------------------------------------------------
@app.command("info")
def info_cmd(
    directory: Optional[Path] = typer.Option(None, "--dir", "-d"),
):
    """
    Inspect a directory and report which VASP files are auto-discovered.

    \b
    Useful when a directory contains multiple *_eps_imag.dat files or when you want
    to confirm which OUTCAR / POSCAR / vasprun.xml will be used before running compute.

    \b
    Examples:
      quantum-metric info
      quantum-metric info --dir /path/to/other/run
    """
    from quantum_metric.io import find_vasp_files

    base_dir = (directory or Path.cwd()).resolve()
    files = find_vasp_files(base_dir)

    table = Table(title=f"VASP files in {base_dir}")
    table.add_column("File", style="cyan")
    table.add_column("Found?", justify="center")
    table.add_column("Path", style="dim")

    for key, path in files.items():
        mark = "[green]✓[/green]" if path else "[red]✗[/red]"
        table.add_row(key, mark, str(path) if path else "—")

    console.print(table)

@app.command("kpath")
def kpath_cmd(
    waveder: Path = typer.Option(
        "WAVEDER", "--waveder", "-w",
        help="Path to VASP WAVEDER (binary, from LOPTICS=.TRUE.).",
    ),
    eigenval: Path = typer.Option(
        "EIGENVAL", "--eigenval", "-e",
        help="Path to VASP EIGENVAL.",
    ),
    kpoints: Path = typer.Option(
        "KPOINTS", "--kpoints", "-k",
        help="Path to line-mode KPOINTS (with high-symmetry labels in '!' comments).",
    ),
    outcar: Optional[Path] = typer.Option(
        "OUTCAR", "--outcar",
        help="Path to OUTCAR (used for reciprocal lattice → correct k-distances). "
             "Optional but recommended.",
    ),
    output: Path = typer.Option(
        "bands_with_gk.png", "--output", "-o",
        help="Output PNG file path.",
    ),
    per_electron: bool = typer.Option(
        False, "--per-electron",
        help="Plot g(k) / NELECT instead of total g(k).",
    ),
    trace_only: bool = typer.Option(
        False, "--trace-only",
        help="Plot only Tr[g(k)], hide individual g_xx, g_yy, g_zz "
             "(cleanest for cubic crystals).",
    ),
    ymin: float = typer.Option(-5.0, "--ymin", help="Band-plot lower y limit (E-Ef, eV)."),
    ymax: float = typer.Option(5.0, "--ymax", help="Band-plot upper y limit (E-Ef, eV)."),
    title: Optional[str] = typer.Option(None, "--title", help="Plot title (e.g. material name)."),
    width: float = typer.Option(8.0, "--width", help="Figure width (inches)."),
    height: float = typer.Option(8.0, "--height", help="Figure height (inches)."),
):
    """
    Plot the band structure with the k-resolved quantum metric g_µν(k) on a
    line-mode k-path.

    Computes g_µν(k) = Re Σ_{n occ, m unocc} <u_n|∂_µ|u_m><u_m|∂_ν|u_n> from
    VASP's WAVEDER (wavefunction-derivative matrix elements) and overlays it
    under the band structure read from EIGENVAL/KPOINTS.

    \b
    Examples:
      quantum-metric kpath                                # all defaults
      quantum-metric kpath --per-electron --title "PbTe"  # per-electron g(k)
      quantum-metric kpath --trace-only -o trace.png      # only Tr[g] curve
      quantum-metric kpath --ymin -7 --ymax 7             # wider band window
      quantum-metric kpath --waveder ../other/WAVEDER     # files from elsewhere
    """
    outcar_path = str(outcar) if outcar and os.path.exists(str(outcar)) else None
    if outcar_path is None and outcar is not None:
        console.print(f"Note: {outcar} not found; using fractional k-distances.")

    out = plot_band_with_metric(
        waveder=str(waveder),
        eigenval=str(eigenval),
        kpoints_file=str(kpoints),
        outcar=outcar_path,
        output=str(output),
        per_electron=per_electron,
        trace_only=trace_only,
        ymin=ymin, ymax=ymax,
        title=title,
        width=width, height=height,
    )
    console.print(f"[green]✓[/green] Saved to {out}")

if __name__ == "__main__":
    app()
