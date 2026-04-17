"""
File I/O for VASP output.

Parses OUTCAR for plasma frequencies, volume, NELECT, NIONS.
Parses POSCAR for lattice vectors.
Parses dielectric function from either eps_imag.dat (sumo-style) or vasprun.xml.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------
@dataclass
class OutcarData:
    """Values extracted from an OUTCAR file."""

    plasma_intra_xx: float  # eV^2, intraband plasma freq squared (xx component)
    plasma_inter_xx: float  # eV^2, interband plasma freq squared (xx component)
    sumrule: float          # eV^2
    volume: float           # Angstrom^3
    nelect: float           # number of valence electrons
    natoms: int             # NIONS
    # Full diagonal for multi-direction mode
    plasma_intra_diag: Optional[np.ndarray] = None  # (3,) xx, yy, zz
    plasma_inter_diag: Optional[np.ndarray] = None


@dataclass
class PoscarData:
    """Values extracted from a POSCAR file."""

    a_len: float                          # |a|  (Angstrom)
    lattice: np.ndarray                   # 3x3 lattice matrix (Angstrom)
    lattice_lengths: np.ndarray           # (3,) |a|, |b|, |c|


@dataclass
class DielectricData:
    """Frequency-dependent imaginary dielectric function eps_2(omega)."""

    energy: np.ndarray                     # (N,) eV
    eps_imag_xx: np.ndarray                # (N,)
    eps_imag_yy: Optional[np.ndarray] = None
    eps_imag_zz: Optional[np.ndarray] = None

    @property
    def has_anisotropic(self) -> bool:
        return self.eps_imag_yy is not None and self.eps_imag_zz is not None


# -----------------------------------------------------------------------------
# OUTCAR parsing
# -----------------------------------------------------------------------------
def read_outcar(path: str | Path) -> OutcarData:
    """
    Parse an OUTCAR file for plasma frequencies, sumrule, volume, NELECT, NIONS.

    VASP writes the plasma frequency squared tensor in the form of a 3x3 matrix;
    we extract the xx component (and the full diagonal for anisotropic analysis).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"OUTCAR not found: {path}")

    text = path.read_text()

    plasma_intra_xx, plasma_intra_diag = _parse_plasma_block(
        text, "plasma frequency squared (from intraband transitions)"
    )
    plasma_inter_xx, plasma_inter_diag = _parse_plasma_block(
        text, "plasma frequency squared (from interband transitions"
    )

    sumrule = _parse_sumrule(text)
    volume = _parse_volume(text)
    nelect = _parse_nelect(text)
    natoms = _parse_nions(text)

    return OutcarData(
        plasma_intra_xx=plasma_intra_xx,
        plasma_inter_xx=plasma_inter_xx,
        sumrule=sumrule,
        volume=volume,
        nelect=nelect,
        natoms=natoms,
        plasma_intra_diag=plasma_intra_diag,
        plasma_inter_diag=plasma_inter_diag,
    )


def _parse_plasma_block(text: str, header: str) -> tuple[float, np.ndarray]:
    """Extract the plasma-frequency-squared tensor block following `header`.

    VASP emits two possible layouts (seen in OUTCAR files):

    (A) 3-column diagonal form, typical for intraband or non-SOC:
            xx  xy  xz
            yx  yy  yz
            zx  zy  zz
        -> diagonal at tensor[i, i]

    (B) 6-column form, seen for interband in SOC calculations:
            <xx_Re> <xx_Im> <xy_Re> <xy_Im> <xz_Re> <xz_Im>
            <yx_Re> <yx_Im> <yy_Re> <yy_Im> <yz_Re> <yz_Im>
            <zx_Re> <zx_Im> <zy_Re> <zy_Im> <zz_Re> <zz_Im>
        -> diagonal at tensor[i, 2*i]

    Returns (xx_component, diagonal_array).
    """
    idx = text.find(header)
    if idx < 0:
        raise ValueError(f"Could not find '{header}' in OUTCAR")

    tail = text[idx:]
    lines = tail.splitlines()

    numeric_rows = []
    for line in lines[1:]:
        # stop if we hit a blank line AFTER collecting some rows
        if not line.strip() and numeric_rows:
            break
        # dashed separator lines — skip
        if set(line.strip()) <= {"-"} and line.strip():
            continue
        parts = line.split()
        try:
            row = [float(x) for x in parts]
        except ValueError:
            # non-numeric line; stop if we've already collected, else keep scanning
            if numeric_rows:
                break
            continue
        # Accept rows with at least 3 numeric fields.
        if len(row) >= 3:
            numeric_rows.append(row)
        if len(numeric_rows) == 3:
            break

    if len(numeric_rows) < 3:
        raise ValueError(f"Could not parse plasma-frequency tensor after '{header}'")

    # Enforce consistent row length; pad or trim to the minimum length for safety.
    row_len = min(len(r) for r in numeric_rows)
    tensor = np.array([r[:row_len] for r in numeric_rows])

    if row_len == 3:
        # Layout A: diagonal is [0,0], [1,1], [2,2]
        diag = np.array([tensor[0, 0], tensor[1, 1], tensor[2, 2]])
    elif row_len >= 6 and row_len % 2 == 0:
        # Layout B: diagonal is [0,0], [1,2], [2,4] (stride 2)
        diag = np.array([tensor[0, 0], tensor[1, 2], tensor[2, 4]])
    else:
        raise ValueError(
            f"Unexpected plasma-frequency tensor row length {row_len} "
            f"after '{header}'. Expected 3 or 6-column layout."
        )

    return float(diag[0]), diag

def _parse_sumrule(text: str) -> float:
    """Extract the sumrule value (total plasma freq squared).

    VASP prints lines like:
        sumrule: sum of plasma frequencies squared should yield:                 925.581
        sumrule: sum of plasma frequencies squared should yield (valence only):  925.581
    We match everything up to the last colon, then grab the first number after it.
    """
    match = re.search(
        r"sumrule:\s+sum of plasma frequencies squared should yield[^:]*:\s*([-\d\.Ee+]+)",
        text,
    )
    if not match:
        raise ValueError("Could not find sumrule in OUTCAR")
    return float(match.group(1))

def _parse_volume(text: str) -> float:
    """Extract the last reported cell volume in Angstrom^3."""
    matches = re.findall(r"volume of cell\s*:\s*([-\d\.E+]+)", text)
    if not matches:
        raise ValueError("Could not find cell volume in OUTCAR")
    return float(matches[-1])


def _parse_nelect(text: str) -> float:
    """Extract NELECT."""
    match = re.search(r"NELECT\s*=\s*([-\d\.E+]+)", text)
    if not match:
        raise ValueError("Could not find NELECT in OUTCAR")
    return float(match.group(1))


def _parse_nions(text: str) -> int:
    """Extract NIONS (number of atoms)."""
    match = re.search(r"NIONS\s*=\s*(\d+)", text)
    if not match:
        raise ValueError("Could not find NIONS in OUTCAR")
    return int(match.group(1))


# -----------------------------------------------------------------------------
# POSCAR parsing
# -----------------------------------------------------------------------------
def read_poscar(path: str | Path) -> PoscarData:
    """Parse POSCAR to get lattice vectors and their lengths.

    Only reads the first 5 lines — we do not care about atom species / positions here.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"POSCAR not found: {path}")

    with path.open() as f:
        f.readline()  # comment
        scale_line = f.readline().strip()
        scale = float(scale_line.split()[0])
        a_vec = np.array([float(x) for x in f.readline().split()[:3]])
        b_vec = np.array([float(x) for x in f.readline().split()[:3]])
        c_vec = np.array([float(x) for x in f.readline().split()[:3]])

    # Negative scale means "scale so that volume = |scale|"; handle the common positive case.
    if scale < 0:
        raw_vol = abs(np.dot(a_vec, np.cross(b_vec, c_vec)))
        factor = (abs(scale) / raw_vol) ** (1.0 / 3.0)
    else:
        factor = scale

    lattice = np.array([a_vec, b_vec, c_vec]) * factor
    lengths = np.linalg.norm(lattice, axis=1)

    return PoscarData(
        a_len=float(lengths[0]),
        lattice=lattice,
        lattice_lengths=lengths,
    )


# -----------------------------------------------------------------------------
# Dielectric function parsing
# -----------------------------------------------------------------------------
def read_dielectric(
    path: str | Path,
    *,
    source: str = "auto",
) -> DielectricData:
    """Read imaginary dielectric function from an eps_imag.dat file or vasprun.xml.

    Parameters
    ----------
    path : path to the file
    source : 'auto' | 'dat' | 'vasprun'
        'auto' infers from the file extension / name.

    For eps_imag.dat (sumo-style): expects columns [energy, xx, yy, zz] or [energy, xx].
    For vasprun.xml: parses the <dielectricfunction> block and extracts eps_imag diagonals.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Dielectric file not found: {path}")

    if source == "auto":
        if path.suffix == ".xml" or path.name.endswith("vasprun.xml"):
            source = "vasprun"
        else:
            source = "dat"

    if source == "dat":
        return _read_dielectric_dat(path)
    if source == "vasprun":
        return _read_dielectric_vasprun(path)
    raise ValueError(f"Unknown source: {source!r}")


def _read_dielectric_dat(path: Path) -> DielectricData:
    """Read a sumo-style eps_imag.dat file."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        raise ValueError(f"{path}: expected at least 2 columns")

    energy = data[:, 0]
    eps_xx = data[:, 1]
    eps_yy = data[:, 2] if data.shape[1] > 3 else None
    eps_zz = data[:, 3] if data.shape[1] > 3 else None

    return DielectricData(
        energy=energy,
        eps_imag_xx=eps_xx,
        eps_imag_yy=eps_yy,
        eps_imag_zz=eps_zz,
    )


def _read_dielectric_vasprun(path: Path) -> DielectricData:
    """Parse <dielectricfunction> from vasprun.xml.

    The block has two <imag><array>...</array></imag> and <real> counterparts.
    Each row of the imag <set> has: energy xx yy zz xy yz zx.
    We take the diagonal (xx, yy, zz).
    """
    tree = ET.parse(path)
    root = tree.getroot()

    # VASP emits one or more <dielectricfunction> blocks. We take the first one
    # (typically the density-density response for LOPTICS).
    df = root.find(".//dielectricfunction")
    if df is None:
        raise ValueError(f"No <dielectricfunction> in {path}")

    imag = df.find("./imag/array/set")
    if imag is None:
        raise ValueError(f"No imaginary part found in {path}")

    rows = []
    for r in imag.findall("r"):
        parts = r.text.split()
        rows.append([float(x) for x in parts])
    arr = np.array(rows)

    return DielectricData(
        energy=arr[:, 0],
        eps_imag_xx=arr[:, 1],
        eps_imag_yy=arr[:, 2],
        eps_imag_zz=arr[:, 3],
    )


# -----------------------------------------------------------------------------
# Auto-discovery within a VASP directory
# -----------------------------------------------------------------------------
def find_vasp_files(directory: str | Path) -> dict[str, Optional[Path]]:
    """Look inside a VASP calculation directory and find standard files.

    Returns a dict with keys: outcar, poscar, vasprun, eps_dat.
    Values are Path objects, or None if not found.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise NotADirectoryError(f"{directory} is not a directory")

    outcar = directory / "OUTCAR"
    poscar = directory / "POSCAR"
    vasprun = directory / "vasprun.xml"

    # Look for any *_eps_imag.dat file (sumo-optplot output)
    eps_candidates = list(directory.glob("*eps_imag*.dat"))
    eps_dat = eps_candidates[0] if eps_candidates else None

    return {
        "outcar": outcar if outcar.exists() else None,
        "poscar": poscar if poscar.exists() else None,
        "vasprun": vasprun if vasprun.exists() else None,
        "eps_dat": eps_dat,
    }
