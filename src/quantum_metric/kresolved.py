"""K-resolved quantum metric from VASP WAVEDER + EIGENVAL.

This module is independent of the library's existing optical-sum-rule path.
It computes g_µν(k) directly from the wavefunction-derivative matrix
elements VASP writes to WAVEDER when LOPTICS=.TRUE. is set.

VASP wiki convention (https://www.vasp.at/wiki/WAVEDER):
    CDER[μ, n, m, k] = <u_n,k | ∂/∂k_μ | u_m,k>   (in Å)

Only occupied↔unoccupied band pairs are stored; pair contributions
between two occupied or two unoccupied states are zero.

The quantum metric is then a direct contraction with no energy denominator:

    g_µν(k) = Re Σ_{n occ, m unocc} <u_n|∂_µ|u_m> <u_m|∂_ν|u_n>

Normalization
-------------
For consistency with the library's bulk SWM value, the per-electron metric is
normalized by the number of BOUND electrons n_bound (not NELECT). n_bound is
obtained from the same f-sum electron-counting used by QMetricCalculator
(intraband plasma frequency from OUTCAR). When the OUTCAR / dielectric data
needed for that count are unavailable, the code falls back to NELECT and warns.

Spin counting
-------------
VASP's CDER is stored per spin channel. For a collinear non-magnetic run
(no SOC) each stored band represents TWO electrons, so the bare single-channel
sum is a factor of 2 too small relative to the optical sum rule. For an SOC
run the spinor wavefunctions already include both spin components, so no extra
factor is needed. The spin-degeneracy factor is determined by reading the
LSORBIT flag from the OUTCAR:
    LSORBIT = .TRUE.  (SOC)        -> factor 1
    LSORBIT = .FALSE. (collinear)  -> factor 2

References:
    Souza, Wilkens & Martin, PRB 62, 1666 (2000).
    Marzari et al., Rev. Mod. Phys. 84, 1419 (2012).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np


__all__ = [
    "WavederData",
    "read_waveder",
    "read_eigenval",
    "compute_kresolved_metric",
]


# -----------------------------------------------------------------------------
# OUTCAR / SOC detection
# -----------------------------------------------------------------------------
def _soc_enabled(outcar_path):
    """Return True if LSORBIT=.TRUE. in the OUTCAR (spin-orbit coupling on).

    Handles both VASP output styles for the value:
      "LSORBIT = .TRUE."   (echoed INCAR tag)
      "LSORBIT =      T    spin-orbit coupling"  (parsed-parameter line)
    """
    try:
        with open(outcar_path) as f:
            for line in f:
                if "LSORBIT" in line:
                    rhs = line.split("=", 1)[1] if "=" in line else line
                    toks = rhs.upper().split()
                    if not toks:
                        continue
                    val = toks[0].strip(".")   # ".TRUE." -> "TRUE", "T" -> "T"
                    return val.startswith("T")
    except OSError:
        pass
    return False


# -----------------------------------------------------------------------------
# WAVEDER parser
# -----------------------------------------------------------------------------

@dataclass
class WavederData:
    """Parsed VASP WAVEDER file (LOPTICS output)."""
    nbands: int
    nb_cder: int
    nkpts: int
    nspin: int
    plasma_node_eV: float   # frequency at which Re ε(ω) = 0 (record 2)
    plasmon: np.ndarray     # 3x3 plasmon tensor (record 3)
    cder: np.ndarray        # (NKPTS, NSPIN, 3, NBANDS, NB_CDER) complex128


def _read_fortran_record(fh):
    head = fh.read(4)
    if not head or len(head) < 4:
        return None
    nbytes = int(np.frombuffer(head, dtype=np.int32, count=1)[0])
    payload = fh.read(nbytes)
    tail = int(np.frombuffer(fh.read(4), dtype=np.int32, count=1)[0])
    if nbytes != tail:
        raise IOError(f"WAVEDER record marker mismatch: head={nbytes} tail={tail}")
    return payload


def read_waveder(path):
    """Parse a VASP WAVEDER binary (Fortran-unformatted, complex64 storage)."""
    with open(path, "rb") as fh:
        ints = np.frombuffer(_read_fortran_record(fh), dtype=np.int32)
        nb_cder, nbands, nkpts, nspin = [int(x) for x in ints]

        plasma_node = float(
            np.frombuffer(_read_fortran_record(fh), dtype=np.float64)[0]
        )
        plasmon = np.frombuffer(_read_fortran_record(fh), dtype=np.float64)
        if plasmon.size == 9:
            plasmon = plasmon.reshape(3, 3)

        raw = _read_fortran_record(fh)
        flat = np.frombuffer(raw, dtype=np.complex64)
        expected = nb_cder * nbands * nkpts * nspin * 3
        if flat.size != expected:
            raise ValueError(
                f"WAVEDER payload size mismatch: got {flat.size}, expected {expected}"
            )
        cder = flat.reshape((nb_cder, nbands, nkpts, nspin, 3), order="F")
        cder = np.transpose(cder, (2, 3, 4, 1, 0)).astype(np.complex128)

    return WavederData(
        nbands=nbands, nb_cder=nb_cder, nkpts=nkpts, nspin=nspin,
        plasma_node_eV=plasma_node, plasmon=plasmon, cder=cder,
    )


# -----------------------------------------------------------------------------
# EIGENVAL parser (slim; only what we need here)
# -----------------------------------------------------------------------------

def read_eigenval(path):
    """Parse VASP EIGENVAL. Returns kpoints, kweights, energies, occupations."""
    with open(path) as f:
        lines = f.read().splitlines()
    nelect, nkpts, nbands = map(int, lines[5].split())

    kpts = np.zeros((nkpts, 3))
    wts = np.zeros(nkpts)
    energies = np.zeros((nkpts, nbands))
    occs = np.zeros((nkpts, nbands))

    idx = 7
    for ik in range(nkpts):
        parts = lines[idx].split()
        kpts[ik] = [float(x) for x in parts[:3]]
        wts[ik] = float(parts[3])
        idx += 1
        for ib in range(nbands):
            vals = lines[idx].split()
            energies[ik, ib] = float(vals[1])
            occs[ik, ib] = float(vals[2]) if len(vals) > 2 else 0.0
            idx += 1
        idx += 1

    return {
        "nelect": nelect,
        "nkpts": nkpts,
        "nbands": nbands,
        "kpoints": kpts,
        "kweights": wts,
        "energies": energies,
        "occupations": occs,
    }


# -----------------------------------------------------------------------------
# n_bound helper (uses the library's f-sum electron count)
# -----------------------------------------------------------------------------

def _get_n_bound(waveder_path, fallback_nelect):
    """Return n_bound for the calculation containing `waveder_path`.

    Runs the library's QMetricCalculator on the WAVEDER's parent directory to
    obtain the f-sum bound-electron count (consistent with the bulk SWM value).
    Falls back to NELECT (with a warning) if the required VASP files are
    missing or the calculator cannot run.
    """
    directory = Path(waveder_path).resolve().parent
    try:
        # Imported lazily to keep this module decoupled when n_bound isn't needed.
        from quantum_metric.calculator import QMetricCalculator
        calc = QMetricCalculator.from_directory(str(directory))
        res = calc.compute()
        return float(res.electrons.n_bound)
    except Exception as exc:
        warnings.warn(
            f"Could not compute n_bound from {directory} ({exc!r}); "
            f"falling back to NELECT={fallback_nelect}. "
            "Per-electron normalization will use NELECT instead of n_bound.",
            RuntimeWarning,
        )
        return float(fallback_nelect)


# -----------------------------------------------------------------------------
# K-resolved quantum metric
# -----------------------------------------------------------------------------

def compute_kresolved_metric(waveder_path, eigenval_path, per_electron=False,
                             spin_factor="auto", n_bound=None, outcar_path=None):
    """Compute g_µν(k) at every k-point.

    Formula (no 1/ΔE² weighting — CDER is the wavefunction derivative,
    not the velocity operator):

        g_µν(k) = Re Σ_{n occ, m unocc} <u_n|∂_µ|u_m> <u_m|∂_ν|u_n>

    Parameters
    ----------
    waveder_path : str
        Path to VASP's WAVEDER file (binary, from LOPTICS=.TRUE.).
    eigenval_path : str
        Path to VASP's EIGENVAL file.
    per_electron : bool, optional
        If True, divide g(k) by the number of BOUND electrons n_bound
        (for consistency with the library's bulk SWM value). The bound count
        is taken from `n_bound` if provided, otherwise computed from the
        WAVEDER's parent directory via QMetricCalculator, otherwise (if that
        fails) NELECT is used as a fallback.
    spin_factor : "auto" | float, optional
        Spin-degeneracy multiplier:
          - "auto" : read LSORBIT from the OUTCAR. SOC (LSORBIT=.TRUE.) -> 1.0;
                     collinear (no SOC) -> 2.0.
          - 1.0    : force no extra factor (SOC spinor runs).
          - 2.0    : force the collinear non-magnetic factor.
          - <float>: any custom multiplier.
        Default "auto".
    n_bound : float or None, optional
        Explicit bound-electron count to use for `per_electron` normalization.
        If None and `per_electron` is True, it is auto-computed (see above).
    outcar_path : str or None, optional
        Path to OUTCAR for SOC detection (spin_factor="auto"). If None, looks
        for an OUTCAR in the WAVEDER's parent directory.

    Returns
    -------
    kpoints : (NKPTS, 3) ndarray
        k-point fractional coordinates.
    g : (NKPTS, 3, 3) ndarray
        Quantum metric tensor at each k-point, in Å² (total if per_electron is
        False; per-bound-electron if per_electron is True).
    metadata : dict
        Energies, occupations, weights, NELECT, spin_factor, n_bound, soc, and
        the normalization actually used (norm_used) — useful downstream.
    """
    w = read_waveder(waveder_path)
    e = read_eigenval(eigenval_path)
    if w.nkpts != e["nkpts"]:
        raise ValueError(
            f"k-point mismatch: WAVEDER={w.nkpts}, EIGENVAL={e['nkpts']}"
        )

    # Reconcile band counts: EIGENVAL.nbands may exceed WAVEDER's band dims
    nb_min = min(e["nbands"], w.nbands, w.nb_cder)
    if not (e["nbands"] == w.nbands == w.nb_cder):
        # Silent truncation; caller can read metadata if curious
        e["energies"] = e["energies"][:, :nb_min]
        e["occupations"] = e["occupations"][:, :nb_min]
        cder = w.cder[:, :, :, :nb_min, :nb_min]
    else:
        cder = w.cder

    # ---- Resolve spin-degeneracy factor from LSORBIT ----------------------
    if outcar_path is None:
        outcar_path = Path(waveder_path).resolve().parent / "OUTCAR"
    soc = _soc_enabled(outcar_path)
    if spin_factor == "auto":
        # SOC (spinor) runs already include both spin components in CDER -> x1.
        # Collinear non-magnetic runs store one channel -> x2.
        sf = 1.0 if soc else 2.0
    else:
        sf = float(spin_factor)

    # ---- Contract CDER into g_µν(k) ---------------------------------------
    g = np.zeros((w.nkpts, 3, 3))
    for ik in range(w.nkpts):
        occ_mask = e["occupations"][ik] > 0.5
        unocc_mask = ~occ_mask
        if not occ_mask.any() or not unocc_mask.any():
            continue
        cd = cder[ik, 0]   # (3, NBANDS, NB_CDER)
        for mu in range(3):
            for nu in range(mu, 3):
                A = cd[mu][occ_mask][:, unocc_mask]
                B = cd[nu][occ_mask][:, unocc_mask]
                val = np.real(np.sum(A * np.conj(B)))
                g[ik, mu, nu] = val
                if mu != nu:
                    g[ik, nu, mu] = val

    # Apply spin-degeneracy factor (both-spins physical value)
    g = g * sf

    # ---- Per-electron normalization (by n_bound, consistent with SWM) -----
    norm_used = None
    if per_electron:
        if n_bound is None:
            n_bound = _get_n_bound(waveder_path, fallback_nelect=e["nelect"])
        g = g / n_bound
        norm_used = n_bound

    e["spin_factor"] = sf
    e["soc"] = soc
    e["n_bound"] = n_bound
    e["norm_used"] = norm_used
    return e["kpoints"], g, e
