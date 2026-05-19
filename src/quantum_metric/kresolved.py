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

References:
    Souza, Wilkens & Martin, PRB 62, 1666 (2000).
    Marzari et al., Rev. Mod. Phys. 84, 1419 (2012).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


__all__ = [
    "WavederData",
    "read_waveder",
    "read_eigenval",
    "compute_kresolved_metric",
]


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
# K-resolved quantum metric
# -----------------------------------------------------------------------------

def compute_kresolved_metric(waveder_path, eigenval_path, per_electron=False):
    """Compute g_µν(k) at every k-point.

    Formula (no 1/ΔE² weighting — CDER is the wavefunction derivative,
    not the velocity operator):

        g_µν(k) = Re Σ_{n occ, m unocc} <u_n|∂_µ|u_m> <u_m|∂_ν|u_n>

    Parameters
    ----------
    waveder_path : str
        Path to VASP's WAVEDER file (binary, from LOPTICS=.TRUE.)
    eigenval_path : str
        Path to VASP's EIGENVAL file
    per_electron : bool, optional
        If True, divide g(k) by NELECT.

    Returns
    -------
    kpoints : (NKPTS, 3) ndarray
        k-point fractional coordinates
    g : (NKPTS, 3, 3) ndarray
        Quantum metric tensor at each k-point, in Å²
    metadata : dict
        Energies, occupations, weights, NELECT — useful for downstream plots
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

    if per_electron:
        g = g / e["nelect"]

    return e["kpoints"], g, e
