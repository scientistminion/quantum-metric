"""Sumo-style band+metric plotting for the k-resolved feature.

Produces a 2-panel figure: band structure on top, g_µν(k) on the bottom,
sharing the high-symmetry x-axis.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np

from quantum_metric.kresolved import compute_kresolved_metric, read_eigenval


__all__ = ["plot_band_with_metric"]


# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------

SUMO_STYLE = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Bitstream Vera Serif",
                   "Times New Roman", "Times", "serif"],
    "mathtext.fontset": "dejavuserif",
    "text.usetex": False,
    "axes.labelsize": 15,
    "axes.titlesize": 16,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 13,
    "ytick.labelsize": 12,
    "xtick.major.size": 0,
    "xtick.minor.size": 0,
    "ytick.major.size": 4,
    "ytick.minor.size": 2,
    "ytick.direction": "in",
    "xtick.direction": "in",
    "axes.grid": False,
    "legend.frameon": False,
    "legend.fontsize": 11,
    "lines.linewidth": 1.3,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


# -----------------------------------------------------------------------------
# Helpers for line-mode KPOINTS / OUTCAR
# -----------------------------------------------------------------------------

def parse_kpath_kpoints(path):
    """Parse a VASP KPOINTS file (line-mode OR explicit/list mode).

    Returns (segment_breaks, label_seq) where:
        segment_breaks : list[int]
            k-point indices (0-based) where high-symmetry points occur.
            The first index is always 0; the last is NKPTS-1.
        label_seq : list[str]
            Corresponding labels (cleaned to plain text, e.g. "Gamma" not "$\\Gamma$").

    Supports:
      - Standard line-mode KPOINTS (segment endpoints with `! Label` comments)
      - Explicit/list mode KPOINTS (every k-point listed, labels only on some)
      - Labels in many styles: Gamma, $\\Gamma$, \\Gamma, GAMMA, Γ, G, Λ, etc.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]

    if len(lines) < 4:
        raise ValueError(f"KPOINTS file too short: {path}")

    # Header line 0: comment
    # Header line 1: count (either points-per-segment for line-mode, or NKPTS for explicit)
    # Header line 2: mode keyword (Line-mode, Reciprocal, Cartesian, Auto, ...)
    # Header line 3+: data
    try:
        header_count = int(lines[1].split()[0])
    except (ValueError, IndexError):
        header_count = 0

    mode_line = lines[2].strip().lower() if len(lines) > 2 else ""
    is_line_mode = mode_line.startswith("l") or "line" in mode_line

    if is_line_mode:
        return _parse_linemode(lines, header_count)
    else:
        return _parse_explicit(lines, header_count)


def _clean_label(raw):
    """Normalize a label string to a clean LaTeX-friendly name.

    Examples:
      "$\\Gamma$"  -> "Gamma"
      "\\Gamma"    -> "Gamma"
      "Γ"          -> "Gamma"
      "GAMMA"      -> "Gamma"
      "X_1"        -> "X_1"
      "$X_1$"      -> "X_1"
    """
    s = raw.strip()
    # Strip $...$ math delimiters
    s = s.strip("$")
    # Strip leading backslash (e.g. \Gamma -> Gamma)
    if s.startswith("\\"):
        s = s[1:]
    # Common Unicode → ASCII Greek
    unicode_to_name = {
        "Γ": "Gamma", "γ": "gamma",
        "Σ": "Sigma", "σ": "sigma",
        "Λ": "Lambda", "λ": "lambda",
        "Δ": "Delta", "δ": "delta",
        "Π": "Pi", "π": "pi",
        "Θ": "Theta", "θ": "theta",
        "Φ": "Phi", "φ": "phi",
        "Ω": "Omega", "ω": "omega",
    }
    for u, name in unicode_to_name.items():
        s = s.replace(u, name)
    # Title-case the canonical Greek names if all-caps
    canonical = {
        "GAMMA": "Gamma", "SIGMA": "Sigma", "LAMBDA": "Lambda",
        "DELTA": "Delta", "THETA": "Theta", "OMEGA": "Omega",
    }
    if s in canonical:
        s = canonical[s]
    return s


def _parse_linemode(lines, n_per_seg):
    """Parse standard VASP line-mode KPOINTS.

    Returns (segment_breaks, label_seq) such that the path can be split into
    n_per_seg-sized contiguous chunks in the EIGENVAL.
    """
    pts = []   # list of labels for each endpoint encountered
    for line in lines[4:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("!")
        coords = parts[0].split()
        label = _clean_label(parts[1]) if len(parts) > 1 else ""
        if len(coords) >= 3:
            pts.append(label)

    if not pts:
        return [0], [""]

    # Collapse adjacent-equal endpoints (segments share boundary points).
    # pts comes as [start1, end1, start2, end2, ...]; collapse pairs of
    # (end_i == start_{i+1}) into one label.
    labels = [pts[0]]
    for i in range(2, len(pts), 2):
        prev = pts[i - 1]
        cur = pts[i]
        labels.append(prev if prev == cur else f"{prev}|{cur}")
    labels.append(pts[-1])

    # Compute segment-break indices in the EIGENVAL (which has n_seg * n_per_seg points).
    n_seg = len(labels) - 1
    breaks = [seg * n_per_seg for seg in range(n_seg)] + [n_seg * n_per_seg - 1]
    return breaks, labels


def _parse_explicit(lines, declared_nkpts):
    """Parse explicit/list-mode KPOINTS.

    Each data line is a k-point with optional `! Label` or trailing-token label.
    Returns (segment_breaks, label_seq) by collecting indices of labeled points.
    """
    # Data starts at line 3 (after comment / count / mode-keyword)
    data_lines = lines[3:]
    labeled = []   # list of (index, label) pairs

    idx = 0
    for raw in data_lines:
        line = raw.strip()
        if not line:
            continue
        # Split off any `!`-style comment first
        if "!" in line:
            data_part, _, comment = line.partition("!")
            data_tokens = data_part.split()
            comment_tokens = comment.split()
        else:
            data_tokens = line.split()
            comment_tokens = []

        # Need at least 3 fractional coords to count as a k-point
        if len(data_tokens) < 3:
            continue
        # First three must parse as floats
        try:
            [float(t) for t in data_tokens[:3]]
        except ValueError:
            continue

        # Pull a label from:
        # 1. a `!` comment, OR
        # 2. trailing non-numeric tokens after the (coords, weight) on the data line
        label = ""
        if comment_tokens:
            label = _clean_label(" ".join(comment_tokens))
        else:
            # Tokens after coords + weight (i.e. position 4+) may be a label
            for tok in data_tokens[4:]:
                try:
                    float(tok)
                except ValueError:
                    # First non-numeric trailing token is the label
                    label = _clean_label(tok)
                    break

        if label:
            labeled.append((idx, label))
        idx += 1

    total = idx  # actual NKPTS from the file

    if not labeled:
        # No labels at all — return endpoints with empty labels
        return [0, max(total - 1, 0)], ["", ""]

    breaks = [p[0] for p in labeled]
    labels = [p[1] for p in labeled]

    # Ensure path includes start (0) and end (total-1)
    if breaks[0] != 0:
        breaks.insert(0, 0)
        labels.insert(0, "")
    if breaks[-1] != total - 1:
        breaks.append(total - 1)
        labels.append("")

    return breaks, labels


def reciprocal_from_outcar(outcar_path):
    """Read the reciprocal lattice (Cartesian, includes 2π) from OUTCAR."""
    with open(outcar_path) as f:
        text = f.read()
    idx = text.rfind("reciprocal lattice vectors")
    if idx < 0:
        return None
    vecs = []
    for line in text[idx:].splitlines()[1:6]:
        parts = line.split()
        if len(parts) >= 6:
            vecs.append([float(parts[3]), float(parts[4]), float(parts[5])])
        if len(vecs) == 3:
            break
    if len(vecs) != 3:
        return None
    return np.array(vecs) * 2 * np.pi


def _format_label(label):
    """Wrap a cleaned label (e.g. 'Gamma', 'X_1', 'L|W') in LaTeX math.

    Labels from parse_kpath_kpoints come already normalized to plain ASCII
    Greek names ("Gamma" not "Γ" or "$\\Gamma$"); this function just adds
    backslashes for Greek and wraps in $...$.
    """
    if not label:
        return ""
    out = label
    for greek in ("Gamma", "Sigma", "Lambda", "Delta", "Pi", "Theta",
                  "Phi", "Omega"):
        out = out.replace(greek, fr"\{greek}")
    return rf"${out}$"


def _build_path_distance(kpts, segment_breaks, recip_matrix):
    """Cumulative Cartesian k-distance, resetting at segment-boundary breaks.

    Parameters
    ----------
    kpts : (N, 3) ndarray
        Fractional k-coordinates.
    segment_breaks : sequence of ints
        Indices where new segments begin (the boundary k-points). The distance
        is held flat across each break to indicate a discontinuity.
    recip_matrix : (3, 3) ndarray or None
        Reciprocal lattice vectors (Cartesian, includes 2π). If None, fractional
        Euclidean distance is used.
    """
    n_kpts = len(kpts)
    dist = np.zeros(n_kpts)
    kcart = kpts @ recip_matrix if recip_matrix is not None else kpts
    # Indices that start a new segment (excluding the very first point)
    seg_starts = set(b for b in segment_breaks if 0 < b < n_kpts - 1)
    for i in range(1, n_kpts):
        if i in seg_starts:
            # Hold distance flat at segment boundary; new segment starts here.
            dist[i] = dist[i - 1]
        else:
            dist[i] = dist[i - 1] + np.linalg.norm(kcart[i] - kcart[i - 1])
    return dist


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

def plot_band_with_metric(
    waveder, eigenval, kpoints_file, output,
    outcar=None,
    per_electron=False,
    fermi=None, ymin=-5, ymax=5,
    title=None,
    width=8.0, height=8.0,
    height_ratio=(2.5, 1.0),
    band_color="#1f77b4",
    trace_only=False,
):
    """Two-panel plot: band structure (top) and g_µν(k) (bottom).

    Parameters
    ----------
    waveder, eigenval, kpoints_file : file paths
    outcar : str or None
        If given, used to read the reciprocal lattice for correct k-distances.
    per_electron : bool
        Plot g(k)/NELECT instead of total g(k).
    trace_only : bool
        Plot only Tr[g] in the bottom panel (cleanest for cubic crystals).
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update(SUMO_STYLE)

    # Compute metric and re-read full EIGENVAL for the band plot
    kpts, gk, e = compute_kresolved_metric(waveder, eigenval, per_electron=per_electron)
    # Re-read EIGENVAL (untruncated) for the band-structure plot
    e_full = read_eigenval(eigenval)

    segment_breaks, label_seq = parse_kpath_kpoints(kpoints_file)
    recip = reciprocal_from_outcar(outcar) if outcar else None
    dist = _build_path_distance(kpts, segment_breaks, recip)

    if fermi is None:
        fermi = e_full["energies"][e_full["occupations"] > 0.5].max()

    fig, axes = plt.subplots(
        2, 1, figsize=(width, height), sharex=True,
        gridspec_kw={"height_ratios": height_ratio, "hspace": 0.05},
    )
    ax_bs, ax_m = axes

    # Bands
    bands = e_full["energies"] - fermi
    for ib in range(bands.shape[1]):
        ax_bs.plot(dist, bands[:, ib], color=band_color, lw=1.4, alpha=0.95)
    ax_bs.axhline(0, color="black", ls="--", lw=0.8, alpha=0.6)
    ax_bs.set_ylabel(r"Energy $-$ E$_{\rm F}$ (eV)")
    ax_bs.set_ylim(ymin, ymax)
    if title:
        ax_bs.set_title(title, pad=10)

    # Metric panel
    if trace_only:
        trace = gk[:, 0, 0] + gk[:, 1, 1] + gk[:, 2, 2]
        ax_m.plot(dist, trace, color="#d62728", lw=1.8, label=r"Tr$\,[g]$")
    else:
        ax_m.plot(dist, gk[:, 0, 0], color="#1f77b4", lw=1.4, label=r"$g_{xx}$")
        ax_m.plot(dist, gk[:, 1, 1], color="#2ca02c", lw=1.4, ls="--", label=r"$g_{yy}$")
        ax_m.plot(dist, gk[:, 2, 2], color="#9467bd", lw=1.4, ls=":", label=r"$g_{zz}$")
        trace = gk[:, 0, 0] + gk[:, 1, 1] + gk[:, 2, 2]
        ax_m.plot(dist, trace, color="#d62728", lw=1.8, label=r"Tr$\,[g]$")

    ylabel = r"$g_{\mu\nu}(\mathbf{k})$  ($\rm\AA^{2}$)" if per_electron \
        else r"$g_{\mu\nu}(\mathbf{k})$  ($\rm\AA^{2}$)"
    ax_m.set_ylabel(ylabel)
    ax_m.set_xlabel("Wavevector")
    ax_m.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5),
        ncol=1, handlelength=1.8, borderpad=0.4, labelspacing=0.6,
    )
    ax_m.set_ylim(bottom=0)

    # High-symmetry tick positions — directly from the parsed breaks
    tick_positions = [dist[i] for i in segment_breaks]
    tick_labels = [_format_label(l) for l in label_seq]

    for ax in (ax_bs, ax_m):
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels if ax is ax_m else [])
        for x in tick_positions[1:-1]:
            ax.axvline(x, color="black", lw=0.7, alpha=0.55)
        ax.set_xlim(dist[0], dist[-1])
        ax.tick_params(which="both", direction="in", top=True, right=True)

    fig.align_ylabels(axes)
    fig.savefig(output)
    plt.close(fig)
    return output
