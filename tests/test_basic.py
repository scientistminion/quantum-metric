"""Basic smoke tests. Uses tmp_path to create synthetic VASP outputs."""

from pathlib import Path

import numpy as np
import pytest

from quantum_metric import (
    QMetricCalculator,
    compute_electron_count,
    read_dielectric,
    read_outcar,
    read_poscar,
)
from quantum_metric.electrons import PREFACTOR_N


# --- Fixture builders -------------------------------------------------------
SYNTHETIC_OUTCAR = """ running on some machine

    NIONS = 3
    NELECT =       26.0000000

   volume of cell :      100.0000

 plasma frequency squared (from intraband transitions) (eV^2)
  --------------------------------------------------------------------------------------
      10.000       0.000       0.000
       0.000      10.000       0.000
       0.000       0.000      10.000

 plasma frequency squared (from interband transitions, int dw w eps(2)(w)) (eV^2)
  --------------------------------------------------------------------------------------
      40.000       0.000       0.000      -0.000       0.000       0.000
       0.000       0.000      40.000       0.000       0.000      -0.000
       0.000      -0.000       0.000       0.000      40.000       0.000
  --------------------------------------------------------------------------------------
 sumrule: sum of plasma frequencies squared should yield (valence only):      50.000
"""

SYNTHETIC_POSCAR = """Test
   1.0
     3.0  0.0  0.0
     0.0  4.0  0.0
     0.0  0.0  5.0
   H
   1
Direct
   0.0  0.0  0.0
"""


def _make_dielectric_dat(path: Path, anisotropic: bool = True):
    """Make a simple eps_imag.dat with a Lorentzian peak."""
    e = np.linspace(0.01, 20.0, 500)
    # Lorentzian around 5 eV
    eps_xx = 10.0 / ((e - 5.0) ** 2 + 0.5)
    if anisotropic:
        eps_yy = 8.0 / ((e - 6.0) ** 2 + 0.5)
        eps_zz = 12.0 / ((e - 4.0) ** 2 + 0.5)
        data = np.column_stack([e, eps_xx, eps_yy, eps_zz])
    else:
        data = np.column_stack([e, eps_xx])
    np.savetxt(path, data)


@pytest.fixture
def vasp_dir(tmp_path):
    (tmp_path / "OUTCAR").write_text(SYNTHETIC_OUTCAR)
    (tmp_path / "POSCAR").write_text(SYNTHETIC_POSCAR)
    _make_dielectric_dat(tmp_path / "test_eps_imag.dat", anisotropic=True)
    return tmp_path


# --- I/O tests --------------------------------------------------------------
def test_read_outcar(vasp_dir):
    d = read_outcar(vasp_dir / "OUTCAR")
    assert d.plasma_intra_xx == 10.0
    assert d.plasma_inter_xx == 40.0
    assert d.sumrule == 50.0
    assert d.volume == 100.0
    assert d.nelect == 26.0
    assert d.natoms == 3


def test_read_poscar(vasp_dir):
    p = read_poscar(vasp_dir / "POSCAR")
    assert p.a_len == 3.0
    np.testing.assert_allclose(p.lattice_lengths, [3.0, 4.0, 5.0])


def test_read_dielectric_dat(vasp_dir):
    d = read_dielectric(vasp_dir / "test_eps_imag.dat")
    assert d.has_anisotropic
    assert d.energy.shape == (500,)


# --- Electron-count tests ---------------------------------------------------
def test_prefactor_value():
    """The hydrogen-relations prefactor 1/(16π a_B³ E_0²) should be ~7.263e-4 Å⁻³ eV⁻²."""
    assert PREFACTOR_N == pytest.approx(7.263e-4, rel=1e-3)


def test_electron_count_synthetic():
    """Check basic f-sum electron count on synthetic data."""
    e = compute_electron_count(
        plasma_intra_ev2=10.0,
        nelect=26.0,
        volume_ang3=100.0,
        natoms=3,
        sumrule_ev2=50.0,
    )
    # n_itinerant = PREFACTOR_N × X_intra × V
    expected_n_it = PREFACTOR_N * 10.0 * 100.0
    assert e.n_itinerant == pytest.approx(expected_n_it, rel=1e-10)
    assert e.n_bound == pytest.approx(26.0 - expected_n_it, rel=1e-10)
    assert e.n_itinerant_per_atom == pytest.approx(e.n_itinerant / 3, rel=1e-10)
    assert e.n_bound_per_atom == pytest.approx(e.n_bound / 3, rel=1e-10)
    assert e.bound_electron_density == pytest.approx(e.n_bound / 100.0, rel=1e-10)
    assert e.itinerant_electron_density == pytest.approx(e.n_itinerant / 100.0, rel=1e-10)
    # Sumrule check: PREFACTOR_N × sumrule × V
    assert e.sumrule_check == pytest.approx(PREFACTOR_N * 50.0 * 100.0, rel=1e-10)


def test_na_bcc_sanity():
    """Sodium bcc should give ~1 itinerant electron per atom (the 3s¹ conduction electron).

    Reference data from a real VASP calculation:
      X_intra = 37.394 eV²,   X_sumrule = 245.701 eV²,
      V = 39.28 Å³,           NELECT = 7  (Na_pv: 2p⁶ 3s¹),  NIONS = 1
    """
    e = compute_electron_count(
        plasma_intra_ev2=37.394,
        nelect=7.0,
        volume_ang3=39.28,
        natoms=1,
        sumrule_ev2=245.701,
    )
    # ~1 itinerant electron per atom (the 3s¹ valence)
    assert e.n_itinerant == pytest.approx(1.067, rel=2e-2)
    # Bound = semicore 2p⁶ minus tiny f-sum truncation correction
    assert e.n_bound == pytest.approx(5.933, rel=2e-2)
    # Sumrule should imply NELECT ≈ 7 — a clean validation of the prefactor
    assert e.sumrule_check == pytest.approx(7.0, rel=2e-2)


# --- End-to-end pipeline ---------------------------------------------------
def test_calculator_from_directory(vasp_dir):
    calc = QMetricCalculator.from_directory(vasp_dir)
    r = calc.compute()

    assert r.material == vasp_dir.name
    assert r.volume == 100.0
    assert r.nelect == 26.0
    assert r.electrons.n_itinerant > 0
    assert r.electrons.n_bound > 0
    assert r.metric.sqrtG_over_A_xx > 0
    assert r.metric.sqrtG_over_A_yy is not None
    assert r.metric.sqrtG_over_A_zz is not None


def test_to_dict_flat(vasp_dir):
    calc = QMetricCalculator.from_directory(vasp_dir)
    r = calc.compute()
    d = r.to_dict()
    # must be flat (no nested dicts)
    for v in d.values():
        assert not isinstance(v, dict)
    # required scalar fields
    assert "sqrtG_over_A_xx" in d
    assert "I_yy" in d   # anisotropic present
    assert "N_itinerant" in d
    assert "N_bound" in d
    assert "sumrule_check_NELECT" in d
    # method/Kai must NOT be present anymore
    assert "Method" not in d
    assert "Kai" not in d


def test_override_files(vasp_dir, tmp_path):
    """Test overriding individual files."""
    # Put OUTCAR in a different place
    other = tmp_path / "elsewhere"
    other.mkdir()
    (other / "OUTCAR").write_text(SYNTHETIC_OUTCAR)

    calc = QMetricCalculator(
        outcar=other / "OUTCAR",
        poscar=vasp_dir / "POSCAR",
        dielectric=vasp_dir / "test_eps_imag.dat",
    )
    r = calc.compute()
    assert r.volume == 100.0
