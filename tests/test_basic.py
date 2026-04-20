"""Basic smoke tests. Uses tmp_path to create synthetic VASP outputs."""

from pathlib import Path

import numpy as np
import pytest

from quantum_metric import (
    QMetricCalculator,
    compute_kai,
    compute_n_itinerant_fsum,
    compute_n_itinerant_kai,
    read_dielectric,
    read_outcar,
    read_poscar,
)


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
def test_kai_ratio():
    k = compute_kai(10.0, 40.0)
    assert k == pytest.approx(0.2, rel=1e-12)


def test_kai_electron_count():
    e = compute_n_itinerant_kai(plasma_intra=10.0, plasma_inter=40.0,
                                nelect=26.0, volume=100.0, natoms=3)
    assert e.kai == pytest.approx(0.2)
    assert e.n_itinerant == pytest.approx(5.2)
    assert e.n_bound == pytest.approx(20.8)


def test_fsum_electron_count():
    # Check both itinerant and bound are computed from their own plasma freqs
    e = compute_n_itinerant_fsum(
        plasma_intra_ev2=10.0,
        plasma_inter_ev2=40.0,
        nelect=26.0,
        volume_ang3=100.0,
        natoms=3,
    )
    assert e.n_itinerant > 0
    assert e.n_bound > 0
    assert e.method == "fsum"
    # The ratio of the two should match the ratio of plasma frequencies squared
    assert e.n_bound / e.n_itinerant == pytest.approx(40.0 / 10.0, rel=1e-10)
    # The residual should be the difference from NELECT
    assert e.fsum_residual == pytest.approx(26.0 - (e.n_itinerant + e.n_bound))

# --- End-to-end pipeline ---------------------------------------------------
def test_calculator_from_directory_kai(vasp_dir):
    calc = QMetricCalculator.from_directory(vasp_dir)
    r = calc.compute(method="kai")

    assert r.material == vasp_dir.name
    assert r.method == "kai"
    assert r.volume == 100.0
    assert r.nelect == 26.0
    assert r.electrons.kai == pytest.approx(0.2)
    assert r.metric.sqrtG_over_A_xx > 0
    assert r.metric.sqrtG_over_A_yy is not None
    assert r.metric.sqrtG_over_A_zz is not None


def test_calculator_from_directory_fsum(vasp_dir):
    calc = QMetricCalculator.from_directory(vasp_dir)
    r = calc.compute(method="fsum")
    assert r.method == "fsum"
    assert r.electrons.kai is None   # not computed in fsum mode
    assert r.metric.sqrtG_over_A_xx > 0


def test_to_dict_flat(vasp_dir):
    calc = QMetricCalculator.from_directory(vasp_dir)
    r = calc.compute()
    d = r.to_dict()
    # must be flat (no nested dicts)
    for v in d.values():
        assert not isinstance(v, dict)
    assert "sqrtG_over_A_xx" in d
    assert "I_yy" in d   # anisotropic present


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
