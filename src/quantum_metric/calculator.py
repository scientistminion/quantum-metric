"""
High-level Calculator that orchestrates the full pipeline.

Typical usage:
    >>> calc = QMetricCalculator.from_directory("./MoS2")
    >>> result = calc.compute()
    >>> print(result)

Override files:
    >>> calc = QMetricCalculator(
    ...     outcar="./a/OUTCAR",
    ...     poscar="./b/POSCAR",
    ...     dielectric="./c/MoS2_eps_imag.dat",
    ... )
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from quantum_metric.electrons import (
    ElectronCount,
    compute_n_itinerant_fsum,
    compute_n_itinerant_kai,
)
from quantum_metric.io import (
    DielectricData,
    OutcarData,
    PoscarData,
    find_vasp_files,
    read_dielectric,
    read_outcar,
    read_poscar,
)
from quantum_metric.metric import (
    DEFAULT_PREFACTOR,
    QuantumMetricResult,
    compute_quantum_metric,
)
from quantum_metric.optics import (
    OpticalIntegralsAllDirections,
    compute_optical_integrals,
)


@dataclass
class QMetricResult:
    """Full result of a quantum-metric calculation for one material."""

    # Inputs summary
    material: str
    method: str                 # "kai" or "fsum"
    volume: float               # Ang^3
    nelect: float
    natoms: int
    a_len: float                # Ang

    # Plasma / sumrule
    plasma_intra_xx: float      # eV^2
    plasma_inter_xx: float      # eV^2
    sumrule: float              # eV^2

    # Optical integrals
    optical: OpticalIntegralsAllDirections

    # Electron counts
    electrons: ElectronCount

    # Quantum metric
    metric: QuantumMetricResult

    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        """Flatten the result into a dict suitable for DataFrame/JSON/TSV."""
        d = {
            "Material": self.material,
            "Method": self.method,
            "Volume_Ang3": self.volume,
            "NELECT": self.nelect,
            "NAtoms": self.natoms,
            "a_len_Ang": self.a_len,
            "PlasmaFreqSquared_Intraband_xx_eV2": self.plasma_intra_xx,
            "PlasmaFreqSquared_Interband_xx_eV2": self.plasma_inter_xx,
            "Sumrule_eV2": self.sumrule,
            # Optical integrals (xx)
            "omega_p_squared_xx": self.optical.xx.omega_p_squared,
            "I_xx": self.optical.xx.I,
            "sigma_int_xx": self.optical.xx.sigma_int,
            "wsigma_xx": self.optical.xx.wsigma,
            "sigma_xx_over_wsquare": self.optical.xx.sigma_over_wsquare,
            # Electrons
            "Kai": self.electrons.kai,
            "N_itinerant": self.electrons.n_itinerant,
            "N_bound": self.electrons.n_bound,
            "N_itinerant_per_atom": self.electrons.n_itinerant_per_atom,
            "N_bound_per_atom": self.electrons.n_bound_per_atom,
            "bound_electron_density": self.electrons.bound_electron_density,
            "fsum_residual": self.electrons.fsum_residual,
            # Metric
            "sqrtG_over_A_xx": self.metric.sqrtG_over_A_xx,
            "prefactor": self.metric.prefactor,
        }
        # Anisotropic fields (only present if dielectric had yy, zz)
        if self.optical.yy is not None:
            d["I_yy"] = self.optical.yy.I
            d["sigma_int_yy"] = self.optical.yy.sigma_int
            d["sqrtG_over_A_yy"] = self.metric.sqrtG_over_A_yy
        if self.optical.zz is not None:
            d["I_zz"] = self.optical.zz.I
            d["sigma_int_zz"] = self.optical.zz.sigma_int
            d["sqrtG_over_A_zz"] = self.metric.sqrtG_over_A_zz
        return d


# -----------------------------------------------------------------------------
class QMetricCalculator:
    """Orchestrator for the quantum-metric pipeline on a single VASP run."""

    def __init__(
        self,
        *,
        outcar: str | Path,
        poscar: str | Path,
        dielectric: str | Path,
        dielectric_source: str = "auto",
        material: Optional[str] = None,
    ):
        self.outcar_path = Path(outcar)
        self.poscar_path = Path(poscar)
        self.dielectric_path = Path(dielectric)
        self.dielectric_source = dielectric_source
        self.material = material or self.poscar_path.parent.name or "unknown"

        # Lazy-load on first access
        self._outcar: Optional[OutcarData] = None
        self._poscar: Optional[PoscarData] = None
        self._dielectric: Optional[DielectricData] = None

    # ------------------------------------------------------------------
    @classmethod
    def from_directory(
        cls,
        directory: str | Path,
        *,
        prefer_vasprun: bool = True,
    ) -> "QMetricCalculator":
        """Auto-discover files inside a VASP calculation directory.

        Order of preference for the dielectric function:
          1. If prefer_vasprun and vasprun.xml exists -> use it (no sumo needed)
          2. Otherwise fall back to *_eps_imag.dat
        """
        directory = Path(directory).resolve()
        files = find_vasp_files(directory)

        if files["outcar"] is None:
            raise FileNotFoundError(f"No OUTCAR found in {directory}")
        if files["poscar"] is None:
            raise FileNotFoundError(f"No POSCAR found in {directory}")

        if prefer_vasprun and files["vasprun"] is not None:
            dielectric = files["vasprun"]
            source = "vasprun"
        elif files["eps_dat"] is not None:
            dielectric = files["eps_dat"]
            source = "dat"
        elif files["vasprun"] is not None:
            dielectric = files["vasprun"]
            source = "vasprun"
        else:
            raise FileNotFoundError(
                f"No dielectric data found in {directory}. "
                "Expected either vasprun.xml or *_eps_imag.dat"
            )

        return cls(
            outcar=files["outcar"],
            poscar=files["poscar"],
            dielectric=dielectric,
            dielectric_source=source,
            material=directory.name,
        )

    # ------------------------------------------------------------------
    @property
    def outcar(self) -> OutcarData:
        if self._outcar is None:
            self._outcar = read_outcar(self.outcar_path)
        return self._outcar

    @property
    def poscar(self) -> PoscarData:
        if self._poscar is None:
            self._poscar = read_poscar(self.poscar_path)
        return self._poscar

    @property
    def dielectric(self) -> DielectricData:
        if self._dielectric is None:
            self._dielectric = read_dielectric(
                self.dielectric_path, source=self.dielectric_source
            )
        return self._dielectric

    # ------------------------------------------------------------------
    def compute(
        self,
        *,
        method: str = "kai",
        prefactor: float = DEFAULT_PREFACTOR,
        e_min: float = 0.0,
        e_max: Optional[float] = None,
    ) -> QMetricResult:
        """Run the full pipeline.

        Parameters
        ----------
        method : 'kai' | 'fsum'
            How to compute N_itinerant / N_bound.
        prefactor : float
            Unit-conversion constant for the metric (default 0.0694).
        e_min, e_max : float
            Integration window on eps_imag (eV).
        """
        o = self.outcar
        p = self.poscar
        d = self.dielectric

        optical = compute_optical_integrals(d, e_min=e_min, e_max=e_max)

        if method == "kai":
            electrons = compute_n_itinerant_kai(
                plasma_intra=o.plasma_intra_xx,
                plasma_inter=o.plasma_inter_xx,
                nelect=o.nelect,
                volume=o.volume,
                natoms=o.natoms,
            )
        elif method == "fsum":
            electrons = compute_n_itinerant_fsum(
                plasma_intra_ev2=o.plasma_intra_xx,
                plasma_inter_ev2=o.plasma_inter_xx,
                nelect=o.nelect,
                volume_ang3=o.volume,
                natoms=o.natoms,
            )
        else:
            raise ValueError(f"Unknown method: {method!r}. Use 'kai' or 'fsum'.")

        metric = compute_quantum_metric(
            I_xx=optical.xx.I,
            I_yy=optical.yy.I if optical.yy is not None else None,
            I_zz=optical.zz.I if optical.zz is not None else None,
            bound_electron_density=electrons.bound_electron_density,
            prefactor=prefactor,
        )

        return QMetricResult(
            material=self.material,
            method=method,
            volume=o.volume,
            nelect=o.nelect,
            natoms=o.natoms,
            a_len=p.a_len,
            plasma_intra_xx=o.plasma_intra_xx,
            plasma_inter_xx=o.plasma_inter_xx,
            sumrule=o.sumrule,
            optical=optical,
            electrons=electrons,
            metric=metric,
        )
