"""
quantum-metric: Compute the quantum metric and optical quantities from VASP output.

High-level API:
    >>> from quantum_metric import QMetricCalculator
    >>> calc = QMetricCalculator.from_directory("./my_vasp_run")
    >>> result = calc.compute()
    >>> print(result.sqrtG_over_A)

Or use individual modules for fine-grained control.
"""

from quantum_metric._version import __version__
from quantum_metric.calculator import QMetricCalculator, QMetricResult
from quantum_metric.io import read_outcar, read_poscar, read_dielectric
from quantum_metric.optics import compute_optical_integrals
from quantum_metric.electrons import (
    compute_kai,
    compute_n_itinerant_kai,
    compute_n_itinerant_fsum,
)
from quantum_metric.metric import compute_quantum_metric

__all__ = [
    "__version__",
    "QMetricCalculator",
    "QMetricResult",
    "read_outcar",
    "read_poscar",
    "read_dielectric",
    "compute_optical_integrals",
    "compute_kai",
    "compute_n_itinerant_kai",
    "compute_n_itinerant_fsum",
    "compute_quantum_metric",
]
