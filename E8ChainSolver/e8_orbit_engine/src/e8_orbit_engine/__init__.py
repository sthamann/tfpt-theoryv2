"""
E8 Orbit Engine
===============

A data-driven engine to derive the quadratic damping function γ(n)
from nilpotent E8 orbits and validate the calibration-free consistency tests.
"""

__version__ = "0.2.0"

from .io import load_orbits, load_rge
from .chain import build_chain, validate_chain
from .chain_search import (
    beam_search_chains,
    save_chain_results,
    explain_chain,
    tokenize_label,
    build_hasse_graph
)
from .fit import (
    fit_gamma, compute_gamma, quad_and_d3, analyze_fit_quality,
    fit_gamma_log_model, fit_gamma_hyperbolic, compare_models
)
from .verify import (
    verify_lnD_cubic, 
    phi_third_diff_from_fit,
    phi_third_diff_from_data,
    verify_lnD_cubic_rolling,
    verify_rge_crossings,
    comprehensive_verification
)
from .report import generate_report

__all__ = [
    "load_orbits",
    "load_rge",
    "build_chain",
    "validate_chain",
    "beam_search_chains",
    "save_chain_results",
    "explain_chain",
    "fit_gamma",
    "compute_gamma",
    "quad_and_d3",
    "analyze_fit_quality",
    "fit_gamma_log_model",
    "fit_gamma_hyperbolic",
    "compare_models",
    "verify_lnD_cubic",
    "verify_lnD_cubic_rolling",
    "phi_third_diff_from_fit",
    "phi_third_diff_from_data",
    "verify_rge_crossings",
    "comprehensive_verification",
    "generate_report"
]
