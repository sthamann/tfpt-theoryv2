import sys
import unittest
from pathlib import Path

import numpy as np

# Ensure `tfpt_suite` package is importable when tests are run from workspace root.
TFPT_SUITE_DIR = Path(__file__).resolve().parents[1]
if str(TFPT_SUITE_DIR) not in sys.path:
    sys.path.insert(0, str(TFPT_SUITE_DIR))

from tfpt_suite.flavor_textures import C_delta, circulant, left_unitary_from_yukawa, scale_y_star_to_match_sigma_max, theta_of_delta  # noqa: E402
from tfpt_suite.modules.pmns_full_pipeline import _pmns_best_convention, _pmns_from_ye_and_kappa, _pmns_pdg  # noqa: E402


class TestFlavorTextures(unittest.TestCase):
    """
    Unit tests for the flavor-texture math utilities.

    These helpers are used by the CKM/PMNS pipelines; failures here typically mean a
    convention drift (phase mode, hermiticity) or a numerical linear algebra regression.
    """

    def test_circulant_layout(self):
        M = circulant(1, 2, 3)
        want = np.array([[1, 2, 3], [3, 1, 2], [2, 3, 1]], dtype=complex)
        self.assertTrue(np.allclose(M, want))

    def test_C_delta_is_hermitian(self):
        delta = 0.61
        C = C_delta(delta, phase_mode="2pi_delta")
        self.assertTrue(np.allclose(C, C.conj().T))

    def test_theta_modes(self):
        delta = 0.5
        self.assertAlmostEqual(theta_of_delta(delta, phase_mode="2pi_delta"), np.pi, places=12)
        self.assertAlmostEqual(theta_of_delta(delta, phase_mode="delta_rad"), 0.5, places=12)
        self.assertAlmostEqual(theta_of_delta(delta, phase_mode="koide_pi_over_12"), np.pi / 12.0, places=12)

    def test_scale_y_star_matches_sigma_max(self):
        base = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], dtype=complex)
        ystar = scale_y_star_to_match_sigma_max(target_y3=0.9, base=base)
        smax = float(np.max(np.linalg.svd(ystar * base, compute_uv=False)))
        self.assertAlmostEqual(smax, 0.9, places=12)

    def test_left_unitary_is_unitary(self):
        Y = np.array([[1.0, 0.2j, 0.0], [0.1, 0.9, 0.3j], [0.0, 0.2, 0.7]], dtype=complex)
        U = left_unitary_from_yukawa(Y)
        self.assertTrue(np.allclose(U.conj().T @ U, np.eye(3), atol=1e-12))

    def test_pmns_takagi_ordering_is_paired_with_masses(self):
        """
        Regression test: PMNS extraction must not silently mislabel columns due to SVD ordering.

        For a symmetric Majorana mν built from a known PDG PMNS matrix, `_pmns_from_ye_and_kappa`
        must keep the singular values paired with the corresponding columns (after reordering),
        so that a permutation scan can recover the expected small θ13.
        """
        Ye = np.eye(3, dtype=complex)  # Ue=I
        v = 246.0
        # A realistic-ish PMNS point (deg) with small θ13
        th12 = np.deg2rad(34.3225)
        th13 = np.deg2rad(8.7437)
        th23 = np.deg2rad(45.5078)
        dcp = np.deg2rad(90.0)
        U = _pmns_pdg(th12, th13, th23, dcp)

        m_eV = np.array([0.0, 0.008602325267042627, 0.05], dtype=float)
        m_GeV = m_eV * 1.0e-9
        mnu_GeV = (U @ np.diag(m_GeV) @ U.T).astype(complex)
        kappa = (mnu_GeV / (v**2)).astype(complex)

        U_rec, masses_eV = _pmns_from_ye_and_kappa(Ye=Ye, kappa=kappa, v_ev_GeV=v)
        fit = _pmns_best_convention(U_pmns=U_rec, mnu_eV=masses_eV, refs_by_ordering={})
        ang = fit["angles"]

        self.assertLess(abs(ang.theta13_deg - np.degrees(th13)), 1e-3)


if __name__ == "__main__":
    unittest.main()

