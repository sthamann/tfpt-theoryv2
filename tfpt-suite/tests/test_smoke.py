import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path


# Ensure `tfpt_suite` package is importable when tests are run from workspace root.
TFPT_SUITE_DIR = Path(__file__).resolve().parents[1]
if str(TFPT_SUITE_DIR) not in sys.path:
    sys.path.insert(0, str(TFPT_SUITE_DIR))

EXPECTED_DEFECT_MULTIPLICITY = 5
EXPECTED_DEFECT_CLASSES = {"bound_color", "bound_weak", "separated_color", "separated_weak", "seam_coupled_color_weak"}
EXPECTED_ALPHA_BRIDGE_CHECKS = {
    "alpha_bridge_leptons_1loop_explicit",
    "alpha_bridge_hadron_policy_declared",
    "alpha_bridge_EW_decoupling_included",
}
B0_TEST_POINTS = 96
FINITE_PIECE_MIN = 1e-6
TEST_MASS_GEV = 2.0
TEST_MU_GEV = 3.0
TEST_V_EV_GEV = 246.0
TEST_OMEGA_A_H2 = 0.12


from tfpt_suite.config import SuiteConfig  # noqa: E402
from tfpt_suite.conventions import GaugeConvention, convert_g1, g1_gut_from_gY, g1_gut_over_gY, gY_from_g1_gut, g_from_alpha  # noqa: E402
from tfpt_suite.matching import (  # noqa: E402
    a0_msbar,
    b0_msbar,
    delta_lambda_h_1loop_buttazzo,
    delta_y_t_ew_hempfling,
    match_gauge,
)
from tfpt_suite.modules.registry import get_module_registry  # noqa: E402
from tfpt_suite.operator_spec_builder import generate_effective_action_r2_operator_spec  # noqa: E402
from tfpt_suite.rg_authority import load_rg_authority_policy, mz_scale_GeV  # noqa: E402
from tfpt_suite.rge_sm import run_sm_gauge_only_2loop_thresholds  # noqa: E402


class TestSmoke(unittest.TestCase):
    """
    Smoke tests for the suite wiring.

    These tests:
    - ensure module discovery is stable (registry contains expected IDs)
    - run a small subset of modules end-to-end into a temporary output directory
    - explicitly disable plotting (headless / CI-friendly)
    """

    def test_registry_has_expected_modules(self):
        reg = get_module_registry()
        for module_id in [
            "core_invariants",
            "two_loop_rg_fingerprints",
            "alpha_precision_audit",
            "effective_action_r2",
            "primordial_spectrum_builder",
            "cosmo_threshold_history",
            "k_calibration",
            "axion_fa_derivation",
            "msbar_matching_map",
            "below_mt_eft_cascade",
            "stability_unitarity_audit",
            "aps_eta_gluing",
            "discrete_consistency_uniqueness",
            "pmns_mechanism_bridge",
            "likelihood_engine",
            "torsion_falsifiability_snr",
        ]:
            self.assertIn(module_id, reg)

    def test_g1_gut_over_gY_ratio(self):
        gY = 0.357
        g1_gut = g1_gut_from_gY(gY)
        self.assertLess(abs((g1_gut / gY) - g1_gut_over_gY()), 1e-12)
        self.assertLess(abs(gY_from_g1_gut(g1_gut) - gY), 1e-15)
        self.assertLess(
            abs(
                convert_g1(
                    g1=g1_gut,
                    from_conv=GaugeConvention.gut_hypercharge(),
                    to_conv=GaugeConvention.sm_hypercharge(),
                )
                - gY
            ),
            1e-15,
        )

    def test_rg_authority_enforced(self):
        """
        RG authority gate:
        - `tfpt_suite/data/rg_authority.json` must be present and parseable
        - legacy SM RG helpers (`rge_sm`) must not be used above MZ by default
        """
        pol = load_rg_authority_policy()
        self.assertEqual(pol.above_MZ, "PyR@TE_2loop_only")
        self.assertIn("SM", pol.models)

        mz = float(mz_scale_GeV())
        old = os.environ.pop("TFPT_ALLOW_LEGACY_SM_RGE_ABOVE_MZ", None)
        try:
            with self.assertRaises(RuntimeError):
                run_sm_gauge_only_2loop_thresholds(
                    mu_start_GeV=mz,
                    mu_end_GeV=2.0 * mz,
                    g_start=(0.5, 0.6, 1.2),
                    apply_alpha3_matching=False,
                )
        finally:
            if old is not None:
                os.environ["TFPT_ALLOW_LEGACY_SM_RGE_ABOVE_MZ"] = old

    def test_matching_is_finite_not_running(self):
        """
        Matching primitives must be pure/algebraic (no hidden RG running).
        """
        couplings = {"gY": 0.35, "g2": 0.65, "g3": 1.17}
        finite = {"alpha3": 1.0e-4}

        out1, meta1 = match_gauge(
            threshold_id="test",
            mu_thr_GeV=100.0,
            direction="up",
            couplings_below=couplings,
            loop_order=1,
            finite_delta_alpha=finite,
        )
        out2, meta2 = match_gauge(
            threshold_id="test",
            mu_thr_GeV=1000.0,  # should not affect the algebraic result
            direction="up",
            couplings_below=couplings,
            loop_order=1,
            finite_delta_alpha=finite,
        )
        self.assertEqual(meta1.status, "matched_with_finite_pieces")
        self.assertEqual(meta2.status, "matched_with_finite_pieces")
        self.assertAlmostEqual(out1["g3"], out2["g3"], places=15)

        # Check the sign convention: finite_delta_alpha encodes (alpha_above - alpha_below) for direction="up".
        alpha3_below = float((couplings["g3"] ** 2) / (4.0 * 3.141592653589793))
        alpha3_above = float(alpha3_below + finite["alpha3"])
        g3_expected = float(g_from_alpha(alpha3_above))
        self.assertAlmostEqual(float(out1["g3"]), g3_expected, places=12)

    def test_operator_spec_action_parse(self):
        action_path = TFPT_SUITE_DIR / "tfpt_suite" / "data" / "microscopic_action_tfpt_v25.json"
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "effective_action_r2_operator_spec.json"
            res = generate_effective_action_r2_operator_spec(microscopic_action_path=action_path, output_path=out_path)
        derivation = res.spec.get("derivation", {})
        action_parse = derivation.get("action_parse", {})
        self.assertTrue(action_parse.get("term_found", False))
        self.assertTrue(action_parse.get("tokens_ok", False))
        self.assertEqual(derivation.get("block_source"), "action_torsion_sector")

    def test_a0_b0_primitives(self):
        a0_val = a0_msbar(mass_GeV=TEST_MASS_GEV, mu_GeV=TEST_MU_GEV)
        expected_a0 = TEST_MASS_GEV**2 * (1.0 - math.log((TEST_MASS_GEV**2) / (TEST_MU_GEV**2)))
        self.assertAlmostEqual(a0_val, expected_a0, places=12)

        b0_val = b0_msbar(
            p_GeV=0.0,
            m1_GeV=TEST_MASS_GEV,
            m2_GeV=TEST_MASS_GEV,
            mu_GeV=TEST_MU_GEV,
            n_points=B0_TEST_POINTS,
        )
        expected_b0 = -math.log((TEST_MASS_GEV**2) / (TEST_MU_GEV**2))
        self.assertAlmostEqual(b0_val, expected_b0, places=9)

    def test_top_higgs_finite_pieces_finite(self):
        mt = 172.76
        mh = 125.25
        mw = 80.379
        mz = 91.1876

        delta_y_t = delta_y_t_ew_hempfling(m_t_GeV=mt, m_h_GeV=mh, mu_GeV=mt, v_ev_GeV=TEST_V_EV_GEV)
        delta_lambda = delta_lambda_h_1loop_buttazzo(
            m_h_GeV=mh,
            m_t_GeV=mt,
            m_w_GeV=mw,
            m_z_GeV=mz,
            mu_GeV=mh,
            v_ev_GeV=TEST_V_EV_GEV,
        )
        self.assertTrue(math.isfinite(delta_y_t))
        self.assertTrue(math.isfinite(delta_lambda))
        self.assertGreater(abs(delta_y_t), FINITE_PIECE_MIN)
        self.assertGreater(abs(delta_lambda), FINITE_PIECE_MIN)

    def test_core_invariants_runs(self):
        reg = get_module_registry()
        mod = reg["core_invariants"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_mass_spectrum_minimal_runs(self):
        reg = get_module_registry()
        mod = reg["mass_spectrum_minimal"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_mass_spectrum_deriver_runs(self):
        reg = get_module_registry()
        mod = reg["mass_spectrum_deriver"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_koide_constraints_runs(self):
        reg = get_module_registry()
        mod = reg["koide_constraints"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_bbn_neff_sanity_runs(self):
        reg = get_module_registry()
        mod = reg["bbn_neff_sanity"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_bbn_consistency_runs(self):
        reg = get_module_registry()
        mod = reg["bbn_consistency"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_baryogenesis_placeholder_runs(self):
        reg = get_module_registry()
        mod = reg["baryogenesis_placeholder"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_baryogenesis_mechanism_runs(self):
        reg = get_module_registry()
        mod = reg["baryogenesis_mechanism"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_gw_background_bounds_runs(self):
        reg = get_module_registry()
        mod = reg["gw_background_bounds"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_gw_background_predictor_runs(self):
        reg = get_module_registry()
        mod = reg["gw_background_predictor"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_g2_and_lamb_shift_proxy_runs(self):
        reg = get_module_registry()
        mod = reg["g2_and_lamb_shift_proxy"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_qed_anomalies_audit_runs(self):
        reg = get_module_registry()
        mod = reg["qed_anomalies_audit"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_arrow_mechanism_runs(self):
        reg = get_module_registry()
        mod = reg["arrow_mechanism"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertTrue(all(chk.passed for chk in res.checks))
        self.assertIn("arrow_mechanism_non_invertible", by_id)
        self.assertIn("arrow_prediction_falsifiable", by_id)

    def test_brst_ghost_deriver_runs(self):
        reg = get_module_registry()
        mod = reg["brst_ghost_deriver"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_boltzmann_transfer_runs(self):
        reg = get_module_registry()
        mod = reg["boltzmann_transfer"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_torsion_condensate_runs(self):
        reg = get_module_registry()
        mod = reg["torsion_condensate"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertTrue(all(chk.passed for chk in res.checks))
        self.assertIn("torsion_condensate_gap_equation_solved", by_id)
        self.assertIn("torsion_operator_spec_beta_R2_loaded", by_id)
        self.assertIn("n_quantization_source", by_id)

    def test_dark_energy_paths_runs(self):
        reg = get_module_registry()
        mod = reg["dark_energy_paths"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertTrue(all(chk.passed for chk in res.checks))
        self.assertIn("ladder_terminal_stage_identified", by_id)

    def test_torsion_observable_spin_fluid_runs(self):
        reg = get_module_registry()
        mod = reg["torsion_observable_spin_fluid"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertTrue(all(chk.passed for chk in res.checks))
        self.assertIn("experiment_specified_with_sensitivity", by_id)

    def test_torsion_observable_designer_runs(self):
        reg = get_module_registry()
        mod = reg["torsion_observable_designer"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_torsion_falsifiability_snr_runs(self):
        reg = get_module_registry()
        mod = reg["torsion_falsifiability_snr"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False, verification_mode="physics")
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("astro_channel_measurable_under_realistic_noise", by_id)
        self.assertEqual(by_id["astro_channel_measurable_under_realistic_noise"].severity, "PASS")
        self.assertIn("noise_psd_frequency_dependent", by_id)
        self.assertEqual(by_id["noise_psd_frequency_dependent"].severity, "PASS")
        self.assertIn("go_no_go_snr_ge_5", by_id)
        self.assertEqual(by_id["go_no_go_snr_ge_5"].severity, "PASS")

    def test_dm_alternative_channels_runs(self):
        reg = get_module_registry()
        mod = reg["dm_alternative_channels"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_flavor_topology_mapper_runs(self):
        reg = get_module_registry()
        mod = reg["flavor_topology_mapper"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_defect_partition_derivation_runs(self):
        reg = get_module_registry()
        mod = reg["defect_partition_derivation"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))
        justification = res.results.get("delta2_justification", {})
        self.assertEqual(justification.get("effective_multiplicity"), EXPECTED_DEFECT_MULTIPLICITY)
        class_ids = {entry.get("class_id") for entry in justification.get("equivalence_classes", []) if isinstance(entry, dict)}
        self.assertTrue(EXPECTED_DEFECT_CLASSES.issubset(class_ids))

    def test_alpha_on_shell_bridge_runs(self):
        reg = get_module_registry()
        mod = reg["alpha_on_shell_bridge"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))
        by_id = {c.check_id: c for c in res.checks}
        for check_id in EXPECTED_ALPHA_BRIDGE_CHECKS:
            self.assertIn(check_id, by_id)
            self.assertTrue(by_id[check_id].passed)

    def test_axion_scenario_matrix_runs(self):
        reg = get_module_registry()
        mod = reg["axion_scenario_matrix"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_axion_dm_pipeline_physics_mode_no_fail(self):
        reg = get_module_registry()
        mod = reg["axion_dm_pipeline"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False, verification_mode="physics")
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_topology_phase_map_runs(self):
        reg = get_module_registry()
        mod = reg["topology_phase_map"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_flavor_joint_objective_scan_runs(self):
        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            # prerequisites (write their results.json to the same output_dir)
            reg["topology_phase_map"].run_and_write(config=cfg)
            reg["ckm_full_pipeline"].run_and_write(config=cfg)
            reg["pmns_full_pipeline"].run_and_write(config=cfg)
            res = reg["flavor_joint_objective_scan"].run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("joint_objective_computed", by_id)
        self.assertIn("ckm_variants_present", by_id)
        self.assertTrue(by_id["ckm_variants_present"].passed)
        # Ensure the module actually produced a non-empty scan table (prevents silent NaN/empty wiring).
        self.assertTrue(bool(res.results.get("ckm", {}).get("variants", [])))
        self.assertIsNotNone(res.results.get("best", None))

    def test_uncertainty_propagator_runs(self):
        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            # prerequisites for reading MC summaries
            reg["msbar_matching_map"].run_and_write(config=cfg)
            reg["ckm_full_pipeline"].run_and_write(config=cfg)
            reg["pmns_full_pipeline"].run_and_write(config=cfg)
            res = reg["uncertainty_propagator"].run_and_write(config=cfg)
        self.assertIn("uncertainty_propagator_runs", {c.check_id for c in res.checks})

    def test_matching_finite_pieces_runs(self):
        reg = get_module_registry()
        mod = reg["matching_finite_pieces"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))
        by_id = {chk.check_id: chk for chk in res.checks}
        self.assertIn("alpha3_matching_invertible", by_id)
        # Regression: the 2-loop matching primitive must be numerically invertible
        # (avoid O(α^5) drift when composing down∘up).
        self.assertEqual(by_id["alpha3_matching_invertible"].severity, "PASS")

    def test_arrow_of_time_proxy_runs(self):
        reg = get_module_registry()
        mod = reg["arrow_of_time_proxy"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_torsion_dm_pipeline_runs(self):
        reg = get_module_registry()
        mod = reg["torsion_dm_pipeline"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_effective_action_r2_runs(self):
        reg = get_module_registry()
        mod = reg["effective_action_r2"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {chk.check_id: chk.passed for chk in res.checks}
        self.assertIn("torsion_operator_specified_for_a2_derivation", by_id)
        self.assertTrue(by_id["torsion_operator_specified_for_a2_derivation"])
        self.assertTrue(by_id.get("R2_scale_formula", False))

    def test_aps_eta_gluing_runs(self):
        reg = get_module_registry()
        mod = reg["aps_eta_gluing"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_discrete_consistency_uniqueness_runs(self):
        reg = get_module_registry()
        mod = reg["discrete_consistency_uniqueness"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=70, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_discrete_consistency_uniqueness_plot_runs_without_warnings(self):
        reg = get_module_registry()
        mod = reg["discrete_consistency_uniqueness"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=70, seed=0, overwrite=True, plot=True)
            res = mod.run_and_write(config=cfg)
            self.assertEqual(res.warnings, [])
            plot_path = Path(str(res.results.get("plot", {}).get("cfe_uniqueness_png", "")))
            self.assertTrue(plot_path.is_file(), f"Expected plot file, got: {plot_path}")

    def test_ckm_full_pipeline_runs(self):
        reg = get_module_registry()
        mod = reg["ckm_full_pipeline"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_pmns_mechanism_bridge_runs(self):
        reg = get_module_registry()
        mod = reg["pmns_mechanism_bridge"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_two_loop_rg_fingerprints_runs(self):
        reg = get_module_registry()
        mod = reg["two_loop_rg_fingerprints"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_unification_gate_runs(self):
        reg = get_module_registry()
        mod = reg["unification_gate"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertIn("unification_gate", {c.check_id for c in res.checks})

    def test_msbar_matching_map_runs(self):
        reg = get_module_registry()
        mod = reg["msbar_matching_map"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_axion_fa_derivation_runs(self):
        reg = get_module_registry()
        mod = reg["axion_fa_derivation"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("f_a_derived_not_quoted", by_id)
        self.assertTrue(by_id["f_a_derived_not_quoted"].passed)

    def test_axion_dm_pipeline_prefers_derived_fa(self):
        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            cfg = SuiteConfig(output_dir=out, mp_dps=80, seed=0, overwrite=True, plot=False)
            reg["axion_fa_derivation"].run_and_write(config=cfg)
            res = reg["axion_dm_pipeline"].run_and_write(config=cfg)
        axion_claim = res.results.get("axion_claim", {})
        self.assertEqual(axion_claim.get("source"), "derived")
        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("c_str_explained_by_topology", by_id)
        self.assertTrue(by_id["c_str_explained_by_topology"].passed)

    def test_global_consistency_test_physics_mode_passes_alpha_gate(self):
        reg = get_module_registry()
        mod = reg["global_consistency_test"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False, verification_mode="physics")
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_likelihood_engine_runs(self):
        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=80, seed=0, overwrite=True, plot=False, verification_mode="physics")
            # prerequisite: global consistency predictions provider
            reg["global_consistency_test"].run_and_write(config=cfg)
            dm_dir = Path(tmp) / "axion_dm_pipeline"
            dm_dir.mkdir(parents=True, exist_ok=True)
            dm_payload = {"results": {"relic_density": {"Omega_a_h2": TEST_OMEGA_A_H2}}, "checks": []}
            (dm_dir / "results.json").write_text(json.dumps(dm_payload), encoding="utf-8")
            res = reg["likelihood_engine"].run_and_write(config=cfg)
        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("likelihood_engine_runs", by_id)
        self.assertTrue(by_id["likelihood_engine_runs"].passed)
        self.assertTrue(by_id.get("covariance_enabled_for_reference_tables").passed)
        self.assertIn("nufit_pmns_grid_plugin_active", by_id)
        self.assertIn("unified_score_p_value_reported", by_id)
        datasets = res.results.get("datasets", [])
        alpha_dataset = next((row for row in datasets if row.get("dataset_id") == "alpha_covariance_minimal"), None)
        self.assertIsNotNone(alpha_dataset)
        self.assertTrue(math.isfinite(float(alpha_dataset.get("chi2", float("nan")))))
        unified = res.results.get("unified_score", {})
        components = unified.get("components", [])
        alpha_comp = next((row for row in components if row.get("sector") == "alpha"), None)
        dm_comp = next((row for row in components if row.get("sector") == "dm"), None)
        self.assertIsNotNone(alpha_comp)
        self.assertIsNotNone(dm_comp)

    def test_bounce_injection_wired_into_boltzmann_transfer(self):
        # Check if CAMB is available (required for full bounce injection testing).
        try:
            import camb  # noqa: F401
            camb_available = True
        except ImportError:
            camb_available = False

        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            cfg = SuiteConfig(output_dir=out, mp_dps=60, seed=0, overwrite=True, plot=False, verification_mode="physics")

            # Minimal synthetic bounce output (keeps this a fast unit test).
            bounce_dir = out / "bounce_perturbations"
            bounce_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "results": {
                    "k_grid": [0.1, 0.3, 1.0, 3.0, 10.0],
                    "T_scalar": [1.0, 1.0, 1.0, 1.0, 1.0],
                    "T_tensor": [1.0, 1.0, 1.0, 1.0, 1.0],
                    "diagnostics": {"k_bounce_s_est_raw": 10.0, "k_bounce_t_est_raw": 10.0},
                },
                "checks": [],
            }
            (bounce_dir / "results.json").write_text(json.dumps(payload), encoding="utf-8")

            reg["k_calibration"].run_and_write(config=cfg)
            reg["primordial_spectrum_builder"].run_and_write(config=cfg)
            res = reg["boltzmann_transfer"].run_and_write(config=cfg)

        by_id = {c.check_id: c for c in res.checks}
        self.assertIn("bounce_feature_injection_wired", by_id)
        # When CAMB is available, bounce injection should pass; otherwise, it returns INFO.
        if camb_available:
            self.assertEqual(by_id["bounce_feature_injection_wired"].severity, "PASS")
            self.assertIn("planck_lowl_evaluated", by_id)
            self.assertIn("planck_lensing_evaluated", by_id)
            signature = res.results.get("signature_policy", {})
            self.assertEqual(signature.get("primary_signature"), "tensor_CMB")
        else:
            self.assertEqual(by_id["bounce_feature_injection_wired"].severity, "INFO")
        self.assertIn("signature_policy_declared", by_id)
        self.assertIn("signature_policy_consistent_with_bounce", by_id)

    def test_k_calibration_prefers_threshold_history(self):
        reg = get_module_registry()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            thr_dir = out / "cosmo_threshold_history"
            thr_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "results": {
                    "reheating": {"N_reheat": 1.234, "T_reheat_GeV": 5.678, "g_star_s_reheat": 98.0},
                    "pivot": {"N_pivot": 52.0},
                }
            }
            (thr_dir / "results.json").write_text(json.dumps(payload), encoding="utf-8")

            cfg = SuiteConfig(output_dir=out, mp_dps=60, seed=0, overwrite=True, plot=False)
            res = reg["k_calibration"].run_and_write(config=cfg)

        assumptions = res.results.get("assumptions", {})
        self.assertAlmostEqual(float(assumptions.get("N_reheat")), 1.234, places=6)
        self.assertAlmostEqual(float(assumptions.get("T_reheat_GeV")), 5.678, places=6)
        self.assertAlmostEqual(float(assumptions.get("g_star_s_reheat")), 98.0, places=6)
        self.assertAlmostEqual(float(assumptions.get("N_inflation_from_transition")), 52.0, places=6)

    def test_stability_unitarity_audit_runs(self):
        reg = get_module_registry()
        mod = reg["stability_unitarity_audit"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))

    def test_below_mt_eft_cascade_runs(self):
        reg = get_module_registry()
        mod = reg["below_mt_eft_cascade"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks))


if __name__ == "__main__":
    unittest.main()

