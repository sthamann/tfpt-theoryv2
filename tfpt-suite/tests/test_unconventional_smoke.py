import os
import sys
import tempfile
import unittest
from pathlib import Path


# Ensure `tfpt_suite` and `tfpt_unconventional` packages are importable when tests are run from workspace root.
TFPT_SUITE_DIR = Path(__file__).resolve().parents[1]
UNCONVENTIONAL_DIR = TFPT_SUITE_DIR / "unconventional"
for p in (TFPT_SUITE_DIR, UNCONVENTIONAL_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


from tfpt_suite.config import SuiteConfig  # noqa: E402
from tfpt_unconventional.modules.registry import get_unconventional_module_registry  # noqa: E402


class TestUnconventionalSmoke(unittest.TestCase):
    def test_registry_has_expected_modules(self):
        reg = get_unconventional_module_registry()
        for module_id in [
            "ux_matching_metamorphic_audit",
            "ux_cosmo_history_sampler",
            "ux_omega_b_aps_bridge",
            "ux_threshold_graph_audit",
            "ux_flavor_holdout_search",
            "ux_gravity_gaugefix_ga",
            "ux_torsion_regime_designer",
        ]:
            self.assertIn(module_id, reg)

    def test_matching_metamorphic_audit_runs(self):
        reg = get_unconventional_module_registry()
        mod = reg["ux_matching_metamorphic_audit"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks), f"Failed checks: {[c for c in res.checks if not c.passed]}")

    def test_omega_b_aps_bridge_runs(self):
        reg = get_unconventional_module_registry()
        mod = reg["ux_omega_b_aps_bridge"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks), f"Failed checks: {[c for c in res.checks if not c.passed]}")

    def test_gravity_gaugefix_ga_runs(self):
        reg = get_unconventional_module_registry()
        mod = reg["ux_gravity_gaugefix_ga"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(all(chk.passed for chk in res.checks), f"Failed checks: {[c for c in res.checks if not c.passed]}")

    def test_torsion_regime_designer_runs(self):
        reg = get_unconventional_module_registry()
        mod = reg["ux_torsion_regime_designer"]
        with tempfile.TemporaryDirectory() as tmp:
            cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=60, seed=0, overwrite=True, plot=False)
            res = mod.run_and_write(config=cfg)
        self.assertTrue(any("proposals" in str(k) for k in res.results.keys()))

    def test_cosmo_history_sampler_runs_fast(self):
        reg = get_unconventional_module_registry()
        mod = reg["ux_cosmo_history_sampler"]

        old_plausible = os.environ.get("TFPT_UX_COSMO_SAMPLES_PLAUSIBLE")
        old_extended = os.environ.get("TFPT_UX_COSMO_SAMPLES_EXTENDED")
        os.environ["TFPT_UX_COSMO_SAMPLES_PLAUSIBLE"] = "30"
        os.environ["TFPT_UX_COSMO_SAMPLES_EXTENDED"] = "30"
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = SuiteConfig(output_dir=Path(tmp), mp_dps=50, seed=0, overwrite=True, plot=False)
                res = mod.run_and_write(config=cfg)
        finally:
            if old_plausible is None:
                os.environ.pop("TFPT_UX_COSMO_SAMPLES_PLAUSIBLE", None)
            else:
                os.environ["TFPT_UX_COSMO_SAMPLES_PLAUSIBLE"] = old_plausible
            if old_extended is None:
                os.environ.pop("TFPT_UX_COSMO_SAMPLES_EXTENDED", None)
            else:
                os.environ["TFPT_UX_COSMO_SAMPLES_EXTENDED"] = old_extended

        # Only basic sanity: module ran and produced expected top-level keys.
        self.assertIn("scenarios", res.results)
        self.assertTrue(any(chk.check_id == "chi_star_computed" for chk in res.checks))


if __name__ == "__main__":
    unittest.main()

