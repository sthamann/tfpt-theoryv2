import tempfile
import unittest
from pathlib import Path

from tfpt_suite.config import SuiteConfig
from theoryv3_suite.modules.registry import get_module_registry
from theoryv3_suite.modules.seed_invariants_audit import SeedInvariantsAuditModule
from theoryv3_suite.modules.yukawa_exponent_index_audit import YukawaExponentIndexAuditModule


class TheoryV3SmokeTests(unittest.TestCase):
    def test_registry_contains_modules(self) -> None:
        registry = get_module_registry()
        expected = {
            "seed_invariants_audit",
            "defect_partition_g5_audit",
            "alpha_backreaction_sensitivity_audit",
            "g5_origin_audit",
            "dark_energy_exponential_audit",
            "dark_energy_norm_half_origin_audit",
            "flavor_pattern_audit",
            "pmns_tm1_audit",
            "yukawa_exponent_index_audit",
            "yukawa_index_mapping_audit",
            "baryon_consistency_audit",
            "axion_dm_audit",
            "g5_crosslink_audit",
            "constant_factory_audit",
        }
        self.assertTrue(expected.issubset(set(registry.keys())))

    def test_seed_invariants_runs(self) -> None:
        module = SeedInvariantsAuditModule()
        with tempfile.TemporaryDirectory() as tmpdir:
            config = SuiteConfig(output_dir=Path(tmpdir), plot=False)
            result = module.run_and_write(config=config)
            self.assertIn("values", result.results)
            self.assertTrue(Path(tmpdir, module.module_id, "results.json").is_file())

    def test_yukawa_index_runs(self) -> None:
        module = YukawaExponentIndexAuditModule()
        with tempfile.TemporaryDirectory() as tmpdir:
            config = SuiteConfig(output_dir=Path(tmpdir), plot=False)
            result = module.run_and_write(config=config)
            rows = result.results.get("rows", [])
            self.assertTrue(isinstance(rows, list))
            self.assertGreaterEqual(len(rows), 1)


if __name__ == "__main__":
    unittest.main()
