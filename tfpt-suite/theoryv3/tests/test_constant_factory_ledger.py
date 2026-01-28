import tempfile
import unittest
from pathlib import Path

from tfpt_suite.config import SuiteConfig
from theoryv3_suite.modules.constant_factory_audit import ConstantFactoryAuditModule


class ConstantFactoryLedgerTests(unittest.TestCase):
    def test_ledger_structure_and_sensitivity(self) -> None:
        module = ConstantFactoryAuditModule()
        with tempfile.TemporaryDirectory() as tmpdir:
            config = SuiteConfig(output_dir=Path(tmpdir), plot=False)
            result = module.run_and_write(config=config)
            ledger = result.results.get("ledger", [])
            self.assertTrue(isinstance(ledger, list))
            by_name = {entry.get("name"): entry for entry in ledger if isinstance(entry, dict)}
            self.assertIn("g_a_gamma_gamma", by_name)
            self.assertIn("beta_rad", by_name)
            self.assertIn("M_over_Mpl", by_name)

            g_entry = by_name["g_a_gamma_gamma"]
            beta_entry = by_name["beta_rad"]
            m_entry = by_name["M_over_Mpl"]

            g_sens = g_entry.get("sensitivity", {})
            beta_sens = beta_entry.get("sensitivity", {})
            m_sens = m_entry.get("sensitivity", {})

            self.assertAlmostEqual(float(g_sens.get("Sc3", 0.0)), 1.0, places=2)
            self.assertAlmostEqual(float(beta_sens.get("Sphi0", 0.0)), 1.0, places=2)
            self.assertAlmostEqual(float(m_sens.get("Sc3", 0.0)), 4.0, places=1)


if __name__ == "__main__":
    unittest.main()
