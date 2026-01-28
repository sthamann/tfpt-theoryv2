# tests

Minimal `unittest` tests for `tfpt-suite`.

These are intentionally **small and fast**:
- they run with `plot=False` (headless / CI-friendly)
- they write artifacts into a temporary directory (no repo pollution)

## What is covered

- `test_smoke.py`
  - verifies module registry wiring (`get_module_registry()`)
  - runs a representative subset of modules end-to-end (including some integration-style checks for flavor/matching scaffolds)
  - asserts the defect-partition derivation exposes a g=5 equivalence-class justification
  - asserts the α on-shell bridge exposes explicit leptonic/hadronic/EW decoupling checks
  - validates the A0/B0 scalar-integral primitives and finite EW pieces (Δy_t, Δλ_H) are finite and non-trivial
  - validates OperatorSpec action-term parsing metadata (`block_source=action_torsion_sector`)
  - ensures the threshold-driven reheating module is registered (`cosmo_threshold_history`)
  - verifies `k_calibration` prefers threshold-derived reheating inputs when present
  - checks that `boltzmann_transfer` reports signature-policy and Planck low-ℓ/lensing hooks
  - validates `axion_fa_derivation`, derived f_a usage, `c_str_explained_by_topology`, torsion gap-equation checks, ladder terminal-stage identification, frequency-dependent noise PSD wiring, experiment-specification checks, NuFIT grid plugin wiring, and unified scorecard reporting (including alpha/DM contributions and relic_density fallback)
  - runtime may be ~O(1 minute) on a typical laptop due to CKM/PMNS pipeline integration tests

- `test_flavor_texture.py`
  - unit tests for the low-level linear algebra helpers used by CKM/PMNS pipelines
  - protects basic invariants like unitarity / hermiticity and deterministic conventions

Run:

```bash
python3 -m unittest discover -s tfpt-suite/tests
```

## Notes

- To verify plot generation, run modules via `run_suite.py` without `--no-plot` (plots are optional by design, and should never make a module fail).

