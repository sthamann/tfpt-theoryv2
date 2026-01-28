# TFPT Suite – Complete Test Documentation

sThis document provides **detailed documentation for every test module** across the three test suites:

1. **tfpt-suite** (63 registry-driven modules)
2. **theoryv3** (14 analysis modules)
3. **unconventional** (7 search/audit modules)
4. **Unit Tests** (6 test files)

---

## Table of Contents

1. [Quick Reference: All Modules](#quick-reference-all-modules)
2. [How to Read Suite Outputs](#how-to-read-suite-outputs)
3. [What "Tests" Mean in TFPT Suite](#what-tests-mean-in-tfpt-suite)
4. [Input Tables and Parameters](#input-tables-and-parameters)
5. [Unit Tests (Python unittest)](#unit-tests-python-unittest)
6. [TFPT Suite – Conventional Modules (63)](#tfpt-suite--conventional-modules-63)
7. [theoryv3 Suite – Analysis Modules (14)](#theoryv3-suite--analysis-modules-14)
8. [Unconventional Suite – Search/Audit Modules (7)](#unconventional-suite--searchaudit-modules-7)
9. [Summary: Explicit Missing Bridges](#summary-explicit-missing-bridges)

---

## Quick Reference: All Modules

### TFPT Suite (63 modules)

| # | Module ID | Category |
|---|-----------|----------|
| 1 | `core_invariants` | Foundation |
| 2 | `ufe_gravity_normalization` | Foundation |
| 3 | `brst_ghost_deriver` | Foundation |
| 4 | `chiral_index_three_cycles` | Foundation |
| 5 | `mass_spectrum_minimal` | Mass Spectrum |
| 6 | `mass_spectrum_deriver` | Mass Spectrum |
| 7 | `bbn_neff_sanity` | Cosmology |
| 8 | `bbn_consistency` | Cosmology |
| 9 | `koide_constraints` | Flavor |
| 10 | `two_loop_rg_fingerprints` | RG/Gauge |
| 11 | `unification_gate` | RG/Gauge |
| 12 | `alpha_precision_audit` | α Sector |
| 13 | `predictions_dashboard` | Summary |
| 14 | `qft_completeness_ledger` | Audit |
| 15 | `anomaly_cancellation_audit` | Audit |
| 16 | `effective_action_r2` | Gravity |
| 17 | `bounce_perturbations` | Cosmology |
| 18 | `gw_background_bounds` | Cosmology |
| 19 | `gw_background_predictor` | Cosmology |
| 20 | `cosmo_threshold_history` | Cosmology |
| 21 | `cosmo_reheating_policy_v106` | Cosmology |
| 22 | `k_calibration` | Cosmology |
| 23 | `primordial_spectrum_builder` | Cosmology |
| 24 | `boltzmann_transfer` | Cosmology |
| 25 | `matching_finite_pieces` | Matching |
| 26 | `msbar_matching_map` | Matching |
| 27 | `below_mt_eft_cascade` | EFT |
| 28 | `stability_unitarity_audit` | Audit |
| 29 | `aps_eta_gluing` | Topology |
| 30 | `discrete_consistency_uniqueness` | α Sector |
| 31 | `defect_partition_derivation` | α Sector |
| 32 | `alpha_on_shell_bridge` | α Sector |
| 33 | `discrete_complexity_minimizer` | Search |
| 34 | `mobius_cusp_classification` | Topology |
| 35 | `mobius_delta_calibration` | Topology |
| 36 | `topology_phase_map` | Topology |
| 37 | `seesaw_block` | Neutrino |
| 38 | `birefringence_tomography` | CMB |
| 39 | `omega_b_conjecture_scan` | Baryon |
| 40 | `axion_fa_derivation` | Axion/DM |
| 41 | `axion_dm_pipeline` | Axion/DM |
| 42 | `axion_scenario_matrix` | Axion/DM |
| 43 | `dark_energy_paths` | Dark Energy |
| 44 | `torsion_condensate` | Torsion |
| 45 | `torsion_bounds_mapping` | Torsion |
| 46 | `torsion_observable_spin_fluid` | Torsion |
| 47 | `torsion_observable_designer` | Torsion |
| 48 | `torsion_falsifiability_snr` | Torsion |
| 49 | `torsion_dm_pipeline` | Torsion |
| 50 | `dm_alternative_channels` | DM |
| 51 | `ckm_full_pipeline` | Flavor |
| 52 | `mobius_z3_yukawa_generator` | Flavor |
| 53 | `pmns_full_pipeline` | Flavor |
| 54 | `flavor_topology_mapper` | Flavor |
| 55 | `flavor_joint_objective_scan` | Flavor |
| 56 | `uncertainty_propagator` | Statistics |
| 57 | `pmns_mechanism_bridge` | Neutrino |
| 58 | `pmns_z3_breaking` | Neutrino |
| 59 | `global_consistency_test` | Summary |
| 60 | `likelihood_engine` | Statistics |
| 61 | `baryogenesis_placeholder` | Baryon |
| 62 | `baryogenesis_mechanism` | Baryon |
| 63 | `g2_and_lamb_shift_proxy` | Precision |
| 64 | `qed_anomalies_audit` | Precision |
| 65 | `arrow_of_time_proxy` | Cosmology |
| 66 | `arrow_mechanism` | Cosmology |

### theoryv3 Suite (14 modules)

| # | Module ID | Purpose |
|---|-----------|---------|
| 1 | `seed_invariants_audit` | π → core invariants consistency |
| 2 | `defect_partition_g5_audit` | δ₂ from g=5 and α(0) closure |
| 3 | `alpha_backreaction_sensitivity_audit` | k sweep (α sensitivity vs backreaction) |
| 4 | `g5_origin_audit` | g=5 from SU(5) holonomy degeneracy |
| 5 | `dark_energy_exponential_audit` | exp(−α⁻¹/2) suppression and ρ_Λ |
| 6 | `dark_energy_norm_half_origin_audit` | n=1/2 from double-cover degree |
| 7 | `flavor_pattern_audit` | λ, δ*, δ_CP, PMNS θ₁₃ checks |
| 8 | `pmns_tm1_audit` | TM1 sum rule for sin²θ₁₂ |
| 9 | `yukawa_exponent_index_audit` | Rational indices for mass ratios |
| 10 | `yukawa_index_mapping_audit` | Map q_ij to charge-squared index sums |
| 11 | `baryon_consistency_audit` | Ω_b identity, η_b proxy, derived H₀ |
| 12 | `axion_dm_audit` | Axion frequency and relic fraction audit |
| 13 | `g5_crosslink_audit` | g=5 consistency across sectors |
| 14 | `constant_factory_audit` | Grouped constant ledger with derivations |

### Unconventional Suite (7 modules)

| # | Module ID | Purpose |
|---|-----------|---------|
| 1 | `ux_matching_metamorphic_audit` | Metamorphic invariants of matching primitives |
| 2 | `ux_cosmo_history_sampler` | k→ℓ bridge feasibility sampling |
| 3 | `ux_flavor_holdout_search` | CKM convention holdout validation |
| 4 | `ux_gravity_gaugefix_ga` | GA search for gauge-fixing choices |
| 5 | `ux_omega_b_aps_bridge` | APS seam → Ω_b coefficient docking |
| 6 | `ux_threshold_graph_audit` | Threshold graph visualization |
| 7 | `ux_torsion_regime_designer` | Falsifiable torsion regime design tool |

### Unit Tests (6 files)

| Suite | File | Purpose |
|-------|------|---------|
| tfpt-suite | `test_smoke.py` | Registry wiring + integration checks |
| tfpt-suite | `test_flavor_texture.py` | Flavor linear algebra helpers |
| tfpt-suite | `test_unconventional_smoke.py` | Unconventional registry smoke test |
| theoryv3 | `test_theoryv3_smoke.py` | theoryv3 registry smoke test |
| theoryv3 | `test_constant_factory_ledger.py` | Constant ledger sensitivity validation |

---

## How to Read Suite Outputs

### Per-Module Artifacts

Each module writes to its output directory:

| Suite | Output Location |
|-------|-----------------|
| tfpt-suite | `tfpt-suite/out/<module_id>/` or `out_physics/` |
| theoryv3 | `tfpt-suite/theoryv3/out/<module_id>/` |
| unconventional | `tfpt-suite/out/unconventional/<module_id>/` |

| File | Content | Purpose |
|------|---------|---------|
| `results.json` | All numerical results + `checks[]` + `warnings[]` | Machine parsing, downstream modules |
| `report.txt` | Human-readable summary | Quick inspection |
| `meta.json` | Runtime, seed, config, SHA256 fingerprints | Reproducibility audit |
| `*.png` | Diagnostic plots | Visual verification |

### Aggregated Reports

| Report | Location | Mode |
|--------|----------|------|
| Engineering PDF | `tfpt-suite/tfpt-test-results.pdf` | Default, focus on determinism |
| Physics PDF | `tfpt-suite/tfpt-test-results-physics.pdf` | Strict, flags deviations |
| theoryv3 PDF | `tfpt-suite/theoryv3/theoryv3-analysis.pdf` | Discrete building blocks |
| Unconventional PDF | `tfpt-suite/unconventional/unconventional-test-results.pdf` | Search/audit tools |

---

## What "Tests" Mean in TFPT Suite

Three verification layers:

### 1. Unit Tests (`unittest`)
Fast Python tests for registry wiring and low-level invariants.
- Location: `tfpt-suite/tests/`, `tfpt-suite/theoryv3/tests/`
- Run: `python -m pytest tfpt-suite/tests/`

### 2. Module Self-Checks
Each module returns `checks[]` with `passed=true/false` and details.
- Newer modules attach `severity` (PASS|WARN|FAIL|INFO)
- Physics mode upgrades large deviations to WARN/FAIL
- Location: `<output_dir>/<module_id>/results.json`

### 3. Reproducibility Checks
`meta.json` plus SHA256 fingerprints for critical inputs.
- Example: `two_loop_rg_fingerprints` stores PyR@TE output hashes

---

## Input Tables and Parameters

### TFPT Internal (No External Data)

| Input | Source | Used By |
|-------|--------|---------|
| c₃ = 1/(8π) | Computed | Core invariants, all modules |
| φ₀ = 1/(6π) + δ_top | Computed | CFE, α sector, flavor |
| β = φ₀/(4π) | Computed | Birefringence, Ω_b |

### External Reference Tables

| File | Content | Used By |
|------|---------|---------|
| `global_reference.json` | Observable means/σ | Dashboard, global χ² |
| `sm_inputs_mz.json` | SM parameters at M_Z | Matching, RG runs |
| `ckm_reference.json` | CKM snapshot | CKM χ² diagnostics |
| `pmns_reference.json` | PMNS snapshot | PMNS χ² diagnostics |
| `lepton_masses_pdg.json` | Lepton pole masses | δ calibration |
| `rge_thresholds_v25.json` | Threshold scales | RG segmentation |
| `birefringence_tomography.json` | CMB birefringence data | Tomography fit |
| `torsion_bounds_vetted.json` | SME bounds | Torsion falsifiability |
| `torsion_regimes.json` | TFPT regime choices | Torsion predictions |
| `torsion_falsifiability_noise_v1.json` | Source+noise policy | Torsion SNR gate |
| `two_loop_rg_fingerprints.json` | Model selector + policy | RG fingerprinting |
| `effective_action_r2_operator_spec.json` | R² operator contract | Heat-kernel evaluation |
| `flavor_texture_v24.json` | Yukawa conventions | CKM/PMNS generation |
| `k_calibration.json` | Cosmology snapshot | k→ℓ mapping |
| `axion_tfpt_v106.json` | Axion parameters | DM pipeline |

---

## Unit Tests (Python unittest)

### tfpt-suite/tests/

#### `test_smoke.py`

**Purpose**: Catch import/schema regressions and registry wiring early.

**What it tests**:
- Module registry can be built (`get_module_registry()`)
- Subset of key modules runs end-to-end with `plot=False`
- Defect-partition derivation exposes g=5 equivalence-class justification
- α on-shell bridge exposes decoupling checks
- A₀/B₀ scalar-integral primitives and finite EW pieces
- OperatorSpec action-term parsing metadata
- Threshold-driven reheating module registration
- `k_calibration` prefers threshold-derived inputs
- `boltzmann_transfer` signature-policy and Planck hooks
- Axion f_a derivation and torsion gap-equation checks

**How to run**:
```bash
python -m pytest tfpt-suite/tests/test_smoke.py -v
```

#### `test_flavor_texture.py`

**Purpose**: Prevent "physically wrong labeling" artifacts in linear algebra helpers.

**What it tests**:
- Unitarity/hermiticity of generated matrices
- Stable conventions across regenerations
- Deterministic labeling

#### `test_unconventional_smoke.py`

**Purpose**: Verify unconventional suite registry wiring.

**What it tests**:
- Unconventional module registry builds correctly
- Representative modules run end-to-end

### theoryv3/tests/

#### `test_theoryv3_smoke.py`

**Purpose**: Verify theoryv3 suite registry and basic execution.

**What it tests**:
- theoryv3 module registry builds correctly
- Representative modules run end-to-end with `plot=False`

#### `test_constant_factory_ledger.py`

**Purpose**: Validate constant factory ledger structure.

**What it tests**:
- Sensitivity fields for key invariants
- Ledger schema compliance

**How to run**:
```bash
python tfpt-suite/theoryv3/run_theoryv3.py run-tests
```

---

## TFPT Suite – Conventional Modules (63)

### Foundation & Core

---

#### 1. `core_invariants` — Fundamental TFPT Constants

**What it tests**: Computation and validation of TFPT's fundamental algebraic invariants.

**Inputs**: None (purely algebraic).

**Calculation**:
- Evaluates closed-form identities with 80 dps precision (mpmath)
- Computes c₃, φ₀, δ_top, β, g_aγγ, M/M_Pl

**Key Results**:

| Quantity | Formula | Value |
|----------|---------|-------|
| c₃ | 1/(8π) | 0.03978873577... |
| φ₀ (tree) | 1/(6π) | 0.05305164769... |
| δ_top | 3/(256π⁴) | 0.00012030447... |
| φ₀ (total) | φ₀_tree + δ_top | 0.05317195217... |
| g_aγγ | -1/(2π) | -0.15915494309... |
| β_rad | φ₀/(4π) | 0.00423128951... |
| β_deg | — | **0.2424°** |
| M/M_Pl | √(8π)·c₃⁴ | **1.256×10⁻⁵** |

**Status**: All PASS. No gaps.

---

#### 2. `ufe_gravity_normalization` — ξ=c₃/φ₀ Normalization

**What it tests**: Gravity-scalar normalization parameter derivation.

**Key Results**:
- ξ_tree = 3/4
- δ_top shift recorded
- OperatorSpec/BRST closure audited

**Status**: PASS

---

#### 3. `brst_ghost_deriver` — BRST/Ghost Closure

**What it tests**: BRST/ghost sector closure from OperatorSpec.

**Key Results**:
- OperatorSpec derivation validated
- Gauge-scan for closure verification

**Status**: PASS

---

#### 4. `chiral_index_three_cycles` — 3-Family Index

**What it tests**: Chiral index theorem for three generations.

**Key Results**: Ind(D)=3 for (1,1,1) fluxes

**Status**: PASS

---

### α Sector

---

#### 5. `discrete_consistency_uniqueness` — CFE Root Uniqueness

**What it tests**: The Closure Fixpoint Equation (CFE) has exactly one positive α root.

**Inputs**: TFPT invariants (c₃, φ₀, b₁).

**Parameters**:
- Baseline: Fixed φ₀
- Self-consistent: φ(α) = φ_tree + δ_top·e^(-kα), with k=2 (double-cover)

**Key Results**:

| Mode | α⁻¹ | Status |
|------|-----|--------|
| Baseline | **137.03650146** | Matches paper |
| Self-consistent (k=2) | **137.03599410** | Matches paper |

- Uniqueness: f(α*) = -3.82×10⁻⁷ < 0 at local minimum → exactly one positive root

**Status**: All 5 PASS.

---

#### 6. `alpha_precision_audit` — α vs Experiment

**What it tests**: TFPT's α prediction against experimental references.

**Inputs**:
- Primary: PDG ᾱ⁵⁻¹(M_Z) = 127.951 ± 0.009
- Secondary: CODATA α⁻¹(0) = 137.035999177 ± 2.1×10⁻⁸

**Key Results**:

| Model | α⁻¹(0) | ppm vs CODATA |
|-------|--------|---------------|
| Baseline (k=0) | 137.03650146 | +3.67 |
| **Canonical (k=2)** | **137.03599410** | **-0.037** |
| Two-defect | 137.03599819 | -0.0072 |

**Status**: All 5 PASS.

---

#### 7. `defect_partition_derivation` — δ₂ Derivation

**What it tests**: Two-defect partition derivation.

**Key Results**: α(0) within 2σ diagnostic

**Status**: PASS

---

#### 8. `alpha_on_shell_bridge` — On-Shell α Connection

**What it tests**: On-shell α connection with decoupling checks.

**Status**: PASS

---

### RG & Gauge Running

---

#### 9. `two_loop_rg_fingerprints` — Gauge Running Fingerprints

**What it tests**: Where α₃(μ) crosses TFPT invariants under 2-loop RG evolution.

**Inputs**:
- PyR@TE3 2-loop beta functions (SHA256 fingerprinted)
- αs(M_Z) = 0.1179 from `sm_inputs_mz.json`

**Key Results**:

| Fingerprint | Value | Scale |
|-------------|-------|-------|
| α₃(1 PeV) | 0.0524 | 10⁶ GeV |
| rel. dev to φ₀ | **1.38%** | — |
| α₃ = φ₀ crossing | — | **791,756 GeV** (~800 TeV) |
| α₃ = c₃ crossing | — | **2.16×10⁸ GeV** (~200 PeV) |

**Status**: All 4 PASS.

---

#### 10. `unification_gate` — Gauge Convergence Indicators

**What it tests**: Gauge coupling unification indicators.

**Key Results**: mismatch_rel ≈ 3.82×10⁻³ < tol=5×10⁻³

**Status**: PASS

---

#### 11. `msbar_matching_map` — PDG→MS̄ Boundary

**What it tests**: MS̄ matching at M_Z.

**Key Results**: g_i(M_Z), y_t(m_t), y_b(m_t)

**Status**: PASS

---

#### 12. `below_mt_eft_cascade` — QCD EFT Below m_t

**What it tests**: QCD EFT running below top threshold.

**Key Results**: αs(μ) at thresholds, m_b running

**Status**: PASS

---

#### 13. `stability_unitarity_audit` — λ(μ) Red-Flags

**What it tests**: Higgs quartic stability and unitarity bounds.

**Key Results**: λ=0 at μ≈1.48×10¹⁰ GeV; max|coupling|<4π

**Status**: PASS

---

#### 14. `matching_finite_pieces` — Finite Matching Pieces

**What it tests**: Finite matching contributions.

**Status**: PASS

---

### Gravity & Inflation

---

#### 15. `effective_action_r2` — Starobinsky R² Coefficient

**What it tests**: Heat-kernel computation of the R² coefficient from OperatorSpec.

**Inputs**: `effective_action_r2_operator_spec.json`

**Key Results**:

| Block | Rank | Statistics | β_R² contribution |
|-------|------|------------|-------------------|
| torsion_trace_vector | 4 | boson | 8.80×10⁷ |
| torsion_axial_vector | 4 | boson | 8.80×10⁷ |
| torsion_tensor_q | 16 | boson | 3.52×10⁸ |
| fp_ghost_vector | 4 | ghost | -3.40×10⁻⁴ |
| **Total** | — | — | **5.278×10⁸** |

- M/M_Pl = **1.2564942083×10⁻⁵**

**Starobinsky Predictions**:

| N | n_s | r | 10⁹·A_s |
|---|-----|---|---------|
| 55 | 0.9636 | 0.00397 | 2.016 |
| **56** | **0.9643** | **0.00383** | **2.090** |
| 57 | 0.9649 | 0.00369 | 2.166 |

**Status**: All 14 PASS. **Gap**: OperatorSpec is a contract, not BRST-derived.

---

### Cosmology

---

#### 16. `bounce_perturbations` — Cosmological Transfer Functions

**What it tests**: Scalar/tensor mode transfer through a non-singular bounce.

**Key Results**:

| Sector | k_bounce (raw) | T(k→∞) | max|T-1| (high k) |
|--------|----------------|--------|-------------------|
| Scalar | 2104 | 1.007 | 0.8% |
| Tensor | 7.8 | 1.004 | 0.7% |

- Wronskian Conservation: max|W/i-1| < 10⁻¹³

**Status**: All 6 PASS.

---

#### 17. `k_calibration` — Bounce Scale to CMB Multipoles

**What it tests**: Maps dimensionless bounce k to CMB multipoles ℓ.

**Key Results**:

| Mapping | ℓ_bounce_s | ℓ_bounce_t |
|---------|------------|------------|
| Naive | ~10⁵⁹ | ~10⁵⁶ |
| With budget (N=56) | ~5.2×10⁶ | ~1.9×10⁴ |

**Status**: All 6 PASS. **Gap**: a₀/a_transition not yet TFPT-threshold-derived.

---

#### 18. `cosmo_threshold_history` — Threshold-Derived Reheating

**What it tests**: Cosmological threshold history for reheating inputs.

**Status**: PASS

---

#### 19. `cosmo_reheating_policy_v106` — Reheating Window

**What it tests**: Reheating policy and a₀/a_transition.

**Status**: PASS

---

#### 20. `primordial_spectrum_builder` — Bounce Injection P(k)

**What it tests**: T(k) → P(k) table construction.

**Status**: PASS

---

#### 21. `boltzmann_transfer` — k̂→ℓ Mapping + CAMB

**What it tests**: Boltzmann transfer with CAMB C_ℓ backend.

**Key Features**: Signature-policy, pivot-normalization, optional Planck plik-lite logL

**Status**: PASS

---

#### 22. `bbn_neff_sanity` — BBN N_eff

**What it tests**: BBN effective neutrino number sanity check.

**Status**: PASS

---

#### 23. `bbn_consistency` — Light Elements Proxy

**What it tests**: BBN light element consistency.

**Status**: PASS

---

#### 24. `gw_background_bounds` — GW Bounds

**What it tests**: Gravitational wave background bounds.

**Status**: PASS

---

#### 25. `gw_background_predictor` — Ω_gw Baseline

**What it tests**: GW background prediction baseline.

**Status**: PASS

---

#### 26. `arrow_of_time_proxy` — Arrow Proxy

**What it tests**: Entropy-production proxy for arrow of time.

**Status**: PASS

---

#### 27. `arrow_mechanism` — Entropy-Production Proxy

**What it tests**: Arrow of time mechanism.

**Status**: PASS

---

### Flavor Sector

---

#### 28. `ckm_full_pipeline` — CKM Matrix Prediction

**What it tests**: CKM matrix from Z₃ Yukawa texture + RG evolution.

**Key Results**:
- χ²(m_t, boundary proxy) ≈ 49.58
- Unitarity: max|V†V-I| = 6.7×10⁻¹⁶

**Status**: All 5 PASS. **Gap**: scheme/scale + below-m_t running needed.

---

#### 29. `pmns_z3_breaking` — Neutrino Mixing Angles

**What it tests**: PMNS angles from Z₃ symmetry breaking.

**Key Results**:

| Angle | Value |
|-------|-------|
| θ₁₃ | **8.74°** |
| θ₁₂ (TM1) | **34.32°** |
| θ₂₃ (LO) | 45.00° |
| δ_CP (LO) | 90.00° |

**Status**: All 3 PASS. **Gap**: δ_CP not yet derived as unique topological output.

---

#### 30. `pmns_full_pipeline` — Full PMNS Prediction

**What it tests**: Complete PMNS at m_t and UV with κ EFT flow.

**Key Results (m_t)**:
- θ₁₂ ≈ 34.32°, θ₁₃ ≈ 8.74°, θ₂₃ ≈ 45.51°, δ_CP ≈ 90°
- Best diagnostic χ² ≈ 32.24
- κ decouples above M_NR3: max|κ|(μ_UV) ~ 2×10⁻¹⁶ GeV⁻¹

**Status**: PASS

---

#### 31. `pmns_mechanism_bridge` — κ→(y_N, M_NR)

**What it tests**: κ EFT to Seesaw parameter mapping.

**Key Results**: κ stability <0.001° angle shift

**Status**: PASS

---

#### 32. `seesaw_block` — Seesaw Anchor

**What it tests**: Seesaw mass scale anchor.

**Key Results**: m_ν₃ ≈ 0.046 eV, M_R anchor

**Status**: 3 PASS

---

#### 33. `mobius_cusp_classification` — Cusp Uniqueness

**What it tests**: Möbius cusp classification.

**Key Results**: {1/3, 2/3, 1} unique for all tested max_den

**Status**: 3 PASS

---

#### 34. `mobius_delta_calibration` — δ* vs δ_M

**What it tests**: Möbius δ calibration.

**Key Results**: 0.157% deviation

**Status**: 3 PASS

---

#### 35. `mobius_z3_yukawa_generator` — Z₃ Yukawa Texture

**What it tests**: Z₃ Yukawa texture generation.

**Key Results**: δ_used = δ*, diagnostic CKM χ²

**Status**: PASS

---

#### 36. `koide_constraints` — Koide Relation

**What it tests**: Koide mass formula constraints.

**Status**: PASS

---

#### 37. `flavor_topology_mapper` — Flavor Topology

**What it tests**: Flavor sector topology mapping.

**Status**: PASS

---

#### 38. `flavor_joint_objective_scan` — Joint Flavor Scan

**What it tests**: Joint CKM/PMNS objective optimization.

**Status**: PASS

---

### Mass Spectrum

---

#### 39. `mass_spectrum_minimal` — Mass Check

**What it tests**: Minimal mass spectrum validation.

**Status**: PASS

---

#### 40. `mass_spectrum_deriver` — Mass Derivation

**What it tests**: Mass spectrum derivation pipeline.

**Status**: PASS

---

### Torsion Sector

---

#### 41. `torsion_bounds_mapping` — Torsion Falsifiability

**What it tests**: TFPT torsion prediction vs SME experimental bounds.

**Key Results**:
- Chosen regime: |S_μ| ≈ 1.44×10⁻⁴² GeV
- Max ratio |S_μ^pred|/|S_μ^bound| ≈ 6.84×10⁻¹²

**Status**: All 6 PASS. **Gap**: current regime is toy benchmark.

---

#### 42. `torsion_condensate` — Λ/⟨K²⟩ Derivation

**What it tests**: Torsion condensate derivation candidate.

**Status**: PASS

---

#### 43. `torsion_observable_designer` — Torsion Observable Table

**What it tests**: Design table for torsion observables.

**Status**: PASS

---

#### 44. `torsion_falsifiability_snr` — Source+Noise SNR Gate

**What it tests**: Torsion SNR with explicit noise policy.

**Key Features**: Magnetar electron-polarization proxy

**Status**: PASS

---

#### 45. `torsion_observable_spin_fluid` — He-3 Torsion Benchmark

**What it tests**: Helium-3 spin-fluid torsion benchmark.

**Status**: PASS (INFO: not measurable under current assumptions)

---

#### 46. `torsion_dm_pipeline` — Optional DM Branch

**What it tests**: Torsion dark matter pipeline.

**Status**: PASS (not required under shipped axion policy)

---

### CMB & Birefringence

---

#### 47. `birefringence_tomography` — Cosmic Birefringence

**What it tests**: TFPT β_inf vs CMB birefringence measurements.

**Key Results**:
- TFPT β_inf = **0.2424°**
- ΔBIC(drift - step) = -2.19 → slight preference for drift

**Status**: 2 PASS

---

### Baryon Sector

---

#### 48. `omega_b_conjecture_scan` — Ω_b Identity

**What it tests**: Is Ω_b = (4π-1)·β_rad under topological assumptions?

**Key Results**:
- Ω_b_pred = **0.0489**
- Planck Ω_b = 0.0493 ± 0.0009
- z-score ≈ 0.42

**Status**: All 5 PASS. **Gap**: needs operator/anomaly-level derivation.

---

#### 49. `baryogenesis_placeholder` — Baryogenesis Gate

**What it tests**: Baryogenesis gate (delegates to mechanism).

**Key Results**: η_b match

**Status**: PASS

---

#### 50. `baryogenesis_mechanism` — Leptogenesis Proxy

**What it tests**: Leptogenesis mechanism proxy.

**Status**: PASS

---

### Dark Matter / Axion

---

#### 51. `axion_fa_derivation` — f_a Derivation

**What it tests**: Axion decay constant derivation.

**Status**: PASS

---

#### 52. `axion_dm_pipeline` — Axion Relic

**What it tests**: Post-inflation RMS + discrete strings/DW factor.

**Key Results**: C_str = 7/3 ⇒ Ωa h² ≈ 0.123

**Status**: PASS

---

#### 53. `axion_scenario_matrix` — Axion Scenario Viability

**What it tests**: Axion scenario viability matrix.

**Status**: PASS

---

#### 54. `dm_alternative_channels` — DM Closure Gate

**What it tests**: Alternative DM channel closure.

**Status**: PASS

---

### Dark Energy

---

#### 55. `dark_energy_paths` — Λ Targets

**What it tests**: Target ledger (Λ, ⟨K²⟩) + bridge to torsion_condensate.

**Status**: PASS (INFO-only)

---

### Topology

---

#### 56. `aps_eta_gluing` — APS Seam Check

**What it tests**: Atiyah-Patodi-Singer eta invariant seam check.

**Key Results**: SF=1 → Δ_Γ=2π

**Status**: 4 PASS

---

#### 57. `topology_phase_map` — Topology Phase Map

**What it tests**: Topological phase mapping.

**Status**: PASS

---

### Audit & Summary

---

#### 58. `anomaly_cancellation_audit` — SM Anomaly Sums

**What it tests**: Standard Model anomaly cancellation.

**Key Results**: All 0; SU(2) doublets=12 (even)

**Status**: 6 PASS

---

#### 59. `qft_completeness_ledger` — Checklist

**What it tests**: QFT block derivation status.

**Status**: PARTIAL

---

#### 60. `global_consistency_test` — Multi-Observable χ²

**What it tests**: Overall consistency against reference observables.

**Key Results**:

| Observable | Prediction | Reference | z-score |
|------------|------------|-----------|---------|
| α⁻¹(0) | 137.0359982 | 137.0359992 | **-46.8** |
| β_deg | 0.2424° | 0.35±0.14° | -0.77 |
| λ_Cabibbo | 0.2245 | 0.2243±0.0005 | +0.32 |
| n_s | 0.9643 | 0.9649±0.0042 | -0.15 |

- χ²_total (no α) = **0.80** (excellent!)

**Status**: 2 PASS, 1 expected FAIL (α within 5σ)

---

#### 61. `predictions_dashboard` — Compact Table

**What it tests**: All key predictions in compact form.

**Status**: 2 PASS

---

#### 62. `likelihood_engine` — Covariance Likelihood Framework

**What it tests**: Dataset spec + plugin stubs for likelihood.

**Status**: PASS

---

#### 63. `uncertainty_propagator` — Uncertainty Framework

**What it tests**: Uncertainty propagation infrastructure.

**Status**: PASS

---

### Precision QED

---

#### 64. `g2_and_lamb_shift_proxy` — g-2 Proxy

**What it tests**: Anomalous magnetic moment proxy.

**Status**: PASS

---

#### 65. `qed_anomalies_audit` — Precision-QED Audit

**What it tests**: Precision QED audit.

**Status**: PASS

---

#### 66. `discrete_complexity_minimizer` — MDL Search

**What it tests**: Minimum description length search.

**Key Results**: Recovers β, λ_C, M/M_Pl expressions

**Status**: PASS

---

## theoryv3 Suite – Analysis Modules (14)

The theoryv3 suite formalizes the "discrete building blocks" thesis with focused analysis modules.

---

#### tv3-1. `seed_invariants_audit`

**What it tests**: π → core invariants consistency.

**Methodology**: Verifies that core TFPT invariants derive from π-based formulas.

**Status**: PASS

---

#### tv3-2. `defect_partition_g5_audit`

**What it tests**: δ₂ from g=5 and α(0) closure.

**Methodology**: 
- Tests g=5 defect partition
- Verifies α(0) closure under defect series

**Plots**: `alpha_defect_series.png`, `g_negative_control.png`

**Status**: PASS

---

#### tv3-3. `alpha_backreaction_sensitivity_audit`

**What it tests**: k sweep (α sensitivity vs backreaction exponent).

**Methodology**: Scans backreaction parameter k ∈ [0, 3] and measures ppm deviation.

**Plot**: `alpha_backreaction_ppm.png`

**Status**: PASS

---

#### tv3-4. `g5_origin_audit`

**What it tests**: g=5 from SU(5) holonomy degeneracy.

**Methodology**: Single origin derivation for g=5.

**Plot**: `g5_origin.png`

**Status**: PASS

---

#### tv3-5. `dark_energy_exponential_audit`

**What it tests**: exp(−α⁻¹/2) suppression and ρ_Λ.

**Methodology**: Tests exponential suppression formula for dark energy scale.

**Plot**: `rho_lambda_candidates.png`

**Status**: PASS

---

#### tv3-6. `dark_energy_norm_half_origin_audit`

**What it tests**: n=1/2 from double-cover degree.

**Methodology**: Derives n=1/2 normalization from topological double-cover.

**Plot**: `dark_energy_norm_origin.png`

**Status**: PASS

---

#### tv3-7. `flavor_pattern_audit`

**What it tests**: λ, δ*, δ_CP, PMNS θ₁₃ checks.

**Methodology**: Validates Möbius + Z₃ flavor formulas.

**Plots**: `flavor_anchors.png`, `mobius_ratios.png`

**Status**: PASS

---

#### tv3-8. `pmns_tm1_audit`

**What it tests**: TM1 sum rule for sin²θ₁₂.

**Methodology**: Checks Tribimaximal-Reactor (TM1) sum rule compliance.

**Status**: PASS

---

#### tv3-9. `yukawa_exponent_index_audit`

**What it tests**: Rational indices for mass ratios.

**Methodology**: Validates rational q indices for Yukawa hierarchies.

**Plot**: `yukawa_q_errors.png`

**Status**: PASS

---

#### tv3-10. `yukawa_index_mapping_audit`

**What it tests**: Map q_ij to charge-squared index sums.

**Methodology**: Tests charge-squared index sum mapping.

**Plot**: `yukawa_index_mapping.png`

**Status**: PASS

---

#### tv3-11. `baryon_consistency_audit`

**What it tests**: Ω_b identity, η_b proxy, derived H₀.

**Methodology**: Cross-checks baryon sector identities.

**Plot**: `baryon_consistency.png`

**Status**: PASS

---

#### tv3-12. `axion_dm_audit`

**What it tests**: Axion frequency and relic fraction audit.

**Methodology**: Validates axion DM target parameters.

**Plot**: `axion_summary.png`

**Status**: PASS

---

#### tv3-13. `g5_crosslink_audit`

**What it tests**: g=5 consistency across sectors.

**Methodology**: Verifies g=5 appears consistently in α, flavor, and cosmology sectors.

**Plot**: `g5_links.png`

**Status**: PASS

---

#### tv3-14. `constant_factory_audit`

**What it tests**: Grouped constant ledger with derivations.

**Methodology**:
- Generates formula_math/formula_text for each constant
- Computes sensitivities
- Validates grammar and crosslinks
- Compares anchor-vs-TFPT views
- EW scale + lepton ladder marked pending

**Outputs**:
- `constant_factory_ledger.json`
- `constant_factory_ledger.tex`

**Plots**: `constant_factory_summary.png`, `constant_factory_sensitivity.png`

**Status**: PASS

---

## Unconventional Suite – Search/Audit Modules (7)

These are **search/audit/design tools** that accelerate closing ToE gaps. They are **not** publication-grade physics derivations.

---

#### ux-1. `ux_matching_metamorphic_audit`

**What it tests**: Metamorphic invariants of matching primitives in `tfpt_suite/matching.py`.

**Methodology**:
- Random samples (seeded)
- Tests: αs threshold up/down ≈ identity; gauge matching roundtrip
- Reports: max abs/rel errors, worst-case samples

**How it helps**: Engineering hardening. Catches silent regressions in scheme/threshold bookkeeping.

**Plot**: `matching_metamorphic_errors.png`

**Status**: PASS if errors < conservative tolerances.

---

#### ux-2. `ux_cosmo_history_sampler`

**What it tests**: k→ℓ bridge feasibility by sampling expansion-history parameters.

**Methodology**:
- Samples: N_inflation, N_reheat, T_reheat (log-uniform priors)
- Computes: a₀/a_transition, ℓ_bounce for each sample
- Reports: fraction with ℓ in CMB window [2, 2500]

**How it helps**: Policy sensitivity study. Determines if bounce features can be CMB-visible.

**Plots**: `cosmo_history_sampler_log10ell.png`, `cosmo_history_sampler_cmb_fraction.png`

---

#### ux-3. `ux_flavor_holdout_search`

**What it tests**: CKM convention choices with holdout validation.

**Methodology**:
- Grid search over discrete choices: delta_source, s13_mode, delta_mode
- Holdout split: fit on 7 CKM entries, evaluate on 2 (V_ub, V_td)
- No continuous parameter tuning

**How it helps**: Safeguard against convention shopping.

**Plot**: `flavor_holdout.png`

---

#### ux-4. `ux_gravity_gaugefix_ga`

**What it tests**: GA search for Laplace-type gauge-fixing choices.

**Methodology**:
- Toy operator family: Δ_μν = -g_μν□ + c_nonmin∇_μ∇_ν + Eg_μν
- Fitness: minimize |c_nonmin|
- Genetic algorithm with seeded randomness

**How it helps**: Scaffolding for BRST derivation.

**Plot**: `gravity_gaugefix_ga.png`

---

#### ux-5. `ux_omega_b_aps_bridge`

**What it tests**: APS seam term → Ω_b coefficient docking.

**Methodology**:
- Uses APS toy model from `aps_eta_gluing`
- Forms K_candidate = 2Δ_Γ - 1 = 4π - 1
- Computes Ω_b_pred = K·β_rad

**How it helps**: Concrete docking point for anomaly/inflow derivation.

**Plot**: `omega_b_aps_bridge.png`

---

#### ux-6. `ux_threshold_graph_audit`

**What it tests**: Declarative threshold graph visualization.

**Methodology**:
- Extracts nodes (thresholds) and edges (segments) from 2-loop runner
- Compares matching-disabled vs matching-enabled
- Records segment metadata

**How it helps**: Reviewer-proof audit trail for threshold handling.

**Plot**: `threshold_graph.png`

---

#### ux-7. `ux_torsion_regime_designer`

**What it tests**: Design tool for falsifiable torsion regimes.

**Methodology**:
- Toy spin-medium proxy: |S_μ| ~ c_spin·ρ_spin/M_eff²
- Computes c_spin_max = |S_μ_bound|/|S_μ_pred(c=1)|
- Proposes regimes compatible with `torsion_regimes.json`

**How it helps**: Proposes nontrivial regimes (magnetar, plasma) for falsification.

**Plot**: `torsion_regime_designer.png`

---

## Summary: Explicit Missing Bridges

The suite explicitly flags these publication-grade gaps:

### 1. Newton Scale / "G Derived"
- `effective_action_r2` is closure-consistent
- **Gap**: OperatorSpec is a contract, not BRST-derived from microscopic action

### 2. k→ℓ Bridge
- `k_calibration` + `cosmo_reheating_policy_v106` provide policy path
- **Gap**: a₀/a_transition not yet derived from TFPT threshold history

### 3. Flavor Publication-Grade Comparison
- CKM: χ² is diagnostic; scheme/scale and below-m_t policy needed
- PMNS: δ_CP is conventional; not yet derived as unique topological output

### 4. Ω_b Derivation
- Identity is conditional; **gap**: operator/anomaly-level derivation needed

### 5. Torsion Falsifiability
- Mapping exists; SNR bookkeeping exists (`torsion_falsifiability_snr`)
- **Gap**: replace σν proxies with publication-grade dataset likelihood

### 6. Dark Matter
- Axion pipeline shows θ_i=φ₀ underproduces DM
- **Gap**: ladder derivation + scenario closure needed

### 7. Dark Energy / Λ
- Target magnitudes explicit; **gap**: deriving ⟨K²⟩ or terminal-stage ladder

### 8. Renormalization/Matching
- Engineering closure exists
- **Gap**: finite pieces + unified covariance-level propagation

---

## Running the Suites

### TFPT Suite (Main)

```bash
# Engineering mode (default)
python tfpt-suite/run_suite.py run-all

# Physics mode (strict)
python tfpt-suite/run_suite.py run-all --physics

# Build PDF report
python tfpt-suite/run_suite.py build-report
```

### theoryv3 Suite

```bash
# Run all theoryv3 modules
python tfpt-suite/theoryv3/run_theoryv3.py run-all

# Build theoryv3 PDF
python tfpt-suite/theoryv3/run_theoryv3.py build-report

# Run tests
python tfpt-suite/theoryv3/run_theoryv3.py run-tests
```

### Unconventional Suite

```bash
# List modules
python tfpt-suite/unconventional/run_unconventional_suite.py list-modules

# Run all
python tfpt-suite/unconventional/run_unconventional_suite.py run-all
```

### Unit Tests

```bash
# Main suite tests
python -m unittest discover -s tfpt-suite/tests

# theoryv3 tests
python tfpt-suite/theoryv3/run_theoryv3.py run-tests
```

---

## Tracking

For implementation status and TODO checklist, see:
- `tfpt-suite/progress.md`
- `tfpt-suite/TOE_CLOSURE_STATUS.md`
- `tfpt-suite/MASTER_TOE_BACKLOG.md`

---

*Last updated: 2026-01-28*
