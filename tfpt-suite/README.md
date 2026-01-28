# TFPT Suite (tfpt-suite)

A **modular, reproducible** Python verification suite for the Topological Field-Phase Theory (TFPT).

---

## Table of Contents

1. [Overview](#overview)
2. [High-Level Architecture](#high-level-architecture)
3. [Directory Structure](#directory-structure)
4. [Quickstart](#quickstart)
5. [Verification Modes](#verification-modes)
6. [Result PDFs Explained](#result-pdfs-explained)
7. [Module Categories](#module-categories)
8. [Output Structure](#output-structure)
9. [Interpreting Results](#interpreting-results)
10. [Related Documentation](#related-documentation)

---

## Overview

This suite turns theoretical claims from the TFPT paper into **deterministic, testable computations**. Every module:

- Takes explicit inputs (invariants, reference data, or upstream module outputs)
- Performs reproducible calculations with documented assumptions
- Emits machine-readable results (`results.json`), human-readable reports (`report.txt`), and optional plots
- Returns structured checks with `PASS/FAIL` status and severity levels

The goal is **reviewer-proof transparency**: every number, assumption, and approximation is traceable.

---

## High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           TFPT Suite                                      │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│  ┌──────────────────────────────────────────────────────────────────┐   │
│  │                    CONVENTIONAL SUITE                             │   │
│  │  tfpt-suite/tfpt_suite/modules/  (registry-driven; see `list-modules`)│   │
│  │                                                                    │   │
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐               │   │
│  │  │   Core      │  │   QFT/RG    │  │  Gravity/   │               │   │
│  │  │ Invariants  │  │  Matching   │  │  Cosmology  │               │   │
│  │  └─────────────┘  └─────────────┘  └─────────────┘               │   │
│  │                                                                    │   │
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐               │   │
│  │  │   Flavor    │  │  Torsion/   │  │  Dashboard  │               │   │
│  │  │  CKM/PMNS   │  │ Birefring.  │  │  & Global   │               │   │
│  │  └─────────────┘  └─────────────┘  └─────────────┘               │   │
│  │                                                                    │   │
│  │  Outputs: tfpt-suite/out/<module_id>/                              │   │
│  │  PDFs:    tfpt-test-results.pdf (engineering)                      │   │
│  │           tfpt-test-results-physics.pdf (strict)                   │   │
│  └──────────────────────────────────────────────────────────────────┘   │
│                                                                           │
│  ┌──────────────────────────────────────────────────────────────────┐   │
│  │                  UNCONVENTIONAL SUITE                             │   │
│  │  tfpt-suite/unconventional/  (registry-driven; see `list-modules`) │   │
│  │                                                                    │   │
│  │  Search / GA / Metamorphic / Holdout / Design Tools               │   │
│  │  Purpose: Accelerate closing explicit ToE gaps                    │   │
│  │  NOT: Publication-grade physics derivations                       │   │
│  │                                                                    │   │
│  │  Outputs: tfpt-suite/out/unconventional/<module_id>/               │   │
│  │  PDF:     unconventional/unconventional-test-results.pdf           │   │
│  └──────────────────────────────────────────────────────────────────┘   │
│                                                                           │
├─────────────────────────────────────────────────────────────────────────┤
│  Shared Infrastructure:                                                   │
│  - tfpt_suite/config.py (SuiteConfig: seed, precision, mode)             │
│  - tfpt_suite/module_base.py (TfptModule interface)                      │
│  - tfpt_suite/data/*.json (input tables, reference data)                 │
│  - Pyrate3/ (PyR@TE3 2-loop RGE outputs)                                 │
└─────────────────────────────────────────────────────────────────────────┘
```

### Key Design Principles

1. **Deterministic**: All randomness is seeded via `SuiteConfig.seed`
2. **Traceable**: Every module emits explicit assumptions and input sources
3. **Modular**: Modules can run independently; dependencies are explicit
4. **Two-layer verification**:
   - Unit tests (`tfpt-suite/tests/`) for fast regression checks
   - Module self-checks (`results.json → checks[]`) for physics-level validation

---

## Directory Structure

```
tfpt-suite/
├── run_suite.py                    # Main CLI runner
├── requirements.txt                # Python dependencies
├── tfpt-test-results.pdf           # Engineering-mode report
├── tfpt-test-results-physics.pdf   # Physics-mode report (stricter)
│
├── tfpt_suite/                     # Core library
│   ├── config.py                   # SuiteConfig (seed, precision, mode)
│   ├── module_base.py              # TfptModule interface
│   ├── report_builder.py           # PDF report generator
│   │
│   ├── modules/                    # 39 verification modules
│   │   ├── registry.py             # Module registry
│   │   ├── core_invariants.py      # c₃, φ₀, β, etc.
│   │   ├── alpha_precision_audit.py
│   │   ├── ... (37 more)
│   │   └── README.md
│   │
│   └── data/                       # Input tables (JSON)
│       ├── sm_inputs_mz.json       # SM parameters at M_Z
│       ├── global_reference.json   # Reference values for comparisons
│       ├── ckm_reference.json      # CKM matrix reference
│       ├── pmns_reference.json     # PMNS matrix reference
│       └── ... (15+ more)
│
├── out/                            # Module outputs (auto-generated)
│   ├── alpha_precision_audit/
│   │   ├── results.json
│   │   ├── report.txt
│   │   ├── meta.json
│   │   └── alpha_k_sensitivity.png
│   ├── ... (one folder per module)
│   └── unconventional/             # Unconventional module outputs
│
├── unconventional/                 # Auxiliary search/audit suite
│   ├── run_unconventional_suite.py # Unconventional CLI
│   ├── unconventional-test-results.pdf
│   ├── tfpt_unconventional/
│   │   └── modules/                # 7 unconventional modules
│   │       └── registry.py
│   └── docs/                       # Per-module documentation
│       └── ux_*/README.md
│
├── tests/                          # Python unit tests
│   ├── test_smoke.py
│   └── test_flavor_texture.py
│
├── theoryv3/                       # Analysis branch (discrete building blocks)
│   ├── run_theoryv3.py              # theoryv3 CLI runner
│   ├── theoryv3-analysis.pdf        # theoryv3 report
│   └── out/                         # theoryv3 outputs (results, plots)
│
├── TESTS_OVERVIEW.md               # Detailed test documentation
├── progress.md                     # Implementation status tracking
├── TOE_CLOSURE_STATUS.md           # Gap analysis + solution paths
└── tasks.md                        # Module requirements
```

---

## Quickstart

```bash
# Setup virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r tfpt-suite/requirements.txt

# Note: `boltzmann_transfer` uses CAMB (`camb` package) to compute CMB C_ℓ spectra.
# When `primordial_spectrum_builder` outputs injected primordial tables (bounce T(k) → P(k)),
# `boltzmann_transfer` consumes them via CAMB's `set_initial_power_table` (with an explicit pivot-normalization policy reported in the module output).
# Optional: set `TFPT_ENABLE_PLANCK_LIKELIHOOD=1` to evaluate Planck 2018 high-ℓ plik-lite + low-ℓ TT/EE + lensing likelihoods (Cobaya-native) and report logL/Δχ².
# If CAMB is unavailable on your platform, the module still runs but only emits the explicit k̂→ℓ mapping.

# List all available modules
python3 tfpt-suite/run_suite.py list-modules

# Run ALL modules (engineering mode, default)
python3 tfpt-suite/run_suite.py run-all

# Run in PHYSICS mode (stricter checks, more warnings)
python3 tfpt-suite/run_suite.py --mode physics run-all

# Run specific modules
python3 tfpt-suite/run_suite.py run --modules core_invariants,alpha_precision_audit

# Build PDF report from existing outputs
python3 tfpt-suite/run_suite.py build-report

# Build a machine-readable + TeX manifest from out/ + out_physics/
# (use --strict to fail on FAIL checks or plot issues)
python3 tfpt-suite/run_suite.py build-manifest --strict

# Export the paper-ready reference ledger table (TeX)
python3 tfpt-suite/run_suite.py export-reference-ledger

# Export the paper-ready prediction ledger (JSON + TeX)
python3 tfpt-suite/run_suite.py export-prediction-ledger
```

### Running the Unconventional Suite

```bash
# List unconventional modules
python3 tfpt-suite/unconventional/run_unconventional_suite.py list-modules

# Run all unconventional modules
python3 tfpt-suite/unconventional/run_unconventional_suite.py run-all

# Run specific unconventional module
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_flavor_holdout_search
```

### Running theoryv3

```bash
# Run theoryv3 analyses
python3 tfpt-suite/theoryv3/run_theoryv3.py run-all

# Build the theoryv3 PDF
python3 tfpt-suite/theoryv3/run_theoryv3.py build-report
```

---

## Verification Modes

The suite supports two verification modes:

### Engineering Mode (default)

```bash
python3 tfpt-suite/run_suite.py --mode engineering run-all
```

- Focus: Deterministic execution, explicit assumptions
- Check semantics: PASS/FAIL based on numerical stability and consistency
- Use case: Development, CI pipelines, baseline verification

### Physics Mode

```bash
python3 tfpt-suite/run_suite.py --mode physics run-all
```

- Focus: Strict physics validation
- Check semantics: Large deviations → WARN/FAIL; remaining publication-grade gaps are typically tracked as INFO (visible in reports without forcing yellow/red).
- Additional checks: α(0) metrology, CKM/PMNS χ² imports
- Current expectation (shipped policy): **no FAILs/WARNs** (PASS/INFO only). Any FAIL/WARN indicates a new regression or a deliberately tightened gate.
- Use case: Pre-publication review, catching "everything green" optics

---

## Result PDFs Explained

### `tfpt-test-results.pdf` (Engineering Mode)

The default report generated with `--mode engineering`. Contains:

- All conventional module outputs present under the selected output directory (module registry is discoverable via `list-modules`)
- Focus on **numerical consistency** and **deterministic execution**
- Check status reflects algorithmic/numerical success
- Best for: Development, debugging, CI validation

### `tfpt-test-results-physics.pdf` (Physics Mode)

Generated with `--mode physics`. Contains:

- Same conventional module set (registry-driven), **stricter interpretation**
- Deviations from experimental references → upgraded to WARN/FAIL
- Extended scorecard with α(0) metrology and flavor χ² imports
- Missing-derivation flags → explicit warnings
- Best for: Pre-submission audit, identifying real physics gaps

### `unconventional/unconventional-test-results.pdf`

Separate report for the unconventional modules (registry-driven; see `unconventional/run_unconventional_suite.py list-modules`):

- **ux_cosmo_history_sampler**: k→ℓ bridge feasibility via sampling
- **ux_flavor_holdout_search**: CKM convention search with holdout validation
- **ux_gravity_gaugefix_ga**: GA search for Laplace-type gauge-fixing choices
- **ux_matching_metamorphic_audit**: Metamorphic tests for matching primitives
- **ux_omega_b_aps_bridge**: APS seam → Ω_b coefficient docking
- **ux_threshold_graph_audit**: Declarative threshold graph visualization
- **ux_torsion_regime_designer**: Design tool for falsifiable torsion regimes

**Key distinction**: Unconventional modules are **search/audit/design tools**, not publication-grade physics derivations. They accelerate closing gaps but don't replace formal derivations.

---

## Module Categories

### 1. Core Invariants & Structural Numbers (5 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `core_invariants` | Fundamental TFPT constants | c₃=1/(8π), φ₀, β, M/M_Pl |
| `ufe_gravity_normalization` | Gravitational coupling normalization | ξ=c₃/φ₀ |
| `chiral_index_three_cycles` | 3-family chiral index | Ind(D)=3 |
| `discrete_complexity_minimizer` | MDL compression search | Minimal expressions |
| `discrete_consistency_uniqueness` | CFE uniqueness proof | α⁻¹ convergence |

### 2. α Sector & CFE (2 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `alpha_precision_audit` | α⁻¹ vs CODATA/PDG comparison | ppm deviation, k-sensitivity |
| `alpha_on_shell_bridge` | On-shell α connection | MS̄→pole bridge |

### 3. RG & Gauge Running (3 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `two_loop_rg_fingerprints` | 2-loop gauge running | α₃(μ)=φ₀ crossing scale |
| `unification_gate` | Gauge coupling convergence | Unification indicators |
| `stability_unitarity_audit` | Vacuum stability red-flags | λ=0 crossing, max couplings |

### 4. QFT Matching & EFT (4 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `msbar_matching_map` | PDG→MS̄ boundary conditions | g_i(M_Z), y_t(m_t) |
| `below_mt_eft_cascade` | QCD EFT below m_t | αs(μ) at thresholds |
| `anomaly_cancellation_audit` | SM anomaly freedom | All sums = 0 |
| `qft_completeness_ledger` | QFT checklist status | Block derivation status |

### 5. Gravity & R² Sector (2 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `effective_action_r2` | Heat-kernel R² coefficient | β_R², Starobinsky M/M_Pl |
| `aps_eta_gluing` | APS spectral flow check | Δ_Γ=2π seam term |

### 6. Cosmology & Bounce (6 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `bounce_perturbations` | Scalar/tensor transfer functions | T(k), Wronskian checks |
| `cosmo_reheating_policy_v106` | Reheating window | a₀/a_transition |
| `k_calibration` | k→ℓ multipole mapping | ℓ_bounce estimates |
| `primordial_spectrum_builder` | Bounce injection into P(k) | P_R(k), P_t(k) tables |
| `boltzmann_transfer` | CAMB C_ℓ backend | TT/TE/EE/BB spectra |
| `dark_energy_paths` | Λ target magnitudes | ⟨K²⟩, φ* targets |

### 7. Flavor Sector: CKM (4 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `mobius_cusp_classification` | Cusp set uniqueness | {1/3, 2/3, 1} |
| `mobius_delta_calibration` | δ* vs δ_M comparison | 0.16% deviation |
| `mobius_z3_yukawa_generator` | Z₃ Yukawa texture | Generated CKM |
| `ckm_full_pipeline` | RG-dressed CKM | χ² vs reference |

### 8. Flavor Sector: PMNS (4 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `pmns_z3_breaking` | Z₃ breaking angles | θ_ij, δ_CP variants |
| `pmns_mechanism_bridge` | κ→(y_N, M_NR) factorization | κ running stability |
| `pmns_full_pipeline` | Full PMNS at m_t/UV | PMNS angles, χ² |
| `seesaw_block` | Seesaw scale anchor | M_R, m_ν₃ order |

### 9. Torsion & Birefringence (2 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `torsion_bounds_mapping` | SME bounds → S_μ limits | Regime vs bounds ratio |
| `birefringence_tomography` | CMB birefringence fit | β_inf model comparison |

### 10. Baryon/DM/Axion (4 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `omega_b_conjecture_scan` | Ω_b=(4π-1)β identity | Coefficient z-score |
| `axion_dm_pipeline` | Axion DM relic abundance | Ω_a h², isocurvature |
| `axion_scenario_matrix` | Axion scenario viability | Pre/post-inflation |
| `baryogenesis_placeholder` | Placeholder for baryogenesis | Status marker |

### 11. Dashboard & Global Tests (4 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `predictions_dashboard` | Compact prediction table | All key predictions |
| `global_consistency_test` | Multi-observable χ² | Total χ², p-values |
| `likelihood_engine` | Covariance likelihood framework | logL from dataset spec |
| `mass_spectrum_minimal` | Mass spectrum check | Particle masses |

### 12. Auxiliary Modules (2 modules)

| Module | Purpose | Key Output |
|--------|---------|------------|
| `gw_background_bounds` | GW background limits | Stochastic bounds |
| `bbn_neff_sanity` | BBN N_eff sanity check | Radiation budget |

---

## Output Structure

Every module produces in `tfpt-suite/out/<module_id>/`:

| File | Content | Use |
|------|---------|-----|
| `results.json` | All numerical results | Machine parsing, downstream modules |
| `report.txt` | Human-readable summary | Quick inspection |
| `meta.json` | Runtime, seed, config | Reproducibility audit |
| `*.png` | Plots (optional) | Visual inspection |

Note: `spec.references` in `results.json` is auto-populated from the module's `results.references` block when present.

### results.json Structure

```json
{
  "module_id": "alpha_precision_audit",
  "title": "Alpha Precision Audit",
  "checks": [
    {
      "name": "primary_mz_within_1sigma",
      "passed": true,
      "severity": "PASS",
      "details": "z = -3.25e-04 ..."
    }
  ],
  "warnings": [],
  "numerical_results": {
    "alpha_inv_0_model2": 137.03599410,
    "ppm_vs_codata": -0.037
  },
  "assumptions": ["k=2 double-cover", "..."]
}
```

---

## Interpreting Results

### Check Status

| Status | Meaning |
|--------|---------|
| **PASS** | Validation successful |
| **FAIL** | Failed check (investigate) |
| **WARN** | Borderline (physics mode) |
| **INFO** | Informational only |

### Key Metrics

| Metric | Ideal | Good | Acceptable |
|--------|-------|------|------------|
| ppm deviation | < 0.1 | < 1 | < 10 |
| z-score | < 1 | < 2 | < 5 |
| χ²/DOF | ≈ 1 | < 2 | < 5 |
| Wronskian |W/i-1| | < 10⁻¹² | < 10⁻⁶ | < 10⁻³ |

### Current Headline Results

| Prediction | TFPT Value | Reference | Status |
|------------|-----------|-----------|--------|
| α⁻¹(0) [k=2] | 137.0359941 | 137.0359992 | -0.037 ppm |
| ᾱ⁵⁻¹(M_Z) | 127.951 | 127.951±0.009 | Excellent |
| β (birefringence) | 0.2424° | 0.35±0.14° | Within 1σ |
| λ_Cabibbo | 0.2245 | 0.2243±0.0005 | +0.32σ |
| n_s (N=56) | 0.9643 | 0.9649±0.0042 | -0.15σ |
| M/M_Pl | 1.256×10⁻⁵ | — | Starobinsky scale |

---

## Related Documentation

| Document | Purpose |
|----------|---------|
| `TESTS_OVERVIEW.md` | Detailed per-module documentation: inputs, theory, calculations, results, gaps |
| `progress.md` | Implementation status + Publication-grade ToE closure TODO |
| `missing.md` | Missing/placeholder data inputs + publication-grade gaps |
| `TOE_CLOSURE_STATUS.md` | Current results analysis + explicit remaining gaps + solution paths |
| `tasks.md` | Module requirements + proof style specifications |
| `tfpt_suite/modules/README.md` | Module implementation notes |
| `tfpt_suite/data/unification_gate_policy.json` | Unification gate policy (mt boundary + optional G8 Δb3 2-loop patch) |
| `unconventional/docs/ux_*/README.md` | Per-module unconventional documentation |

---

## Notes on Scope and "Proofs"

Several paper "missing proofs" are implemented as **verification modules**:

- Heat-kernel operator evaluation → `effective_action_r2`
- Numerical transfer functions → `bounce_perturbations`
- k→ℓ calibration → `k_calibration` + `cosmo_reheating_policy_v106`
- Bounce→CMB wiring (P(k) injection + CAMB C_ℓ) → `primordial_spectrum_builder` + `boltzmann_transfer`

For two-loop running: `two_loop_rg_fingerprints` **generates** gauge running from PyR@TE3 outputs under `Pyrate3/` with explicit model selection + SHA256 fingerprinting (fail-fast against silent model swaps).

This suite implements computable pipelines following the paper's stated assumptions (K1–K4 etc.) and produces:

- Explicit numerical outputs tied to TFPT invariants (c₃, φ₀)
- Internal consistency checks (limits, invariants, stability)
- Certificates/search logs for "uniqueness within bounded rule sets"

---

## Recent Updates

- 2026-01-27: Likelihood engine now backfills alpha_bar5_inv_MZ from global_consistency terms and reads DM relic density from axion outputs; k_calibration widens the ell-range diagnostic grid for tensor CMB coverage.
- 2026-01-27: Updated `TOE_CLOSURE_STATUS.md` physics snapshot values (CKM χ²/p and physics score p) to match `out_physics` artifacts.
- 2026-01-27: Added two-defect sector enumeration + g=5 justification chain to `defect_partition_derivation`.
- 2026-01-27: Updated α on-shell bridge to use PDG 2024 inputs with explicit MSbar↔on-shell EW decoupling pieces.
- 2026-01-27: Added explicit below-MZ QED policy (`below_mz_policy.json`) and wired it into `below_mt_eft_cascade`.
- 2026-01-27: Implemented PDG Eq. 10.13 Δα̂−Δα formula in `matching.py` and surfaced it in `matching_finite_pieces`.
- 2026-01-27: Added explicit EW finite pieces for top Yukawa and Higgs quartic (Hempfling & Kniehl Eq. 2.6; Buttazzo App. A.1) and reported them in `matching_finite_pieces`.
- 2026-01-27: OperatorSpec generation now parses the torsion-sector action term and `brst_ghost_deriver` verifies action-derived block matching.
- 2026-01-27: Added threshold-driven reheating module (`cosmo_threshold_history`) and wired `k_calibration` to use its derived T_reh/N_reh.
- 2026-01-27: Added optional Planck low-ℓ/lensing likelihood hooks in `boltzmann_transfer` and declared the CMB vs small-scale signature policy in `k_calibration.json`.
- 2026-01-27: `likelihood_engine` now supports an optional NuFIT PMNS grid plugin (flat χ² table + interpolation) when a grid file is provided.
- 2026-01-27: Added `axion_fa_derivation` to derive f_a from the E8 ladder/block constants and wired axion modules to prefer derived outputs.
- 2026-01-27: `torsion_condensate` now uses a discrete gap equation from torsion-sector β_R2 plus spectral-flow quantization and reports a log10 ρ_Λ sigma policy.
- 2026-01-27: `dark_energy_paths` now identifies the E8 ladder terminal stage (n≈30) consistent with ρ_Λ and reports the φ_n sequence.
- 2026-01-27: `torsion_falsifiability_snr` now uses frequency-dependent PTA-style PSD noise inputs instead of σν-only proxies.
- 2026-01-27: `torsion_observable_spin_fluid` now declares a concrete He-3 experiment spec with required sensitivity.
- 2026-01-27: `arrow_mechanism` now includes a torsion-flux spectral-flow mechanism and a falsifiable prediction stub.
- 2026-01-27: `likelihood_engine` now reports a unified scorecard χ²/p-value combining α, flavor, CMB, and DM contributions.
- 2026-01-27: Added L1–L5 visualization plots (alpha defect z-score, gauge unification, CKM/PMNS residuals, k→ℓ feasibility, global radar).
- 2026-01-27: Updated CKM/PMNS reference tables to PDG 2024 and NuFIT 5.3 (2024) snapshots; validated global_reference (BK18 r bound + snapshot dates).
- 2026-01-27: Added `theoryv3/` analysis branch with discrete-invariant audits, plots, and a dedicated PDF report.
- 2026-01-27: Added `tfpt_suite/data/references.json` reference ledger and `reference_ledger.py` loader for versioned datasets.
- 2026-01-27: `topology_phase_map` now emits explicit holonomy-class → (δ, δ_CP) candidates with invariance checks.
- 2026-01-27: `flavor_joint_objective_scan` now includes a policy-driven mass-ratio penalty in the joint χ².
- 2026-01-27: CKM/PMNS pipelines now optionally consume `msbar_matching_map` mt-boundary outputs and report scale-policy checks.

---

*Last updated: 2026-01-27*
