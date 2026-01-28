# tfpt_suite

Core library for the TFPT verification suite (modules, data helpers, and shared utilities).

## Recent Updates
- 2026-01-27: `likelihood_engine` now backfills alpha_bar5_inv_MZ from global_consistency terms and reads DM relic density from `relic_density` when present; `k_calibration` widens the ell-range diagnostic grid to cover tensor CMB targets.
- 2026-01-27: Added explicit EW finite-piece helpers in `matching.py` (Hempfling & Kniehl Δy_t and Buttazzo λ^(1) via A0/B0), plus optional finite-delta hooks in `match_yukawa` and `match_quartic`.
- 2026-01-27: OperatorSpec generation now parses the torsion-sector action term to derive block structure (action-parse metadata + block-source audit).
- 2026-01-27: Added `cosmo_threshold_history` (threshold-derived reheating) and wired `k_calibration` to prefer its outputs.
- 2026-01-27: Added optional Planck low-ℓ/lensing likelihood hooks in `boltzmann_transfer` and declared the CMB vs small-scale signature policy in `k_calibration.json`.
- 2026-01-27: `likelihood_engine` now supports an optional NuFIT PMNS grid plugin (flat χ² table + interpolation) when a grid file is provided.
- 2026-01-27: Added `axion_fa_derivation` (E8 ladder + block constants) and wired axion modules to prefer derived f_a/m_a outputs.
- 2026-01-27: `torsion_condensate` now uses a discrete gap equation from torsion-sector β_R2 plus spectral-flow quantization and reports a log10 ρ_Λ sigma policy.
- 2026-01-27: `dark_energy_paths` now identifies the E8 ladder terminal stage (n≈30) consistent with ρ_Λ and reports the φ_n sequence.
- 2026-01-27: `torsion_falsifiability_snr` now uses frequency-dependent PTA-style PSD noise inputs instead of σν-only proxies.
- 2026-01-27: `torsion_observable_spin_fluid` now declares a concrete He-3 experiment spec with required sensitivity.
- 2026-01-27: `arrow_mechanism` now includes a torsion-flux spectral-flow mechanism and a falsifiable prediction stub.
- 2026-01-27: `likelihood_engine` now reports a unified scorecard χ²/p-value combining α, flavor, CMB, and DM contributions.
- 2026-01-27: Added L1–L5 visualization plots (alpha defect z-score, gauge unification, CKM/PMNS residuals, k→ℓ feasibility, global radar).
- 2026-01-27: Defect-partition logic now enumerates two-defect sectors and exposes a g=5 justification chain used by the delta2 derivation module.
- 2026-01-27: Alpha running inputs now include explicit MSbar↔on-shell EW decoupling pieces (PDG Eq. 10.13) and top decoupling policy.
- 2026-01-27: Below-MZ QED policy moved to `tfpt_suite/data/below_mz_policy.json` and enforced in `below_mt_eft_cascade`.
- 2026-01-27: Implemented PDG Eq. 10.13 Δα̂−Δα formula in `matching.py` and surfaced it in `matching_finite_pieces`.
- 2026-01-27: Updated CKM/PMNS reference tables to PDG 2024 and NuFIT 5.3 (2024) snapshots; validated global_reference (BK18 r bound + snapshot dates).
- 2026-01-27: `topology_phase_map` now emits explicit holonomy-class → (δ, δ_CP) candidates with invariance checks.
- 2026-01-27: `flavor_joint_objective_scan` now includes a policy-driven mass-ratio penalty in the joint χ².
- 2026-01-27: CKM/PMNS pipelines now optionally consume `msbar_matching_map` mt-boundary outputs and report scale-policy checks.
# `tfpt_suite` package

Internal Python package for `tfpt-suite/`.

- `constants.py`: canonical TFPT invariants and derived quantities (paper v2.5)
- `conventions.py`: single source of truth for suite-wide conventions (e.g. hypercharge normalization g1_GUT vs gY) + conversion helpers
- `module_base.py`: standard module interface + output schema
- `io.py`: deterministic writers (JSON/TXT) and metadata
- `heat_kernel.py`: reusable heat-kernel a2 bookkeeping (Laplace-type operator specs)
- `sm_inputs.py`: SM MZ inputs + algebraic gauge-coupling initialization (no RG running)
- `rg_authority.py`: RG-authority policy loader + enforcement context (blocks legacy SM RG helpers above MZ except whitelisted EFT modules)
- `pyrate_boundary_runner.py`: authoritative mt boundary helper (short MZ→mt run using PyR@TE-generated betas; 2-loop + 1-loop diagnostic route)
- `rge_sm.py`: legacy SM RG helpers (kept for below‑MZ EFT audits and debugging; guarded by `rg_authority`)
- `rge_pyrate_2loop.py`: PyR@TE-driven 2-loop RGE engine (complex Yukawas + threshold bookkeeping + optional G8 Δb3 2-loop patch)
- `matching.py`: threshold matching primitives (minimal, explicit; placeholders for finite pieces beyond the implemented QCD α_s threshold matching)
- `pyrate_pythonoutputs.py`: central loader for `tfpt_suite/data/pyrate_pythonoutputs.json` (fail-fast registry of PyR@TE `PythonOutput/` locations + expected model names)
- `cosmo_scale_map.py`: explicit scale-factor mapping helpers used by `k_calibration` (a0/a_transition from inflation+reheating assumptions; k_hat→k(Mpc^-1)→ℓ mapping)
- `e8_ladder.py`: shared E8 cascade ladder helper (varphi_n) used by block-derivation modules
- `axion_inputs.py`: helper to resolve axion inputs (derived f_a preferred if available)
- `mobius_z3_yukawa_generator.py`: deterministic Möbius/Z3 quark Yukawa generator (v1.07SM-style CKM “cold pass”), including an explicit `QuarkPhaseMap` helper to keep \(δ_{\mathrm{Mobius}}\), \(θ_{\mathrm{texture}}\), and \(δ_{\mathrm{CKM}}\) distinct
- `suite_manifest.py`: scans `out/` + `out_physics/` and emits `suite_manifest.json` + paper-ready TeX summaries (single source of truth for module/check/plot counts)
- `export_reference_ledger.py`: exports `out/paper_reference_ledger.tex` from `data/references.json` (with stable citation keys)
- `modules/`: individual verification/proof modules

