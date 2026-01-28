# `tfpt_suite.modules`

This folder contains the **module implementations** referenced by `tfpt-suite/tasks.md`.

Each module:

- is deterministic given `SuiteConfig` (precision + seed)
- writes results under `tfpt-suite/out/<module_id>/`
  - `meta.json`
  - `results.json`
  - `report.txt`
  - `*.png` (optional plots; only when `SuiteConfig.plot=True`)

Module discovery is centralized in `registry.py`.

Notes:
- Plot generation is **best-effort** and must never make a module fail; plot errors are surfaced via `warnings` in `results.json`.
- The central PDF report builder embeds any `png/jpg/jpeg` found in each `out/<module_id>/` folder automatically.
- `unification_gate` reads `tfpt_suite/data/unification_gate_policy.json` (mt-boundary route + optional G8 Δb3 2-loop patch) to keep the unification scan discrete and reproducible.

## Implemented modules (quick list)

- `core_invariants`
- `ufe_gravity_normalization`
- `chiral_index_three_cycles`
- `two_loop_rg_fingerprints`
- `alpha_precision_audit`
- `predictions_dashboard`
- `qft_completeness_ledger`
- `anomaly_cancellation_audit`
- `effective_action_r2`
- `bounce_perturbations`
- `cosmo_reheating_policy_v106`
- `k_calibration`
- `primordial_spectrum_builder`
- `msbar_matching_map`
- `below_mt_eft_cascade`
- `stability_unitarity_audit`
- `aps_eta_gluing`
- `discrete_consistency_uniqueness`
- `mobius_cusp_classification`
- `mobius_delta_calibration`
- `seesaw_block`
- `birefringence_tomography`
- `omega_b_conjecture_scan`
- `axion_fa_derivation`
- `axion_dm_pipeline`
- `dark_energy_paths`
- `discrete_complexity_minimizer`
- `torsion_bounds_mapping`
- `global_consistency_test`
- `likelihood_engine`
- `ckm_full_pipeline`
- `pmns_full_pipeline`
- `pmns_mechanism_bridge`
- `pmns_z3_breaking`
- `brst_ghost_deriver`
- `mass_spectrum_minimal`
- `mass_spectrum_deriver`
- `koide_constraints`
- `bbn_neff_sanity`
- `bbn_consistency`
- `unification_gate`
- `cosmo_threshold_history`
- `matching_finite_pieces`
- `defect_partition_derivation`
- `alpha_on_shell_bridge`
- `topology_phase_map`
- `flavor_topology_mapper`
- `flavor_joint_objective_scan`
- `uncertainty_propagator`
- `axion_scenario_matrix`
- `gw_background_bounds`
- `gw_background_predictor`
- `boltzmann_transfer`
- `torsion_observable_spin_fluid`
- `torsion_observable_designer`
- `torsion_falsifiability_snr`
- `torsion_condensate`
- `torsion_dm_pipeline`
- `dm_alternative_channels`
- `baryogenesis_placeholder`
- `baryogenesis_mechanism`
- `g2_and_lamb_shift_proxy`
- `qed_anomalies_audit`
- `arrow_of_time_proxy`
- `arrow_mechanism`

Notes:
- `ufe_gravity_normalization` records the **Einstein-limit docking normalization** for κ/ξ from TFPT invariants (per eliminating_k.tex) and cross-checks that the closure-level BRST/ghost OperatorSpec chain is present (generated from `microscopic_action_tfpt_v25.json`). Deeper first-principles derivations from the full torsionful connection action remain future work, but are no longer represented as an unconditional WARN gate by default.
- `brst_ghost_deriver` now parses the torsion-sector action term to derive the block structure and verifies `operator_blocks_match_action`.
- `cosmo_threshold_history` derives T_reh/N_reh from the TFPT threshold ladder and writes a deterministic reheating summary consumed by `k_calibration` when available.
- `boltzmann_transfer` now supports optional Planck low-ℓ and lensing likelihoods (Cobaya, opt-in) and reports the CMB vs small-scale signature policy from `k_calibration.json` (the default policy is now **conditional** / gate-dependent, not an unconditional “tensor-CMB” claim).
- `axion_fa_derivation` derives f_a from the E8 ladder + block constants (n=10 PQ block), and downstream axion modules prefer the derived output when present.
- `torsion_condensate` now solves a discrete gap equation using β_R2 from the torsion operator spec and integer spectral-flow indices (APS η-gluing), with an explicit log10 ρ_Λ sigma policy.
- `dark_energy_paths` now identifies the ladder terminal stage by extrapolating the E8 φ_n sequence to n≈30 and records consistency with ρ_Λ.
- `torsion_falsifiability_snr` now consumes a frequency-dependent PTA-style PSD noise model (no longer σν-only proxies).
- `torsion_observable_spin_fluid` now declares a concrete He-3 comagnetometer experiment spec and required sensitivity (check `experiment_specified_with_sensitivity`).
- `arrow_mechanism` now includes a torsion-flux spectral-flow mechanism and a falsifiable prediction stub (`arrow_mechanism_non_invertible`, `arrow_prediction_falsifiable`).
- `likelihood_engine` now supports an optional NuFIT PMNS grid plugin (flat χ² table + interpolation) when a grid file is provided.
- `likelihood_engine` now reports a unified scorecard χ²/p-value combining α, flavor, CMB, and DM sector contributions.
- `likelihood_engine` now backfills alpha_bar5_inv_MZ from global_consistency terms and reads DM relic density from `axion_dm_pipeline` `relic_density` outputs.
- Plot updates: `alpha_precision_audit` adds `alpha_defect_zscore_overview.png`, `two_loop_rg_fingerprints` adds `gauge_unification_running.png`, `ckm_full_pipeline` adds `ckm_residuals_sigma.png`, `pmns_full_pipeline` adds `pmns_residuals_sigma.png`, `k_calibration` adds `k_to_ell_feasibility.png`, `boltzmann_transfer` adds `cl_tt_ee_comparison.png` and `cl_ratio_tt_ee.png`, and `global_consistency_test` adds `global_radar.png`.
- `chiral_index_three_cycles` makes the Appendix-J “3 boundary cycles → 3 families” mechanism explicit and exports a discrete **Wilson-line phase atom** set as a docking point for future topology→CP-phase derivations.
- `topology_phase_map` now emits explicit **holonomy-class → (δ, δ_CP)** candidate maps from Wilson-line atoms and includes invariance checks (cycle relabeling + complex conjugation), keeping the phase set finite and reproducible.
- `flavor_joint_objective_scan` now aggregates CKM + PMNS χ² with a declared mass-ratio penalty (per `flavor_texture_v24.json` policy) and verifies filter-only topology usage.
- `ckm_full_pipeline` and `pmns_full_pipeline` now optionally consume `msbar_matching_map` mt-boundary outputs and report explicit scale-policy checks (native mt vs reference).
- `defect_partition_derivation` now emits a two-defect sector enumeration (bound/separated/seam-coupled) and a g=5 justification chain as part of its results JSON.
- `alpha_on_shell_bridge` now consumes PDG 2024 inputs and an explicit MSbar↔on-shell EW decoupling shift (Δα̂−Δα) instead of a hidden remainder.
- `below_mt_eft_cascade` now consumes `tfpt_suite/data/below_mz_policy.json` for explicit below-MZ QED threshold policy.
- `matching_finite_pieces` now recomputes the PDG Eq. 10.13 Δα̂−Δα shift for the QED/EW finite-piece audit trail.
- `matching_finite_pieces` now also reports the 1-loop EW finite pieces for top Yukawa and Higgs quartic (Hempfling & Kniehl Eq. 2.6; Buttazzo App. A.1 via A0/B0).
- `cosmo_reheating_policy_v106` implements the v1.06 **ΔN reheating window** as a deterministic policy layer and emits a canonical \(a_0/a_t\) estimate under explicit assumptions.
- `msbar_matching_map` makes the **PDG-style inputs → MSbar boundary** bookkeeping explicit (gauge couplings at MZ and mt, yt(mt) from mt pole→MSbar(mt), derived yb/yτ under explicit approximations), reports deterministic αs(MZ) sensitivity, and (optionally) a deterministic Monte‑Carlo uncertainty hook + QCD‑focused IR diagnostics. This is the “boring but decisive” bridge needed for scheme/scale-consistent flavor comparisons.
- `below_mt_eft_cascade` is an explicit **audit trail** for the **below‑mt EFT cascade**: it reports αs across mc/mb thresholds (above/below) and a piecewise LO QCD running of mb(μ), making “nf steps + matching” assumptions inspectable.
- `stability_unitarity_audit` is a **red-flag detector** for vacuum stability (tracks λ(μ) and estimates a λ=0 crossing scale) and perturbativity (max coupling vs 4π), using the same 2-loop PyR@TE threshold engine as the flavor pipelines.
- `two_loop_rg_fingerprints` uses an explicit model-selection config (`tfpt_suite/data/two_loop_rg_fingerprints.json`) which selects a `pyrate_pythonoutput_kind` from the central registry (`tfpt_suite/data/pyrate_pythonoutputs.json`). It writes a fail-fast model fingerprint (expected model name + PythonOutput/YAML SHA256) to prevent silent SM↔E8 swaps in plots. The same config also carries an optional **gravity α³ runner patch** (κ-vector), and the report includes a baseline-vs-gravity diff summary when enabled.
- `ckm_full_pipeline` now performs **mt→μUV-only** RG evolution using **PyR@TE-generated 2-loop beta functions** integrated in-suite (so complex Yukawas are supported). It builds its mt-boundary quark Yukawas from the deterministic **Möbius/Z3 v1.07SM-style generator** (`tfpt_suite/mobius_z3_yukawa_generator.py`) with δ fixed from τ/μ and explicit CKM angle/phase conventions. The report prints an explicit **phase map** \(δ_{\mathrm{Mobius}}, θ_{\mathrm{texture}}, δ_{\mathrm{CKM}}\) to avoid convention drift. For its diagnostic χ² it uses a **minimal |V_ij| subset** declared by `tfpt_suite/data/ckm_reference.json/reference.chi2_keys` (default: `Vus,Vub,Vcb,Vtd`) to avoid double-counting unitarity constraints; the full |V_ij| matrix is still printed for transparency. If `tfpt_suite/data/flavor_texture_v24.json` sets `ckm_variants.phase_selection_rule_mode="filter_only"` and `topology_phase_map` outputs are available, CKM δ_CP is taken **only** from the finite topology candidate set and the gate `phase_set_derived_not_searched` is emitted as PASS.
- `pmns_full_pipeline` closes the “neutrino EFT flow” loop: it runs a coupled system with the Weinberg operator **κ** (1-loop EFT beta) plus 2-loop RG evolution for gauge/Yukawas, and implements stepwise threshold matching by activating \(N_{R1,2,3}\) at their scales (κ → 0 above \(M_{NR3}\)). It reports a diagnostic **χ² vs a reference table** (`tfpt_suite/data/pmns_reference.json`) to prevent “wrong labeling” artifacts. When `neutrino_mechanism.kappa_source="mechanism_bridge"`, the module additionally performs a **discrete** scan over topology-derived candidates (θ and `delta0`) and integer `eps` multipliers (configured in `tfpt_suite/data/flavor_texture_v24.json`). With `topology_phase_atoms.wiring.phase_selection_rule_mode="filter_only"` it fixes the PMNS convention via **mass-splitting canonicalization** (`pmns_convention_policy="mass_splitting_canonical"`) and emits the gate `no_convention_shopping_possible_under_fixed_rule` as PASS (no χ²-driven permutation choice).
- `pmns_mechanism_bridge` bridges the **phenomenology layer** (`pmns_z3_breaking`) to a **mechanism-layer** κ(mt): it reconstructs κ(mt) consistent with the selected Z3-breaking variant and a minimal neutrino spectrum choice, then factorizes κ into a consistent \((y_N, M_{NR})\) setup (tree level) and reports κ decoupling above \(M_{NR3}\).
- `mobius_delta_calibration` is the bridge from “Möbius deformation δ” in the older TFPT flavor notes to a deterministic numeric anchor (δ from τ/μ vs δ⋆ from varphi0).
- `seesaw_block` is an anchor-level block from paper v1.06 (MR scale and mν3 estimate); the full pipeline lives in `pmns_full_pipeline`.
- `k_calibration` maps the bounce scale to an approximate CMB multipole ℓ under explicit assumptions; it can now optionally **derive** \((N, N_{\rm reh})\) from a v1.06 reheating policy (see `assumptions.reheating_policy_v106` in `tfpt_suite/data/k_calibration.json`) and uses an expanded k_hat grid for the ℓ-coverage diagnostic.
- `effective_action_r2` supports a closure-level **gauge-parameter scan** diagnostic and an optional **SM matter a2 contributions** bookkeeping toggle via `TFPT_R2_INCLUDE_MATTER=1` (report shows raw sum + gravity-sector renormalization condition fixing TFPT’s \(M/\bar M_{\mathrm{Pl}}\)).
- `birefringence_tomography` now splits **data families**: CMB points at z≈1100 are treated as **offset checks** (no tomography), while step-vs-drift comparisons are performed on **low-z** points only.
- `bounce_perturbations` reports Wronskian conservation **for all k** and uses a per-k IC selector that avoids starting near ω²≈0 turning points (mitigates “cosmetic” Wronskian passes).
- `primordial_spectrum_builder` deterministically injects bounce transfer features into a primordial spectrum table \(P_\mathcal{R}(k), P_t(k)\) (baseline power law times \(|T|^2\)) using the `k_calibration` expansion-budget policy for \(a_0/a_{\rm tr}\). `boltzmann_transfer` consumes this table via CAMB’s `set_initial_power_table`.
- `boltzmann_transfer` can optionally evaluate the Planck 2018 high-ℓ plik-lite native likelihood (Cobaya) and report logL/Δχ² when `TFPT_ENABLE_PLANCK_LIKELIHOOD=1` is set.
- `axion_dm_pipeline` turns the TFPT axion numbers (five_problems.tex) into a **world-contact closure**: misalignment relic density and an explicit isocurvature tension check (engineering-level; flags missing ladder derivation + scenario dependence).
- `dark_energy_paths` turns “two TFPT-consistent paths” (five_problems.tex) into explicit **numerical targets** for \(\langle K^2\rangle\) (UFE condensate) and \(\varphi_\*\) (cascade magnitude), i.e. an engineering constraint module (not yet a prediction).
- Some modules remain intentionally **assumption-explicit scaffolds/proxies** (e.g. `dark_energy_paths` as a target ledger): remaining publication‑grade gaps are tracked as **INFO** (visible in the PDF/report text) without forcing the physics run yellow/red unless a concrete numeric gate is violated. Publication-grade follow-ups include **full external likelihood plugins** (Planck, NuFIT grid χ²), deeper connection-level BRST/operator derivations, and EW/QED matching with covariance propagation.

