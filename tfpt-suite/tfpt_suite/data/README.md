# `tfpt_suite.data`

Data files used by `tfpt-suite` modules.

## Files

- `torsion_bounds_vetted.json`: vetted axial torsion bounds (Sun-centered A_μ components) used by `torsion_bounds_mapping` (Kostelecký–Russell–Tasson 2008, Eq. (6)).
- `torsion_bounds.json`: legacy/fallback bounds file (kept for backwards compatibility; currently matches the vetted table).
- `torsion_regimes.json`: explicit TFPT-side torsion prediction regimes used by `torsion_bounds_mapping` to make the module falsifiable (vacuum today vs nontrivial source/benchmark regimes).
- `torsion_falsifiability_noise_v1.json`: explicit source+noise policy for the torsion falsifiability SNR module (turns “measurable?” into a reproducible SNR calculation; includes a frequency-dependent PTA-style PSD proxy for magnetar timing noise).
- `global_reference.json`: vetted reference means/σ (with definitions/units/sources) for `global_consistency_test` (still no covariance; treat as a dashboard unless you supply a full likelihood).
- `global_reference_cov.json`: optional minimal covariance artifact (currently a conservative 2×2 α-sector example) used by `global_consistency_test` for an auxiliary covariance-based α gate.
- `global_reference_minimal.json`: minimal citeable external reference table (CODATA + Planck 2018 base-ΛCDM core parameters).
- `references.json`: consolidated reference ledger (dataset_id, version, value, sigma, units) used by theoryv3 and audits (includes CKM/PMNS anchors, alpha_s/sin2_thetaW at MZ, quark/lepton masses, and Planck/CODATA baselines).
- `likelihood_datasets_v1.json`: unified likelihood dataset schema (means + covariance/bounds + nuisance-policy) used by the likelihood-engine module(s) to move beyond independent-Gaussian scorecards. The `nuisance_policy.planck_2018` block declares fixed `A_planck` and pivot normalization.
- `bbn_reference.json`: reference table for BBN light-element checks (`bbn_consistency`: Y_p, D/H, N_eff; engineering-level falsifiability anchor until a publication-grade likelihood is added).
- `k_calibration.json`: assumptions + fallback bounce diagnostics used by `k_calibration` to map bounce k-scales to CMB multipoles (includes `assumptions.k_hat_grid_for_ell_range` for the ℓ-coverage diagnostic, independent of the transfer-function k-grid, now widened to cover tensor-CMB targets). Also carries a deterministic reheating policy block under `assumptions.reheating_policy_v106` and a CMB vs small-scale `assumptions.signature_policy`.
- `cosmo_reheating_policy_v106.json`: v1.06 reheating-window policy inputs (w_reh, g*, c_end, scan range) used by `cosmo_reheating_policy_v106`.
- `effective_action_r2_operator_spec.json`: machine-readable OperatorSpec contract consumed by `effective_action_r2` (torsion multi-block decomposition 4+4+16 plus a FP ghost block; generated deterministically from `microscopic_action_tfpt_v25.json` with action-term parse metadata).
- `microscopic_action_tfpt_v25.json`: canonical microscopic action/field-content spec (UFE action + SM sector) referenced by QFT-completeness modules; intended as the single source of truth for action-level inputs.
- `sm_inputs_mz.json`: minimal SM inputs at μ=MZ (including `alpha_s_sigma` for uncertainty propagation and basic mass thresholds `mc_GeV, mb_GeV, mt_GeV`) used to initialize RG dressing in `ckm_full_pipeline` and other RG-dressed modules.
- `lepton_masses_pdg.json`: lepton pole masses (e, μ, τ) used by flavor-anchor modules (e.g. Möbius δ calibration).
- `ckm_reference.json`: CKM |V_ij| table for `ckm_full_pipeline` (PDG 2024 global-fit values with explicit snapshot metadata and `reference.chi2_keys` subset for diagnostic χ²).
- `pmns_reference.json`: PMNS mixing reference values for `pmns_full_pipeline` (NuFIT 5.3 2024 with SK-atm; sin²θ12, sin²θ13, sin²θ23, δCP for NO/IO; used only for a diagnostic χ² after PMNS canonicalization).
- `nufit_pmns_grid_v53_sk_atm_NO.txt.xz`: optional NuFIT PMNS χ² grid table (flat columns for sin²θ12, sin²θ13, sin²θ23, δCP, Δm²21, Δm²3ℓ + χ²) consumed by the likelihood-engine plugin when enabled.
- `birefringence_tomography.json`: birefringence analysis config + shipped real CMB birefringence constraints (global rotation) for `birefringence_tomography`.
- `two_loop_rg_fingerprints.json`: explicit model-selection + fingerprint config for `two_loop_rg_fingerprints` (selects a `pyrate_pythonoutput_kind` from `pyrate_pythonoutputs.json` to avoid duplicating paths; fail-fast to prevent silent SM/E8 swaps; also carries the optional gravity α³ runner patch settings).
- `unification_gate_policy.json`: policy inputs for `unification_gate` (mt-boundary route, Δb3 discrete scan, optional 2-loop G8 patch assumption).
- `pyrate_pythonoutputs.json`: central registry of PyR@TE `PythonOutput/` directories + expected model names + YAML sources used by RG runners (`rge_pyrate_2loop` and pipelines).
- `alpha_running_pdg.json`: external SM/QED inputs for the secondary ᾱ^(5)(MZ) consistency check (alpha(0) → alpha(MZ) via Δα), including an explicit MSbar↔on-shell shift (PDG Eq. 10.13) and explicit top decoupling.
- `below_mz_policy.json`: explicit below-MZ EFT policy for QED running (EW integrated out at MZ; charged-fermion thresholds for α_em running).
- `flavor_texture_v24.json`: explicit conventions for the Z3 “circulant + diagonal” Yukawa texture (defines phase mode, δ source, and the fixed-coefficients rule used by the suite) plus a discrete `ckm_variants` list (used by `ckm_full_pipeline`) and a `neutrino_mechanism` block (used by `pmns_full_pipeline`). Also carries the `topology_phase_atoms` docking block referencing `chiral_index_three_cycles`; topology_phase_map now emits an explicit holonomy-class → (δ, δ_CP) candidate map used as a discrete filter. The `joint_objective_policy` section declares weights and mass-ratio penalty settings for `flavor_joint_objective_scan`.
- `chiral_index_three_cycles.json`: flux/holonomy sketch inputs (ν1,ν2,νT + example SM hypercharges) used by `chiral_index_three_cycles`.
- `axion_tfpt_v106.json`: TFPT axion numbers quoted in five_problems.tex (f_a, m_a) plus the explicit cosmology scenario policy used by `axion_dm_pipeline` / `axion_scenario_matrix` (pre- vs post-inflation) including a discrete post-inflation strings/domain-walls enhancement policy (`strings_domain_walls_factor_policy`, default \(C_{\rm str}=7/3\)). Also carries QCD inputs for the standard m_a formula.
- `rge_thresholds_v25.json`: threshold scales and scheme metadata for mt→UV RG pipelines (MSigma, MG8, MNR1..3, MPhi).

These are meant to be **editable inputs** (e.g. update with the latest SME tables / literature values).

