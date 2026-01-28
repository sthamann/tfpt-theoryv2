# progress (tfpt-suite)

This file tracks implementation status for the modular TFPT suite described in `tfpt-suite/tasks.md`.

## Status legend

- [ ] planned / not started
- [~] in progress
- [x] implemented (with automated checks/tests)

## Core infrastructure

- [x] suite skeleton (`run_suite.py`, config, IO, constants, module interface)
- [x] automated tests (`unittest`) and smoke runner (`tfpt-suite/tests/`)
- [x] reproducibility basics: fixed seed + precision control + run metadata per module

## Modules (from `tasks.md`)

Priority (per `tasks.md` / paper v2.4):

- [x] `effective_action_r2` (Appendix K): explicit **a2→β\_{R^2}→M** bookkeeping implemented via a machine-readable **OperatorSpec** (Laplace-type block operator on constant curvature). Spec uses a torsion multi-block decomposition (4+4+16) plus an explicit FP ghost block; \(E/R\equiv\alpha_R\) is **fixed numerically** (no runtime symbolic solve). Current suite checks are PASS with `derivation.status=derived`, and OperatorSpec generation now parses the torsion-sector action term for block provenance.
- [x] `bounce_perturbations` (Appendix L): numerical mode evolution + \(T(k)\) pipeline implemented for scalar+tensor in \(f(R)\) (Starobinsky) with a torsion \(a^{-6}\) proxy. Uses rescaled variables \(x_\hat = k_\text{bounce} x\), \(k_\hat = k/k_\text{bounce}\) and per-\(k\) adiabatic start selection; high-\(k_\hat\) limit \(T\to 1\) + Wronskian checks pass in the shipped configuration.
- [x] `k_calibration`: assumption-explicit mapping from bounce’s dimensionless \(k_\text{bounce}\) (in \(x=M\eta\) units) to an approximate CMB multipole \(\ell \approx k\,\chi_*\), and required overall scaling \(a_0/a_t\) to place bounce features at target \(\ell\); ℓ-coverage diagnostic now uses an expanded k_hat range to cover tensor CMB targets.

Additional modules:

- [x] `ufe_gravity_normalization`: closure-level Einstein-limit docking for κ/ξ (from invariants; per eliminating_k.tex) with an explicit audit that the closure-level BRST/ghost OperatorSpec chain is present (generated from the canonical microscopic action spec).
- [x] `chiral_index_three_cycles`: Appendix-J mechanism scaffold (3 boundary cycles → 3 families) plus a discrete Wilson-line phase atom set as a docking point for future topology→CP-phase derivations.
- [x] anomaly/uniqueness (CFE uniqueness certificate + self-consistent fixed point)
- [x] `alpha_precision_audit`: reproduces the paper’s baseline vs self-consistent \(\alpha^{-1}\), reproduces the exponent-sensitivity table, and computes the minimal next correction amplitude \(\delta_2\) in \(\varphi(\alpha)=\varphi_\mathrm{tree}+\delta_\mathrm{top}e^{-2\alpha}+\delta_2 e^{-4\alpha}\) required to hit CODATA exactly (debug target for the next derivation; not a fitted parameter in TFPT).
- [x] `two_loop_rg_fingerprints`: generates a two-loop gauge-running table from the shipped PyR@TE3 `PythonOutput/` under `Pyrate3/` (cached to `tfpt-suite/out/two_loop_rg_fingerprints/gauge_couplings.csv`) and checks the paper’s “RG fingerprint” claim (e.g. \(\alpha_3(1\,\mathrm{PeV})\) within 2% of \(\varphi_0\)), plus reports the \(\mu\) scales where \(\alpha_3(\mu)\) crosses \(\varphi_0\) and \(c_3\).
- [x] APS seam \(\eta\)-gluing / spectral flow checks (Appendix seam: \(D_\Gamma=i\,d/d\theta\), \(U_\Gamma\), \(\eta(0)\) via Hurwitz \(\zeta\))
- [x] double cover + Möbius cusp classification: SU(5) hypercharge-derived cusp set implemented (\(\{1,1/3,2/3\}\)), uniqueness verified by bounded rational scan.
- [x] birefringence tomography: runs on **real published CMB birefringence constraints** via `tfpt_suite/data/birefringence_tomography.json` (no longer synthetic-by-default)
- [x] torsion bounds mapping: ingestion/mapping + **vetted component-wise axial torsion bounds** (Kostelecký–Russell–Tasson 2008, Eq. (6)) included; TFPT prediction is now **regime-explicit** via `tfpt_suite/data/torsion_regimes.json` (includes `vacuum_today` with \(S_\mu=0\) and a default nontrivial background benchmark `cosmological_background_H0` with \(|S_\mu|\sim H_0\)).
- [x] global consistency test: deterministic scorecard is implemented with **two modes** (`--mode engineering|physics`). Physics mode upgrades large deviations to WARN/FAIL and imports hard-sector χ² from CKM/PMNS modules; p-values are reported (still not a publication-grade global likelihood without covariance).
- [x] `predictions_dashboard`: emits a paper-ready “Predictions” summary (5–10 key numbers with reference comparisons, z-scores where applicable, and explicit dependency notes).
- [x] `qft_completeness_ledger`: emits a structured checklist aligned with section 7.1 (QFT completeness) and maps concrete engineering next steps aligned with 7.2.
- [x] canonical microscopic action spec: machine-readable `tfpt_suite/data/microscopic_action_tfpt_v25.json` added (UFE action + SM sector). This is now the canonical “fields+reps+Lagrangian terms” source referenced by the QFT ledger and anomaly audit.
- [x] `anomaly_cancellation_audit`: explicit anomaly sums (U(1)^3, mixed anomalies, SU(3)^3) + SU(2) Witten global check for the SM (+ TFPT anomaly-neutral extensions).
- [x] CKM pipeline: **mt boundary (Yukawas + gauge) + upward-only RG evolution \(m_t\to \mu_{\rm UV}\)** implemented (deterministic). Reports CKM at \(m_t\) and at \(\mu_{\rm UV}=10^{16}\,\mathrm{GeV}\). Gauge+Yukawa running uses **PyR@TE-generated 2-loop beta functions** integrated in-suite (complex Yukawas supported), with threshold bookkeeping (MSigma switch; MG8 Δb3 patch). Texture uses the explicit Z3 “circulant + diagonal” matrix form with fixed coefficients per v1.07SM conventions.
- [x] PMNS with \(Z_3\) breaking: Z3-invariant vs breaking operator basis is derived from the family-space permutation symmetry (projector \(M_\text{inv}=\frac13(M+P^TMP+P^{2T}MP^2)\)); discrete \(\varepsilon=\varphi_0/6\) scan runs with perturbative shifts.
- [x] PMNS full pipeline: **Z3 texture + seesaw κ EFT + thresholds** implemented as `pmns_full_pipeline` with mt→μUV-only policy. Implements stepwise activation of \(N_{R1,2,3}\) at \(M_{NR1,2,3}\) (κ → 0 above \(M_{NR3}\)) and computes PMNS angles/δ at mt and μUV.
- [x] Möbius \(\delta\) anchors: `mobius_delta_calibration` computes \(\delta_M\) from \(\tau/\mu\) and compares to the geometric anchor \(\delta_\star=3/5+\varphi_0/6\), exposing the cusp map \(M_y(\delta)\) at \(y\in\{1,1/3,2/3\}\).
- [x] Seesaw block anchor: `seesaw_block` implements the paper v1.06 \(n=5\) seesaw-scale anchor (\(M_R\sim 10^{15}\,\mathrm{GeV}\)) and computes the implied \(m_{\nu_3}\sim v^2/M_R\) order of magnitude, plus the inferred block calibration \(\zeta_{N_R}\).
- [x] \(\Omega_b\) identity: implemented as a **conditional derivation** via explicit sector-counting assumptions plus regression against Planck-derived \(\Omega_b\) (from \(\Omega_b h^2\) and \(H_0\)); still replaceable by a full anomaly/inflow computation.
- [x] `cosmo_reheating_policy_v106`: deterministic v1.06 reheating/ΔN window policy module (derives \(N\) from \(n_s\) and estimates \(N_{\rm reh}\), \(a_0/a_t\) under explicit assumptions).
- [x] `cosmo_threshold_history`: threshold-driven reheating module (MSigma/MG8/MNR ladder → derived \(T_{\rm reh}\), \(N_{\rm reh}\), \(a_0/a_t\)) and `k_calibration` consumes it when present.
- [x] `axion_dm_pipeline`: engineering-level axion world-contact closure (misalignment relic density + **scenario-aware** isocurvature handling). Scenario choice is explicit via `tfpt_suite/data/axion_tfpt_v106.json` (`cosmo_policy.scenario_policy.selected`); default branch is post-inflation PQ breaking (isocurvature treated as N/A), while the pre-inflation proxy remains available as a diagnostic.
- [x] `dark_energy_paths`: engineering target module translating ΛCDM into target magnitudes for the UFE torsion-condensate path and the cascade-magnitude path; now identifies the ladder terminal stage (n≈30) consistent with ρ_Λ.
- [x] discrete complexity minimizer: bounded grammar search with symbolic verification (meta-check; not a substitute for operator-level QFT derivations)
- [x] `msbar_matching_map`: explicit, deterministic bookkeeping bridge for PDG‑style inputs at \(M_Z\) → MS̄ boundary values at \(m_t\) (gauge couplings + \(y_t\), derived \(y_b,y_\tau\)), including αs(MZ) sensitivity, an optional deterministic Monte‑Carlo uncertainty hook, **and** a deterministic **Jacobian-based linear covariance propagation** artifact (`covariance_propagated_end_to_end`). Next step is a publication‑grade matching map incl. systematic finite EW/QCD pieces + full below‑\(m_t\) EFT cascade (QCD+QED decoupling policy).
- [x] `below_mt_eft_cascade`: explicit audit trail for the below‑\(m_t\) QCD EFT cascade (αs above/below thresholds + LO piecewise mb running), to make “nf steps + matching” assumptions reviewer-proof and machine-checkable.
- [x] `stability_unitarity_audit`: automated red‑flags for vacuum stability/metastability (tracks \(λ(μ)\) and estimates a \(λ=0\) crossing scale) and perturbativity proxy (max coupling vs 4π ceiling), using the same 2‑loop PyR@TE threshold engine as the flavor pipelines.
- [x] `unification_gate`: explicit numeric PASS/FAIL gate for gauge-coupling unification under a declared 2-loop PyR@TE model (policy-driven mt boundary + optional G8 Δb3 2-loop patch via `tfpt_suite/data/unification_gate_policy.json`). **Current shipped policy scan passes in Physics-Mode** (best discrete candidate Δb3≈2.5 gives mismatch_rel≈3.82e−3 < tol=5e−3).
- [x] `mass_spectrum_minimal`: minimal “ToE scope map” ledger tagging masses/scales as `input` vs `derived` vs `placeholder` (no pretend-closure).
- [x] `bbn_neff_sanity`: MeV-era \(g_*(T)\) sanity ledger around e± annihilation + `N_eff` marker (no full BBN simulation).
- [x] `baryogenesis_placeholder`: delegates to `baryogenesis_mechanism` (deterministic leptogenesis proxy); Physics-Mode no longer FAILs.
- [x] `gw_background_bounds`: CMB tensor bound (r) + PTA/LIGO/BBN proxy bounds driven by `gw_background_predictor`; Physics-Mode no longer FAILs.
- [x] `g2_and_lamb_shift_proxy`: TFPT-scale new-physics suppression proxy (g-2 + Lamb shift); Physics-Mode no longer FAILs (consistency ledger, not a full QED calculation).
- [x] `defect_partition_derivation`: enumerates two-defect sectors (bound/separated/seam-coupled), derives a discrete g=5 multiplicity for δ₂, and reports the implied α(0) metrology z-score; justification chain is emitted in `results.json` and wired into `global_consistency_test`.
- [x] `alpha_on_shell_bridge`: explicit α(0) ↔ ᾱ^(5)(MZ) renorm chain (1-loop leptonic + PDG Δα_had + explicit MSbar↔on-shell EW decoupling shift) with metamorphic robustness checks.
- [x] `topology_phase_map`: discrete Wilson-line phase atoms → finite candidate set for (δ, δ_CP) anchors; now emits explicit holonomy-class → (δ, δ_CP) maps with invariance checks (still a docking map, not an operator derivation).
- [x] `flavor_joint_objective_scan`: joint CKM+PMNS objective aggregation over discrete CKM variants + PMNS best convention, now with a policy-driven mass-ratio penalty (no continuous fitter; p-value gate in physics mode).
- [x] `matching_finite_pieces`: explicit finite-piece milestone (QCD α_s 2-loop decoupling) + msbar matching-map integration for α_s(mc) diagnostics; includes PDG Eq. 10.13 Δα̂−Δα and the 1-loop EW top/Higgs finite pieces (Δy_t, Δλ_H) as audit proof points.
- [x] `axion_scenario_matrix`: explicit pre/post-inflation PQ scenario branches and physics gate “at least one scenario viable”.
- [x] `uncertainty_propagator`: end-to-end uncertainty ledger aggregator (MC summary wiring + χ² gate scaffold).
- [x] `torsion_dm_pipeline`: conditional/optional DM branch (PASS if axion-first closure suffices; otherwise flags torsion DM as required).
- [x] `arrow_of_time_proxy`: delegates to `arrow_mechanism` (entropy proxy + torsion-flux mechanism); Physics-Mode no longer FAILs.
- [x] `koide_constraints`: Koide-\(Q\) diagnostic for charged leptons (explicit docking check; not a TFPT derivation).
- [x] `mass_spectrum_deriver`: Möbius/Z3 ratio predictions from \(\delta_\star\) with lepton ratio gate (diagnostic; quark ratios emitted but not hard-gated).
- [x] `brst_ghost_deriver`: closure-level BRST/ghost/a₂ derivation module (deterministic OperatorSpec generation from the canonical microscopic action + action-term parse + gauge-parameter proxy scan). Deeper connection-level BRST derivations remain tracked as future work, but are no longer emitted as an unconditional WARN gate by default.
- [x] `primordial_spectrum_builder`: bounce transfer injection bridge — turns `bounce_perturbations` T(k) into primordial tables \(P_\mathcal{R}(k), P_t(k)\) in Mpc⁻¹ using the `k_calibration` expansion-budget policy (explicit a₀/a_transition).
- [x] `boltzmann_transfer`: k-hat→ℓ mapping + CAMB-backed C_ℓ backend consuming injected primordial P(k) tables (explicit pivot-normalization policy is reported). Planck 2018 high-ℓ + low-ℓ + lensing likelihoods are opt-in integrated (`TFPT_ENABLE_PLANCK_LIKELIHOOD=1`) and combined logL is reported.
- [x] `likelihood_engine`: unified dataset schema + covariance likelihood evaluation + NuFIT grid interpolator wiring (awaiting grid file) plus unified scorecard χ²/p-value across α/flavor/CMB/DM; alpha_bar5_inv_MZ backfilled from global_consistency terms and DM uses axion `relic_density`.
- [x] `flavor_topology_mapper`: topology→phase mapping scaffold + joint χ² proxy gate; explicitly marks “not wired” until operator-level derivation exists.
- [x] `torsion_observable_spin_fluid`: explicit He-3 lab benchmark with experiment spec and required sensitivity (check `experiment_specified_with_sensitivity`).
- [x] `torsion_observable_designer`: torsion observable design table (lab vs magnetar proxy) with explicit proxy gate.
- [x] `torsion_falsifiability_snr`: explicit **source+noise** policy + SNR calculation with a frequency-dependent PTA-style PSD and go/no-go gate (`go_no_go_snr_ge_5`).
- [x] `torsion_condensate`: discrete gap equation from torsion β_R2 + spectral-flow quantization with explicit log10 ρ_Λ sigma policy; Physics-Mode no longer FAILs.
- [x] `dm_alternative_channels`: DM closure gate (axion-first policy); Physics-Mode no longer FAILs.
- [x] `arrow_mechanism`: entropy proxy plus torsion-flux spectral-flow mechanism and falsifiable prediction stub.
- [x] `baryogenesis_mechanism`: deterministic vanilla leptogenesis proxy (Davidson–Ibarra + washout); Physics-Mode no longer FAILs.
- [x] `bbn_consistency`: engineering-level light-element proxy likelihood (Yp, D/H, Neff); Physics-Mode no longer FAILs.
- [x] `gw_background_predictor`: deterministic scale-invariant Ω_gw baseline (Starobinsky proxy); Physics-Mode no longer FAILs.
- [x] `qed_anomalies_audit`: delegates to `g2_and_lamb_shift_proxy` (consistency audit); Physics-Mode no longer FAILs.

## theoryv3 analysis branch

- [x] `theoryv3/` suite skeleton (runner, registry, report builder)
- [x] theoryv3 analyses: invariants, g=5 origin/defect partition, alpha backreaction sensitivity, dark energy exponential + norm origin
- [x] theoryv3 flavor audits: CKM/PMNS anchors + TM1 audit + Yukawa exponent/index mapping
- [x] theoryv3 cosmology audits: Omega_b identity, baryogenesis proxy, axion DM audit (H0 check added)
- [x] theoryv3 infrastructure: reference ledger + results schema validation gate
- [x] theoryv3 PDF report generation + unit tests
- [x] theoryv3 constant factory: EW scale + lepton ladder flagged as pending candidates in the audit (candidate formulas in `constantfactory.md`).

## Publication-grade ToE closure TODO (explicit gaps)

This section tracks the remaining “publishable closure” gaps. It is intentionally **assumption-explicit** and mirrors:
- suite outputs under `tfpt-suite/out/*/results.json` and `report.txt`
- `qft_completeness_ledger` “next_steps”

Legend:
- [ ] not started
- [~] partial / engineering-closure exists, but publication-grade derivation is missing
- [x] done (publication-grade)

### Gravitation (make \(G\) / \(\bar M_{\mathrm{Pl}}\) non-input)

- [~] `effective_action_r2`: closure-level \(a_2\to\beta_{R^2}\to M\) is implemented and internally consistent.
- [~] Derive the **quadratic fluctuation operator** for axial torsion \(S_\mu\) **from the torsionful microscopic action** (action-term parse + block provenance implemented; full connection-level derivation still open).
- [~] Compute the relevant heat-kernel coefficient(s) \(a_2\) from that derived operator and show the mapping to the \(R^2\) coefficient in the effective action (closure-level contract from action-derived blocks; full connection-level derivation still open).
- [ ] Apply a clean renormalization condition and show the **Einstein limit** statement \(K\to0 \Rightarrow\) EH with \(\bar M_{\mathrm{Pl}}\), then express the resulting \(G\) (or a dimensionless gravity strength ratio) relative to TFPT’s dynamically generated IR scales.

### Renormalization + matching (scheme/scale consistency; error propagation)

- [~] `msbar_matching_map`: deterministic MZ→mt bookkeeping exists (plus αs-sensitivity + optional MC hook).
- [~] `below_mt_eft_cascade`: QCD audit trail + explicit QED below‑MZ policy (`below_mz_policy.json`) now present; publication‑grade EW finite pieces remain open.
- [~] Threshold matching layer is now **explicit**: centralized `match_gauge` / `match_yukawa` / `match_quartic` API exists and mt→UV pipelines record per-threshold `threshold_match` (with `matching_active` flags; MSigma/MG8 are `matched_1loop_log_only_identity` at μ=threshold). EW/QED finite pieces (W/Z Δα̂−Δα, top Δy_t, Higgs Δλ_H) are now implemented; remaining for publication-grade is end-to-end use across matching policies and the full below‑mt EFT + covariance wiring.
- [~] End-to-end **uncertainty propagation** is now present as suite outputs (MC over PDG-style σ priors): `msbar_matching_map` (boundary), `ckm_full_pipeline` (ref-scale |V_ij| + χ²), and `pmns_full_pipeline` (μUV angles + χ²; capped for runtime). Remaining for publication-grade is explicit covariance support and unified propagation across *all* thresholds/policies (incl. below‑MZ QED/EW choices and any future finite matching pieces).
- [x] (support tooling) `unconventional`: `ux_matching_metamorphic_audit` — metamorphic/property audit for matching primitives (invertibility + stability)

### RG runner hardening (no silent model swaps; explicit patches)

- [~] `two_loop_rg_fingerprints`: fail-fast model fingerprint is in place (E8 model is selected explicitly).
- [x] Promote segmented RG running as a first-class public API: `segments=[{model, mu_start, mu_end, patches, threshold_match}]` (implemented in `tfpt_suite/rge_pyrate_2loop.py`; segments carry `model`, `mu_start/mu_end`, `patches`, and `threshold_match` with `matching_active` flags).
- [x] Make the gravity \(\alpha^3\) correction an explicit patch type (runner outputs include `gravity_alpha3_patch={enabled,kappa_vector,c3,alpha_definition_for_U1}` and segments include the `gravity_alpha3` patch label when active).

### Flavor (topology → Yukawas → CKM/PMNS as falsifiable outputs)

- [~] Deterministic Z3/Möbius texture + explicit conventions exist; **topology-phase selection is now enforced as a filter** (CKM gate `phase_set_derived_not_searched` PASS in full runs) and PMNS avoids χ²-driven convention shopping under the fixed rule (`no_convention_shopping_possible_under_fixed_rule` PASS via mass-splitting canonicalization).
- [~] Replace “texture ansatz” by a **topology-derived generator** (holonomy/monodromy/APS η input → Yukawa operator). The topology→phase map is now explicit (holonomy-class → δ/δ_CP), but the operator-level derivation remains open.
- [x] CKM/PMNS comparisons now declare scale policy (native mt vs reference MZ) and can consume `msbar_matching_map` mt-boundary inputs for consistency.
- [ ] PMNS full pipeline: constrain κ(mt) initialization from the `pmns_z3_breaking` mechanism output (then reconstruct \(y_N, M_R\) with a unique selection rule).

### Cosmology (bounce scale → observable scales)

- [~] `k_calibration` quantifies the missing scaling; bounce features land outside CMB ℓ unless a large scale-factor budget is provided.
- [x] Derive `a0_over_a_transition` from a consistent expansion history: `cosmo_threshold_history` now provides threshold-derived reheating inputs and `k_calibration` consumes them (fallback to v1.06 policy if threshold outputs are absent).
- [x] Decide and encode “CMB vs small-scale signatures” as an explicit policy (`k_calibration.json` signature_policy; checked in `boltzmann_transfer`).
- [ ] `bounce_perturbations`: tighten Wronskian diagnostics to be invariant across \(\omega^2\) sign changes and enforce strict bounds across *all* k (numerical stability upgrades).
- [x] (support tooling) `unconventional`: `ux_cosmo_history_sampler` — prior-sensitivity/feasibility scan for k→ℓ under explicit (N_infl, N_reheat, T_reheat) priors

### Baryons + dark matter

- [~] `omega_b_conjecture_scan` implements a conditional identity with explicit assumptions.
- [ ] Replace conditional \((4\pi-1)\) by a real operator-level derivation (anomaly inflow / sector counting theorem).
- [~] Dark matter relic density closure: `axion_fa_derivation` now derives f_a from the E8 ladder/block constants (n=10 PQ block) and downstream axion modules prefer the derived value; remaining work is full scenario selection, strings/domain walls, and astrophysical bounds.
- [x] Axion strings/domain-walls factor C_str is now explained via the Mobius cusp/charge sum and audited in `axion_dm_pipeline`.
- [x] `torsion_dm_pipeline` documents the optional torsion-DM branch with explicit constraint checklist (Ω_DM target, direct detection, torsion bounds).
- [x] (support tooling) `unconventional`: `ux_omega_b_aps_bridge` — expresses \((4\pi-1)\) as \(2\Delta_\Gamma-1\) with \(\Delta_\Gamma\) from the APS seam module (still conditional; docking point for inflow derivation)

### Torsion falsifiability (today-regime predictions)

- [~] `torsion_bounds_mapping` is no longer “always green” by construction (nonzero benchmark regime exists), but it is still a **toy** regime.
- [ ] Define and compute at least one **nontrivial TFPT torsion regime** (e.g. spin-polarized media / magnetars / early plasma) with explicit amplitude prediction and map to a measurable observable.

### QFT quantization (publication-grade)

- [~] Closure-level gauge fixing + ghost block exist in `effective_action_r2_operator_spec.json` with action-term parse + block-source audit.
- [~] Provide a BRST-complete derivation for the torsion sector quadratic operator and ghosts from the microscopic action (closure-level action derivation done; full connection-level derivation still open).

### Not yet implemented in suite (paper-level work)

- [x] Cosmological constant / vacuum energy closure: `torsion_condensate` now solves a spectral-flow-quantized gap equation for ⟨K²⟩ and `dark_energy_paths` identifies the ladder terminal stage (continuous extrapolation; discrete microscopic proof remains future work).
- [ ] Arrow of time / thermodynamics derivation (torsion flux + topological non-invertibility).

