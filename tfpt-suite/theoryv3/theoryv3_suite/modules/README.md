# theoryv3_suite.modules

Theoryv3 analysis modules focused on discrete building blocks and closure patterns.

## Modules

- `seed_invariants_audit`: pi -> core invariants consistency
- `defect_partition_g5_audit`: delta2 from g=5 and alpha(0) closure
- `alpha_backreaction_sensitivity_audit`: k sweep (alpha sensitivity vs backreaction exponent)
- `g5_origin_audit`: g=5 from SU(5) holonomy degeneracy (single origin)
- `dark_energy_exponential_audit`: exp(-alpha_inv/2) suppression and rho_L
- `dark_energy_norm_half_origin_audit`: n=1/2 from double-cover degree
- `flavor_pattern_audit`: lambda, delta_star, delta_cp, PMNS theta13 checks
- `pmns_tm1_audit`: TM1 sum rule for sin2 theta12
- `yukawa_exponent_index_audit`: rational indices for mass ratios
- `yukawa_index_mapping_audit`: map q_ij to charge-squared index sums
- `baryon_consistency_audit`: Omega_b identity, eta_b proxy, derived H0
- `axion_dm_audit`: axion frequency and relic fraction audit
- `g5_crosslink_audit`: g=5 consistency across sectors
- `constant_factory_audit`: grouped constant ledger (formula_math/formula_text + sensitivities + grammar checks + crosslinks + anchor-vs-TFPT views) with derivations, status tags, and reference comparisons (EW scale + lepton ladder marked pending as candidate formulas)

## Conventions

- ASCII-only reports (avoid Unicode in text outputs)
- Deterministic outputs and explicit assumptions
- Plots are best-effort and must not fail module execution
