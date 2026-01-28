# progress (theoryv3)

Status tracking for theoryv3 analyses and tests.

## Status legend

- [ ] planned / not started
- [~] in progress
- [x] implemented (with automated checks/tests)

## Core

- [x] theoryv3 suite skeleton (runner, registry, utilities)
- [x] PDF report builder (theoryv3-analysis.pdf)
- [x] unit tests (theoryv3/tests)

## Analyses

- [x] seed invariants audit (pi -> c3, varphi0, beta_rad)
- [x] defect partition g=5 audit (delta2 -> alpha_inv_0)
- [x] alpha backreaction sensitivity audit (k sweep)
- [x] g=5 origin audit (single-origin SU(5) holonomy degeneracy)
- [x] dark energy exponential audit (phi_star, rho_L)
- [x] dark energy norm origin audit (n=1/2 from double cover)
- [x] flavor pattern audit (lambda, delta_star, delta_cp, PMNS theta13)
- [x] PMNS TM1 audit (sin2 theta12 from TM1 sum rule)
- [x] yukawa exponent index audit (q_ij rationalization)
- [x] yukawa index mapping audit (q_ij to charge-squared sums)
- [x] baryon consistency audit (Omega_b, eta_b, derived H0)
- [x] axion DM audit (nu, Omega_a h^2)
- [x] g=5 crosslink audit (delta2, unification patch, gamma0)
- [x] constant factory audit (ledger + sensitivities + grammar + crosslinks)

## Infrastructure

- [x] references.json ledger (centralized datasets + versions)
- [x] results.json schema v1 + validate-artifacts gate
- [x] constant ledger schema v1 + validator hook
- [x] symbol registry for constant formulas + report appendix split
