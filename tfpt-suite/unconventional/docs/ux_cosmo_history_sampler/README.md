# `ux_cosmo_history_sampler`

## Purpose

This module explores the **k→ℓ bridge feasibility** by sampling explicit expansion-history parameters and computing the resulting bounce multipole location.

It is intentionally **not a “fit”**. The search space is defined by explicit priors (ranges), and the output is a **robustness statement**:

- under these declared priors, does the bounce feature land in the CMB ℓ window?
- if not, what ℓ range is robustly implied (small-scale prediction)?

## Why this matters for “ToE closure”

`tfpt-suite/progress.md` flags the k→ℓ bridge as a key open closure item:

- `k_calibration` quantifies the missing overall scaling
- what is still missing is a **publication-grade expansion history policy/module** that deterministically yields \(a_0/a_{\mathrm{transition}}\)
- the project must decide whether TFPT predicts **CMB-visible** bounce features or **small-scale-only** features

This module is the “unconventional” complement to `k_calibration`:

- `k_calibration` is deterministic for a fixed policy input
- this module probes **policy sensitivity** and feasibility under explicit priors

## What the module computes

Inputs (aligned with the main suite):

- bounce scale \(k_{\mathrm{bounce}}^{s,t}\):
  - prefers `tfpt-suite/out/bounce_perturbations/results.json` if present
  - otherwise uses the fallback values in `tfpt_suite/data/k_calibration.json`
- cosmology snapshot (flat ΛCDM) from `tfpt_suite/data/k_calibration.json`
- scale-factor mapping in `tfpt_suite/cosmo_scale_map.py`:
  - entropy mapping after reheating
  - optional reheating expansion `exp(N_reheat)`
  - inflation expansion `exp(N_inflation_from_transition)`

For each scenario (set of priors), it samples:

- `N_inflation_from_transition`
- `N_reheat`
- `T_reheat_GeV` (log-uniform prior)

Additionally, the module includes a **policy-driven scenario** (`policy_v106_threshold_driven`) that reduces free knobs:

- \(N_{\mathrm{pivot}}\) is derived from `global_reference_minimal.json` via the Starobinsky relation \(N\simeq 2/(1-n_s)\)
- \(N_{\mathrm{reh}}\) is derived from \(T_{\mathrm{reh}}\) via the v1.06 reheating \(\rho\)-scaling model (`cosmo_reheating_policy_v106.json`), including the declared \(\Delta N\) floor filter

and computes:

- \(a_0/a_{\mathrm{transition}}\)
- \(\ell_{\mathrm{bounce}}^{s,t}\)
- fraction of samples with \(\ell\) in a diagnostic “CMB window” (default: \([2,2500]\))
- best candidates closest to target ℓ values (in log-space)

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_cosmo_history_sampler
```

Outputs:

- `tfpt-suite/out/unconventional/ux_cosmo_history_sampler/results.json`
- `tfpt-suite/out/unconventional/ux_cosmo_history_sampler/report.txt`
- Plots (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_cosmo_history_sampler/cosmo_history_sampler_log10ell.png`
  - `tfpt-suite/out/unconventional/ux_cosmo_history_sampler/cosmo_history_sampler_cmb_fraction.png`

## Interpreting results

- If **plausible priors** yield **nonzero CMB hits**:
  - then bounce visibility is not excluded
  - the next step is to **encode the expansion history policy explicitly** (reheating model, g\*(T), thresholds, consistency with BBN) so the k→ℓ bridge becomes publication-grade
- If even **extended priors** yield **zero CMB hits**:
  - then the theory is making a strong statement: bounce features are likely **small-scale-only**
  - downstream work should target the correct observables (e.g. μ-distortions, small-scale structure, etc.) instead of CMB ℓ

## Next steps (towards publication-grade)

Implement a deterministic “cosmology history” module (main suite, not unconventional) that:

- computes \(a_0/a_{\mathrm{transition}}\) from an explicit chain:
  - reheating temperature model + entropy conservation
  - g\*(T) history and threshold steps
  - BBN consistency constraints (minimum \(T_{\mathrm{reh}}\))
- declares which transition is used (bounce end vs inflation start vs pivot exit)
- outputs a reproducible uncertainty band if priors are still needed (explicit MC, deterministic seed)

