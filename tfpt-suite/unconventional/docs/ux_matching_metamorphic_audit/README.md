# `ux_matching_metamorphic_audit`

## Purpose

This module implements **metamorphic / property-style audits** for the low-level matching primitives in:

- `tfpt_suite/matching.py`

It is an **engineering hardening layer**: it does *not* add physics content, but it reduces the risk of silent regressions and makes the “scheme/threshold bookkeeping” reviewer-proof.

## Why this matters for “ToE closure”

In `tfpt-suite/progress.md`, a key publication-grade gap is **renormalization + threshold matching**:

- finite matching pieces at thresholds (beyond “identity at μ=M”)
- uncertainty propagation across thresholds/policies
- stable, explicit APIs (`match_gauge`, `match_yukawa`, `match_quartic`)

Before adding more physics, we want strong evidence that the matching plumbing behaves as intended.

This module addresses that by verifying **metamorphic invariants** such as:

- **Invertibility (approx/exact)**:
  - αs quark threshold: `down(up(alpha3)) ≈ alpha3` (series-truncation-level)
  - gauge matching with finite δα: `down(up(g_i, δα)) = g_i` (should be exact up to float rounding)
- **Stability**:
  - matching never produces negative α for the sampled δα range

## What the module computes

- A deterministic set of random samples (seeded by `SuiteConfig.seed`)
- For each property, it computes:
  - max absolute error
  - max relative error
  - the worst-case sample (inputs + outputs)
- The module returns explicit `checks[]` that fail if errors exceed conservative tolerances.

## How to run

From repo root:

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_matching_metamorphic_audit
```

Outputs:

- `tfpt-suite/out/unconventional/ux_matching_metamorphic_audit/results.json`
- `tfpt-suite/out/unconventional/ux_matching_metamorphic_audit/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_matching_metamorphic_audit/matching_metamorphic_errors.png`

## How to interpret results

- If **`match_gauge_finite_delta_up_down_identity` fails**:
  - treat this as a **blocking bug** for any serious matching work
  - fix the sign conventions / conversion logic before adding finite pieces
- If **`alpha3_threshold_up_down_near_identity` fails**:
  - this indicates the current αs decoupling approximation has drifted beyond its intended truncation accuracy
  - tighten the implementation, or adjust tolerances only with a clearly justified loop-order/truncation argument

## Next steps (when closing the matching gap)

- Add finite matching pieces at the thresholds where they are known and relevant (EW/QCD, heavy new fields).
- Extend the metamorphic audit suite with:
  - round-trip invariance for whole RG segments (run up → match → run down)
  - “metamorphic reparameterization” tests (unit scaling, convention conversions)
  - covariance-aware uncertainty propagation checks once covariance is introduced

