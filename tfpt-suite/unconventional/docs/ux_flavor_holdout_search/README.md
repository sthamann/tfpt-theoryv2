# `ux_flavor_holdout_search`

## Purpose

This module implements the **“holdout principle”** from `tfpt-suite/unconventional/tasks.md` (Section C):

> search only over **discrete** convention choices, evaluate with a **holdout split**, and report both.

It is a safeguard against *silent convention shopping* in the flavor sector:

- if you allow yourself to try many discrete phase conventions, you can accidentally overfit the CKM table “for free”.
- the holdout protocol forces you to keep some observables unseen during ranking.

## What it searches

A finite discrete grid of convention choices:

- `delta_source ∈ {delta_star, tau_mu}`
- `s13_mode ∈ {A_lam3_over_3, A_lam3_times_1_minus_delta}`
- `delta_mode ∈ {pi_times_delta, pi_times_1_minus_delta, 2pi_times_delta, koide_pi_over_12}`

No continuous parameters are tuned.

## Holdout split

- **Fit keys (7)**: `Vud, Vus, Vcd, Vcs, Vcb, Vts, Vtb`
- **Holdout keys (2)**: `Vub, Vtd`

Rationale: the smallest entries are most sensitive and easiest to “convention-fix”.

## PMNS holdout (added)

The module also performs a **PMNS holdout** to avoid silently “breaking” \(\delta_{CP}\) via convention shopping:

- It runs `pmns_full_pipeline` once to obtain the **discrete PMNS convention candidates** (the 6 column permutations).
- It evaluates a holdout split against `tfpt_suite/data/pmns_reference.json`:
  - **Fit keys (3)**: `sin2_theta12, sin2_theta13, sin2_theta23`
  - **Holdout key (1)**: `delta_cp_deg`

This keeps the search **discrete-only** (no continuous fitting) while surfacing whether a “good” fit is robust when \(\delta_{CP}\) is held out.

## Scale policy

Reference is labeled at \(M_Z\) in:

- `tfpt-suite/tfpt_suite/data/ckm_reference.json`

The module therefore runs:

- mt → \(M_Z\) (short down-run) using `tfpt_suite/rge_pyrate_2loop.py`

This keeps the comparison at least *scale-aware* (still diagnostic; the reference is a low-energy fit snapshot, not a running MSbar parameter).

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_flavor_holdout_search
```

Outputs:

- `tfpt-suite/out/unconventional/ux_flavor_holdout_search/results.json`
- `tfpt-suite/out/unconventional/ux_flavor_holdout_search/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_flavor_holdout_search/flavor_holdout.png`

## Interpreting results

You should look for candidates that:

- have low **χ²_fit**,
- do not catastrophically fail **χ²_holdout**,
- and have low “complexity” (explicitly reported as a simple integer proxy).

This does *not* replace the real ToE gap:

- deriving the topology/holonomy → Yukawa operator map (and phase branches) at operator level.

But it provides a robust audit trail while that derivation is pending.

