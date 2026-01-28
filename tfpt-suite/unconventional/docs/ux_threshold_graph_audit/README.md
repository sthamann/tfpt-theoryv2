# `ux_threshold_graph_audit`

## Purpose

This module addresses `tfpt-suite/unconventional/tasks.md` (Section B) deliverable:

- **a declarative threshold graph**: which theory/model is active in which \(\mu\)-interval, and what “matching actions” happen at boundaries.

The main suite’s 2-loop runner already implements this declaratively via:

- `tfpt_suite/rge_pyrate_2loop.py`
  - `threshold_rules` (what happens at each threshold)
  - `segments` (piecewise integration with explicit segment metadata)

This unconventional module extracts that into a standalone, inspectable artifact and compares:

- **matching disabled** (segments labeled “continuous_by_assumption”)
- **matching enabled** (explicit `match_gauge/match_yukawa/match_quartic` outcomes recorded; 1-loop is identity at μ=threshold)

## Inputs

- Threshold table: `tfpt-suite/tfpt_suite/data/rge_thresholds_v25.json`
- SM boundary at \(m_t\): derived deterministically from `tfpt-suite/tfpt_suite/data/sm_inputs_mz.json` via `tfpt_suite/pyrate_boundary_runner.py` (PyR@TE-driven MZ→mt, RG-authority compliant)
- PyR@TE beta sources from `tfpt-suite/tfpt_suite/data/pyrate_pythonoutputs.json`

Additionally (Gate-2 matching audit):

- below-MZ heavy-quark thresholds `mc_GeV, mb_GeV` from `sm_inputs_mz.json`, used to verify finite 2-loop \(\alpha_s\) decoupling when `apply_alpha3_matching=True`.

## What it computes

- **Nodes**: thresholds (sorted), e.g. MSigma, MG8, MNR1..3
- **Edges**: segments \([\mu_\mathrm{start},\mu_\mathrm{end}]\) with:
  - active model (`sm_tfpt_2loop_v25` or `e8_sigma_yN_2loop`)
  - active patches (e.g. `delta_b3_g8`, `gravity_alpha3` if enabled)
  - threshold transition at segment start (if any)
  - matching record (enabled/disabled)
- **Threshold status** (matching ON): `finite_pieces_implemented` vs `identity_matching`, based on recorded `match_gauge`/`match_yukawa`/`match_quartic` outcomes.

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_threshold_graph_audit
```

Outputs:

- `tfpt-suite/out/unconventional/ux_threshold_graph_audit/results.json`
- `tfpt-suite/out/unconventional/ux_threshold_graph_audit/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_threshold_graph_audit/threshold_graph.png`

## Interpretation

- If **matching is disabled**, the runner will (correctly) mark thresholds as “continuous_by_assumption”.
- If **matching is enabled**, the runner records explicit matching events even if the numerical step is zero at 1-loop.
- The report now lists each threshold with a **finite-piece status tag** (`finite_pieces_implemented` vs `identity_matching`) to make publication‑grade gaps explicit.
- The module also reports a **below-MZ finite matching diagnostic**: \(\alpha_s(m_c)\) with matching ON vs OFF (2-loop decoupling), to confirm the “finite pieces” milestone is actually wired into the EFT cascade.

This is the right scaffolding for publication-grade matching:

- add finite matching pieces where applicable
- keep the graph/segment metadata stable and reviewer-proof

