# `ux_omega_b_aps_bridge`

## Purpose

The main suite currently contains an explicitly **conditional** Ω\_b identity module:

- `tfpt_suite/modules/omega_b_conjecture_scan.py`

There, the coefficient \(4\pi-1\) is justified via a **sector-counting assumptions block**.

This unconventional module provides a more “docked” candidate bridge:

> express \(4\pi-1\) as a function of an already-computed *topological seam quantity* in the suite.

Concretely it uses the APS seam toy model:

- `tfpt_suite/modules/aps_eta_gluing.py`

and its seam term:

\[
  \Delta_\Gamma := 2\pi \cdot \mathrm{SF}(U_\Gamma)
\]

For the minimal nontrivial class \(m=1\), the suite computes \(\mathrm{SF}=1\) and hence \(\Delta_\Gamma=2\pi\).

Then the coefficient becomes:

\[
  K_{\mathrm{candidate}} := 2\Delta_\Gamma - 1 = 4\pi - 1.
\]

## Why this matters for “ToE closure”

In `tfpt-suite/progress.md`, Ω\_b is flagged as:

- implemented as a conditional identity
- still missing an **operator/anomaly-level derivation** (inflow / index / η-gluing logic)

This module does **not** claim to finish that derivation. What it does provide is a concrete docking point:

- Ω\_b coefficient ↔ APS seam term \(\Delta_\Gamma\)

That makes it much clearer what a future anomaly/inflow calculation would have to reproduce.

## What the module computes

- Recomputes the APS seam spectral flow for \(m=1\) (periodic spin) using the same toy model as `aps_eta_gluing`.
- Forms \(K_{\mathrm{candidate}} = 2\Delta_\Gamma - 1\).
- Computes:
  - \( \Omega_b^{\mathrm{pred}} = K_{\mathrm{candidate}} \beta_{\mathrm{rad}}\)
  - Planck-derived reference Ω\_b and the implied coefficient \(K := \Omega_b/\beta_{\mathrm{rad}}\)
  - a diagnostic z-score (agreement indicator, not a proof)
- Emits an explicit assumptions list; the key bridge step remains **conditional**.

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_omega_b_aps_bridge
```

Outputs:

- `tfpt-suite/out/unconventional/ux_omega_b_aps_bridge/results.json`
- `tfpt-suite/out/unconventional/ux_omega_b_aps_bridge/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_omega_b_aps_bridge/omega_b_aps_bridge.png`

## Next steps (upgrade path)

To turn this into a publication-grade Ω\_b derivation:

- Replace the bridge postulate “\(K=2\Delta_\Gamma-1\)” by an explicit **operator/anomaly/inflow computation**.
- Keep the *same* final numeric regression target \(K\approx 11.65\), but ensure the chain is:
  - microscopic action → operator → η/Index/inflow → seam contribution → Ω\_b coefficient
- Maintain the suite’s policy: all assumptions must remain explicit, and diagnostics must not be promoted to derivations silently.

