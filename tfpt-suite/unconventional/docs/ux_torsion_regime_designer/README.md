# `ux_torsion_regime_designer`

## Purpose

`tfpt-suite/progress.md` flags torsion falsifiability as a ToE closure gap:

- the suite has a vetted bounds ingestion (`torsion_bounds_mapping`)
- but the “today torsion” prediction is still a **toy regime**
- what’s missing is a **nontrivial, computable present-day regime** (spin-polarized matter / magnetars / plasma)

This module is an *unconventional design tool* to propose such regimes in a way that is:

- explicit about assumptions
- directly comparable to vetted bounds
- compatible with the existing `torsion_bounds_mapping` regime evaluation

## Toy model used (explicit, not derived)

The module uses a minimal spin-medium proxy:

\[
|S_\mu| \sim c_{\rm spin}\;\frac{\rho_{\rm spin}}{M_{\rm eff}^2},
\qquad
\rho_{\rm spin} \equiv (\text{polarization})\cdot(\text{spin\_per\_particle})\cdot n
\]

where:

- \(n\) is a number density (converted to GeV\(^3\))
- \(M_{\rm eff}\) is either:
  - `Mpl_reduced` (very conservative), or
  - `TFPT_M` (uses the TFPT Starobinsky scale \(M\))
- \(c_{\rm spin}\) is a dimensionless factor (not derived here)

The key design output is:

\[
c_{\rm spin,max} := \frac{|S_\mu|_{\rm bound}}{|S_\mu|_{\rm pred}(c_{\rm spin}=1)},
\]

so you can see **how much suppression** is required to stay within bounds for a given physical scenario.

## Outputs

- tightest vetted component-wise bound (from `torsion_bounds_vetted.json`)
- a few “anchor scenarios” (lab matter, nuclear matter)
- random scenario exploration to find interesting candidates
- a simple **observable proxy** per scenario: \(\Delta\nu \approx 2|b|\cdot(\mathrm{GeV}\to\mathrm{Hz})\) with \(|b|\approx(3/4)|S|\) (design-phase; not a full experiment model)
- **ready-to-copy JSON regime proposals** compatible with `tfpt_suite/data/torsion_regimes.json`

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_torsion_regime_designer
```

Outputs:

- `tfpt-suite/out/unconventional/ux_torsion_regime_designer/results.json`
- `tfpt-suite/out/unconventional/ux_torsion_regime_designer/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_torsion_regime_designer/torsion_regime_designer.png`

## Next steps (towards publication-grade)

1. Decide which physical environment is the intended falsification target (magnetar core, lab spin medium, early plasma…).
2. Replace the toy spin-medium mapping by a TFPT-derived coupling (operator-level).
3. Add an explicit observable mapping (spin precession, birefringence in strong fields, timing signatures, etc.).
4. Integrate the chosen regime into `torsion_regimes.json` and make it part of the main suite’s falsification dashboard.

