# `ux_gravity_gaugefix_ga`

## Purpose

This module implements the **“GA search scaffold”** described in `tfpt-suite/unconventional/tasks.md` (Section A).

It searches (in a **toy operator family**) for gauge-fixing choices that make the quadratic operator **Laplace-type** (minimal), i.e. remove nonminimal \(\nabla_\mu\nabla_\nu\) structure.

## What it is (and is not)

- **It is**: a deterministic genetic-algorithm harness over a small discrete+continuous genome that minimizes a *structural nonminimality measure*.
- **It is not**: the publication-grade BRST derivation of the torsion-sector quadratic operator from the microscopic action. That remains the real ToE gap.

## Toy model used

We model a vector-type quadratic operator as:

\[
  \Delta_{\mu\nu} = -g_{\mu\nu}\Box + c_{\mathrm{nonmin}}\,\nabla_\mu\nabla_\nu + E\,g_{\mu\nu}.
\]

The module uses a toy parameterization:

\[
  c_{\mathrm{nonmin}} = 1 - \gamma/\xi,
\]

with:

- \(\xi>0\): gauge parameter
- \(\gamma\): discrete normalization choice of the gauge-fixing functional

Laplace-type is achieved when \(c_{\mathrm{nonmin}}=0\), i.e. \(\xi=\gamma\).

This is the minimal place where a GA is useful: as soon as you have a real operator family with many discrete choices, you want a deterministic search layer that optimizes **structure** (not experimental fit).

## Outputs

The module reports:

- best genome (family, \(\gamma\), \(\xi\))
- best nonminimal coefficient \(c_{\mathrm{nonmin}}\)
- a small top-10 snapshot of the final GA population
- “docking context” info about `effective_action_r2_operator_spec.json` (for narrative continuity)

## How to run

```bash
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_gravity_gaugefix_ga
```

Outputs:

- `tfpt-suite/out/unconventional/ux_gravity_gaugefix_ga/results.json`
- `tfpt-suite/out/unconventional/ux_gravity_gaugefix_ga/report.txt`
- Plot (when plotting is enabled):
  - `tfpt-suite/out/unconventional/ux_gravity_gaugefix_ga/gravity_gaugefix_ga.png`

## Upgrade path (towards publication-grade)

Once a real torsion-sector quadratic operator family is parameterized, extend the genome to include:

- explicit gauge-fixing functional basis choices
- field redefinitions that reduce derivative mixing
- discrete ghost-sector choices (BRST-grade)

Then replace the toy nonminimality measure by an operator-level criterion:

- penalty for off-diagonal \(\nabla_\mu\nabla_\nu\) terms
- penalty for non-selfadjoint pieces
- bonus for a clean Laplace-type block form \(-\nabla^2 + E\)

and add metamorphic checks:

- gauge-parameter variation should not change physical \(a_2\) invariants (once the full derivation exists)

