# Unconventional TFPT Suite (`tfpt-suite/unconventional/`)

This folder hosts an **auxiliary “unconventional” suite** that *docks to* the main `tfpt-suite/` infrastructure, but is intentionally **kept separate** from the publication-facing verification modules in `tfpt_suite/modules/`.

The purpose is to implement **search / audit / conjecture-engine tooling** (genetic algorithms, brute enumeration, counterexample search, metamorphic/property tests, etc.) that can accelerate closing the explicit ToE gaps tracked in:

- `tfpt-suite/progress.md` (Publication-grade ToE closure TODO)
- `tfpt-suite/unconventional/tasks.md` (idea catalog + guardrails)

## What this suite is (and is not)

- **It is**: a deterministic, reproducible way to generate **candidates**, **stress tests**, and **robustness statements** for open derivation gaps (e.g. gauge-fixing choices, scale-history feasibility, matching consistency).
- **It is not**: a substitute for formal operator/QFT derivations. Anything “physics-final” still needs a paper-grade derivation and/or falsifiable predictions in the main suite.

## Guardrails (anti “fit / p-hack”)

These are enforced by design and documented per module:

- **No direct experiment-fitting as a fitness function**. Search objectives are structural (minimality, invariances, consistency, stability). If an external comparison is used, it is clearly marked as *diagnostic* and (where applicable) placed behind holdout-style evaluation.
- **Assumptions are explicit**. Every module emits a structured list of assumptions and the exact conventions used.
- **Deterministic outputs**. All randomness is seeded via `SuiteConfig.seed`.

## How it docks to the main suite

- Reuses `tfpt_suite.config.SuiteConfig` and the `TfptModule` interface (`tfpt_suite/module_base.py`).
- Writes artifacts under a separate output root (default: `tfpt-suite/out/unconventional/`).
- Has its own runner and registry:
  - `run_unconventional_suite.py`
  - `tfpt_unconventional/modules/registry.py`

## How to run

From repo root:

```bash
# List unconventional modules
python3 tfpt-suite/unconventional/run_unconventional_suite.py list-modules

# Run all unconventional modules
python3 tfpt-suite/unconventional/run_unconventional_suite.py run-all

# Run a subset
python3 tfpt-suite/unconventional/run_unconventional_suite.py run --modules ux_matching_metamorphic_audit,ux_omega_b_aps_bridge
```

Outputs land in:

- `tfpt-suite/out/unconventional/<module_id>/` (`results.json`, `report.txt`, `meta.json`, optional plots)

## Module documentation

Each unconventional module has an accompanying design/usage note under:

- `tfpt-suite/unconventional/docs/<module_id>/README.md`

