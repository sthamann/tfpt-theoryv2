# theoryv3_suite

Internal package for the theoryv3 analysis branch.

## Contents

- `modules/`: analysis modules (each writes results.json, report.txt, plots)
- `runtime.py`: centralized numeric context (mp.dps + deterministic RNG seed)
- `utils.py`: shared helpers for loading prior outputs and computing metrics
- `symbols.yaml`: symbol registry for constant factory formula grammar checks

## Output schema

Theoryv3 modules follow the same output schema as the main suite:

```
out/<module_id>/
  meta.json
  results.json
  report.txt
  *.png
```

## Notes

- Modules use deterministic formulas and prefer existing `tfpt-suite/out_physics` results when available.
- No changes are made to the main suite; theoryv3 is analysis-only.
- constant_factory_audit emits grouped constant tables plus a ledger (formula_math/formula_text, sensitivities, grammar checks, crosslinks, anchor-vs-TFPT views) for PDF rendering.
