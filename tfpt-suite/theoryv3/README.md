# theoryv3 (tfpt-suite)

Focused analysis branch that formalizes the "discrete building blocks" thesis:
- pi as seed -> core invariants
- defect partition (g=5) -> delta2 -> alpha(0)
- g=5 origin from SU(5) holonomy degeneracy
- alpha backreaction sensitivity (k sweep)
- exponential suppression -> dark energy scale
- n=1/2 normalization from double cover
- Mobius + Z3 flavor formulas -> CKM/PMNS anchors
- Yukawa exponent indices -> charge-squared index sums
- TM1 sum rule check for sin2 theta12
- baryon sector identities -> Omega_b and eta_b
- axion DM target (frequency, relic fraction)
- constant factory -> grouped constant tables + ledger with formula_math/formula_text, sensitivities, grammar checks, crosslink signatures, and anchor-vs-TFPT views

Outputs are deterministic, report-ready, and meant for a single consolidated PDF
with an English introduction, a test map, per-test A-E sections (goal, inputs,
expected result, method, results/JSON), and plot verification.

## Quickstart

```bash
# Requires tfpt-suite/out_physics/ to exist

# Run all theoryv3 modules (engineering mode)
python3 tfpt-suite/theoryv3/run_theoryv3.py run-all

# Build the theoryv3 PDF report
python3 tfpt-suite/theoryv3/run_theoryv3.py build-report

# Validate artifacts against schema
python3 tfpt-suite/theoryv3/run_theoryv3.py validate-artifacts

# Build an auto-generated manifest (JSON + TeX table) for paper integration
python3 tfpt-suite/theoryv3/theoryv3_manifest.py

# Run theoryv3 tests (fast, temporary outputs)
python3 tfpt-suite/theoryv3/run_theoryv3.py run-tests
```

## Directory structure

```
tfpt-suite/theoryv3/
├── run_theoryv3.py            # CLI runner for theoryv3 analyses
├── report_builder.py          # PDF report builder (theoryv3)
├── theoryv3_suite/
│   ├── modules/               # theoryv3 analysis modules
│   ├── utils.py               # shared helpers (IO, parsing, math)
│   └── README.md
├── tests/                     # unit tests for theoryv3
└── out/                       # generated outputs (results.json, report.txt, plots)
```

## Outputs

- Module outputs: `tfpt-suite/theoryv3/out/<module_id>/`
- PDF report: `tfpt-suite/theoryv3/theoryv3-analysis.pdf`
- JSON artifacts: schema v1 (validated via `validate-artifacts`)
- Constant factory ledger: `out/constant_factory_audit/constant_factory_ledger.{json,tex}`
- Symbol registry: `tfpt-suite/theoryv3/theoryv3_suite/symbols.yaml`
- Paper integration artifacts (auto-generated):  
  `tfpt-suite/theoryv3/out/theoryv3_manifest_{summary,table}.tex` + `theoryv3_manifest.json`  
  (summary defines totals including PASS/WARN/FAIL check counts via `\\TFPTTheoryVThreeChecks{Pass,Warn,Fail}Total`)

## Data sources

Theoryv3 modules prefer current numbers from:
- `tfpt-suite/out_physics/<module_id>/results.json` (if present)
- otherwise fall back to the deterministic formulas in `tfpt_suite`
- reference datasets from `tfpt-suite/tfpt_suite/data/references.json`

## Notes

- This branch is analysis-focused: it does not modify the main suite.
- Plots are optional and gated by `--no-plot`.
- `constant_factory_audit` treats the EW-scale + lepton ladder as pending candidates until a discrete coefficient justification is available.
