# E8 Orbit Engine - Usage Guide

## Quick Start

### Strict ΔD=2 Mode (Theory-Conformant)

Generate a 27-step chain with purely ΔD=2 transitions:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --allowed-steps "2" \
  --beam-width 64 --top-k 3
```

**Expected Results:**
- 27 orbits from D=60 to D=8
- 0 large jumps (ΔD=4)
- CV(Δ³ln(D)) ≈ 1.75
- Report: `results/orbit_engine_report.html`
- Plots: `results/plots/`

### Exploratory Mode (ΔD ∈ {2,4})

Allow both ΔD=2 and ΔD=4 transitions for shorter chains:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --allowed-steps "2,4" \
  --beam-width 64 --top-k 3
```

**Expected Results:**
- ~14 orbits (shorter chain)
- ~13 large jumps (ΔD=4)
- Better CV(Δ³ln(D)) ≈ 1.47 (but less theoretical resolution)

### Baseline Chain (Hardcoded)

Use the predefined baseline chain:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --baseline
```

### Without HTML Report

Skip HTML report generation for faster execution:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --allowed-steps "2" \
  --no-html
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--allowed-steps` | "2,4" | Comma-separated allowed ΔD steps (e.g., "2" for strict ΔD=2) |
| `--beam-width` | 32 | Beam search width (higher = more thorough search) |
| `--top-k` | 10 | Number of top chains to save |
| `--label-threshold` | 1.0 | Maximum label distance for edges (1.0 = no filter) |
| `--max-height-diff` | 999 | Maximum height difference for edges (999 = no filter) |
| `--html/--no-html` | --html | Generate HTML report (default: on) |
| `--baseline` | False | Use hardcoded baseline chain instead of searching |

## Output Files

```
results/
├── orbit_engine_report.html     # Comprehensive HTML report
├── e8_chain_summary.json        # Chain and fit summary
├── data/
│   └── e8_chain_clean.csv       # Final chain data
├── chains/
│   ├── ranking.json             # Chain search ranking
│   └── top10/                   # Top chains from search
│       ├── rank_1.csv
│       ├── rank_2.csv
│       └── ...
├── plots/                       # Generated plots
│   ├── lnD_vs_n.png
│   ├── gamma_vs_n.png
│   ├── residuals.png
│   ├── d3_lnD.png
│   └── d3_phi.png
└── tables/
    └── e8_chain.csv             # Full data table
```

## Model Comparison

The engine compares four models for γ(n):

1. **Quadratic (n≥1)**: γ(n) ≈ a + bn + cn² (diagnostic only)
2. **Log-exact (anchored)**: γ(n) = λ·ln((60-2n)/(58-2n)) with λ fixed
3. **Log-exact (fitted)**: Same but with λ fitted
4. **Hyperbolic**: γ(n) ≈ A/(B-n) (best fit typically)

## Testing

Run regression tests:

```bash
cd e8_orbit_engine
python3 -c "
from tests.test_chain_strict_d2 import *
test_strict_d2_chain_has_27_steps_and_no_jumps()
test_normalization_anchor_gamma0()
test_a4_a1_has_corrected_dimension()
test_baseline_chain_properties()
print('All tests passed!')
"
```

## Key Results

### Strict ΔD=2 Chain
- **Length**: 27 orbits
- **D sequence**: 60, 58, 56, ..., 10, 8
- **CV(Δ³ln(D))**: ~1.75
- **Best model**: Hyperbolic (R² ≈ 1.000)

### Mixed ΔD={2,4} Chain
- **Length**: ~14 orbits
- **Large jumps**: ~13
- **CV(Δ³ln(D))**: ~1.47 (better but less resolution)

### Key Physics
- **Normalization**: γ₀ = 0.834 (anchored)
- **First step**: s₀ = ln(248) - ln(60) = 1.419
- **Scaling**: λ = 0.834/s₀ = 0.588
- **Topological constraint**: γ₂ ≈ γ₀/(8π²) fails (~87% deviation)
