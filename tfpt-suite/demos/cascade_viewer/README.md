# TFPT Cascade Viewer

Interactive visualization of the E8→SM symmetry breaking cascade from the TFPT (Topological Field-Point Theory) suite.

## Important Design Principle

**This demo does NOT compute anything.** It exclusively reads artifacts produced by the TFPT analysis suite. All values displayed come directly from:

- `tfpt_suite/data/*.json` (input data files)
- `out/*/results.json` (module output files)

If the analysis files change, the visualization automatically reflects those changes. If an analysis file is missing, the demo crashes hard — this is intentional to ensure the visualization always represents the actual state of the theory.

## Quick Start

```bash
cd tfpt-suite/demos/cascade_viewer

# Install dependencies
npm install

# Run in development mode
npm run dev

# Or build and run
npm run build
npm start
```

Then open http://localhost:3333 in your browser.

## Architecture

```
cascade_viewer/
├── config/
│   └── sources.yaml      # SINGLE SOURCE OF TRUTH for data paths
├── loader/               # Data loading (reads only, no computation)
│   ├── cascade.ts        # RG thresholds + two-loop fingerprints
│   ├── alpha.ts          # α precision audit
│   ├── spectrum.ts       # APS η-gluing + R² blocks + bounce
│   ├── torsion.ts        # Torsion bounds mapping
│   └── global.ts         # Global χ² consistency
├── state/
│   └── tfpt_store.ts     # Centralized state (all views read from here)
├── views/                # SVG rendering (no data access)
│   ├── LadderView.ts     # Cascade ladder visualization
│   ├── AlphaMonitor.ts   # CFE fixpoint + k-sensitivity
│   ├── SpectrumHeatmap.ts # Spectral flow + operator blocks
│   ├── TorsionBounds.ts  # SME bounds + prediction
│   └── GlobalRadar.ts    # χ² radar + scorecard
└── index.ts              # Server (ONLY file that does I/O)
```

## Views and Their Data Sources

### 1. Cascade Ladder

**What it shows:** Vertical ladder of E8 cascade stages with clickable nodes.

**Data sources:**
- `tfpt_suite/data/rge_thresholds_v25.json` — Threshold scales (MSigma, MG8, MNR1-3, MPhi)
- `out/two_loop_rg_fingerprints/results.json` — α₃ running, crossings at φ₀ and c₃

**Paper reference:** Section 4 (E8 Cascade), Table 2

### 2. Alpha Monitor

**What it shows:** CFE self-consistency iteration with k-sensitivity curve.

**Data source:**
- `out/alpha_precision_audit/results.json` — k-grid, α⁻¹ values, ppm deviation, convergence

**Paper reference:** Section 3.2, Eq. (3.7)

**Formula displayed:**
```
CFE: α³ - 2c₃³α² - 8b₁c₃⁶ ln(1/φ) = 0
Backreaction: φ(α) = φ_tree + δ_top exp(-kα)
```

### 3. Spectrum Heatmap

**What it shows:** APS spectral flow table, R² operator block contributions, transfer functions.

**Data sources:**
- `out/aps_eta_gluing/results.json` — Spectral flow SF(m), η invariants
- `out/effective_action_r2/results.json` — Operator blocks, β_R² contributions, M/M_Pl
- `out/bounce_perturbations/results.json` — T_scalar(k), T_tensor(k)

**Paper reference:** Appendix K (APS gluing), Appendix L (bounce perturbations)

### 4. Torsion Bounds

**What it shows:** SME-style torsion bounds vs TFPT prediction.

**Data sources:**
- `out/torsion_bounds_mapping/results.json` — Inferred S_μ bounds, prediction
- `tfpt_suite/data/torsion_regimes.json` — Regime definitions

**Paper reference:** Section 6, Appendix M

**Key relationship:** `b_μ = (3/4) · S_μ` (minimal coupling mapping)

### 5. Global Radar

**What it shows:** χ² scorecard with z-score radar and contribution bars.

**Data source:**
- `out/global_consistency_test/results.json` — Per-observable χ², z-scores, p-value

**Paper reference:** Table 1, Section 7

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Main HTML page with all views |
| `/api/store` | GET | Full store data as JSON |
| `/api/summary` | GET | Store summary (loaded counts) |
| `/api/set-level/:n` | POST | Set active cascade level for view sync |

## Acceptance Criteria

This demo is considered correct when:

1. **Re-running the TFPT suite automatically updates all visuals** — No manual parameter adjustment needed.

2. **Removing an analysis file causes a hard crash** — The demo must not silently fall back to defaults.

3. **Every displayed value traces back to a specific file** — Hover tooltips and source labels provide full traceability.

4. **View synchronization works** — Clicking a cascade stage updates ACTIVE_LEVEL, which all views respect.

## Traceability Table

| View Component | Source File | Analysis Script |
|----------------|-------------|-----------------|
| Cascade stages | `rge_thresholds_v25.json` | — (data file) |
| α₃ crossings | `two_loop_rg_fingerprints/results.json` | `modules/two_loop_rg_fingerprints.py` |
| k-sensitivity | `alpha_precision_audit/results.json` | `modules/alpha_precision_audit.py` |
| Spectral flow | `aps_eta_gluing/results.json` | `modules/aps_eta_gluing.py` |
| R² blocks | `effective_action_r2/results.json` | `modules/effective_action_r2.py` |
| Transfer T(k) | `bounce_perturbations/results.json` | `modules/bounce_perturbations.py` |
| S_μ bounds | `torsion_bounds_mapping/results.json` | `modules/torsion_bounds_mapping.py` |
| χ² terms | `global_consistency_test/results.json` | `modules/global_consistency_test.py` |

## Philosophy

> This is not a demo folder. This is a visible proof that TFPT is a running machine.

The visualization layer is intentionally thin:
- **Loaders** only parse and restructure — no math
- **Store** only holds state — no computation
- **Views** only render — no data access
- **Server** only serves — no business logic

This ensures that what you see is exactly what the theory produces.
