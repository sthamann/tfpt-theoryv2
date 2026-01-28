# Topological Fixed Point Theory (TFPT)

**A parameter-free derivation of the fine-structure constant and other fundamental quantities from topology, geometry, and quantum consistency.**

**Authors:** Stefan Hamann & Alessandro Rizzo

**Website:** [www.fixpoint-theory.com](https://www.fixpoint-theory.com)

**Contact:** sh@sh-future.de

**Version:** 2.5 — January 27, 2026

---

## Table of Contents

1. [Overview](#overview)
2. [The Core Derivation](#the-core-derivation)
3. [The 8 Axioms](#the-8-axioms)
4. [Key Results](#key-results)
5. [Repository Structure](#repository-structure)
6. [Getting Started](#getting-started)
7. [Components](#components)
8. [Predictions & Falsification](#predictions--falsification)
9. [Reproduce It Yourself](#reproduce-it-yourself)
10. [Publications](#publications)
11. [License](#license)

---

## Overview

TFPT derives the fine-structure constant α and other fundamental constants from **pure mathematical constraints**—topology, geometry, and quantum consistency—rather than experimental measurements.

**Key characteristics:**
- **Zero free parameters** — The kernel is parameter-free; all policies are explicit and auditable
- **Discrete grammar** — c₃, φ₀, b₁ from topology; k=2 (double cover), g=5 (holonomy degeneracy)
- **Sub-ppm precision** — Two-defect result achieves z = 1.87 vs CODATA 2022

### The Three Fundamental Invariants

| Symbol | Definition | Value | Origin |
|--------|------------|-------|--------|
| **c₃** | 1/(8π) | 0.039788... | 11D Chern-Simons quantization |
| **φ₀** | φ_tree + δ_top | 0.053172... | Möbius fiber geometry |
| **b₁** | 41/10 | 4.1 (exact) | SM spectral trace coefficient |

Where:
- φ_tree = 1/(6π) ≈ 0.0530516 (η-invariant gluing on orientable double cover)
- δ_top = 3/(256π⁴) ≈ 0.000120 (topological correction from orbifold sectors)

---

## The Core Derivation

### The 5-Step Derivation

From topology to α — each step builds deterministically on the previous:

| Step | Domain | What It Fixes | Result |
|------|--------|---------------|--------|
| **1. Topology** | 11D Chern-Simons | Topological coupling | c₃ = 1/(8π) |
| **2. Geometry** | Double cover + η-gluing | Geometric scale | φ₀ = 1/(6π) + 3/(256π⁴) |
| **3. SM Input** | Spectral trace | UV coefficient | b₁ = 41/10 |
| **4. UV Log** | Reference scale | Log coefficient | K(α) = (b₁/2π) ln(1/φ₀(α)) |
| **5. Closure** | Fixed-point equation | Self-consistent α | CFE → α⁻¹ = 137.035999 |

### The Cubic Fixed-Point Equation (CFE)

```
α³ − 2c₃³α² − 8b₁c₃⁶ ln(1/φ₀(α)) = 0
```

With the self-consistent backreaction:

```
φ₀(α) = φ_tree + δ_top·e^(−kα) + δ₂·e^(−2kα)
```

Where k=2 (double cover) and δ₂ = (g/4)·δ_top² with g=5 (SU(5) holonomy degeneracy).

### Three-Tier Improvement

| Level | Formula | α⁻¹ | Deviation |
|-------|---------|-----|-----------|
| **Baseline** | Fixed φ₀ | 137.036501 | +3.67 ppm |
| **One-defect** | φ₀(α) with k=2 | 137.0359941 | −0.037 ppm |
| **Two-defect (g=5)** | + δ₂ term | 137.035999216 | z = 1.87 vs CODATA |

---

## The 8 Axioms

TFPT rests on 8 minimal assumptions (7 established, 1 new postulate):

### Geometric (4)

| ID | Axiom | Status | Output |
|----|-------|--------|--------|
| **A1** | 11D Chern-Simons Base | Established | c₃ = 1/(8π) |
| **A2** | Orientable Double Cover | New Postulate | k = 2 multiplicity |
| **A3** | Möbius Fiber Topology | Established | φ_tree = 1/(6π) |
| **A4** | Discrete Curvature Quantization | Established | δ_top = 3/(256π⁴) |

### QFT (2)

| ID | Axiom | Status | Output |
|----|-------|--------|--------|
| **A5** | Fujikawa Anomaly Measure | Established | Anomaly normalization |
| **A6** | Background-Field Ward Identities | Established | RG scheme independence |

### Closure (2)

| ID | Axiom | Status | Output |
|----|-------|--------|--------|
| **A7** | Fixed-Point Stationarity | Established | Cubic equation for α |
| **A8** | GR Limit Requirement | Established | Consistency with GR |

---

## Key Results

### Primary Result: Fine-Structure Constant

| Quantity | TFPT (g=5) | CODATA 2022 | Status |
|----------|------------|-------------|--------|
| **α⁻¹(0)** | 137.035999216 | 137.035999177(21) | z = 1.87 |
| **ᾱ⁻¹₍₅₎(M_Z)** | 127.940519 | 127.93 ± 0.008 | z ≈ 1.31 |

### Observable Predictions

| Prediction | Value | Measurement | Status |
|------------|-------|-------------|--------|
| **Cosmic birefringence β** | 0.2424° | 0.30° ± 0.11° (Planck PR4) | Solid (0.5–1.7σ) |
| **Cabibbo angle λ** | 0.2245 | 0.2243 ± 0.0005 | Solid (−0.24%) |
| **RG fingerprint α₃(1 PeV)** | ≈ φ₀ | — | Solid (1.38% dev) |
| **Inflation r** | ≈ 0.0038 | — | Open (R² completion) |
| **Axion mass** | 65.2 μeV | — | Testable (ν = 15.76 GHz) |
| **PMNS sin²θ₁₃** | 0.0231 | 0.0220 ± 0.0007 | Testable (TM1) |

### v2.5 Crosslink Layer

| Parameter | Value | Origin |
|-----------|-------|--------|
| **g** (holonomy degeneracy) | 5 | SU(5) hypercharge eigenspace |
| **δ₂** (two-defect weight) | (5/4)δ_top² ≈ 1.81×10⁻⁸ | Defect partition enumeration |
| **γ(0)** (damping anchor) | g/(g+1) = 5/6 | E₈ cascade |

---

## Repository Structure

```
TFPT-Github/
├── readme.md                    # This file
│
├── latex/                       # LaTeX paper sources
│   ├── tfpt-theory-fullv25.tex  # Latest standalone paper (v2.5)
│   ├── tfpt-theory-fullv24.tex  # Previous version (v2.4)
│   ├── tfpt-v25-snippets.tex    # Auto-extracted numerics
│   ├── progress.md              # Paper revision history
│   └── README.md                # Build instructions
│
├── tfpt-suite/                  # Verification suite (~90 modules)
│   ├── run_suite.py             # CLI runner
│   ├── tfpt_suite/              # Core Python library
│   │   ├── modules/             # 39+ verification modules
│   │   └── data/                # Reference data (JSON)
│   ├── out/                     # Engineering-mode outputs
│   ├── out_physics/             # Physics-mode outputs
│   ├── theoryv3/                # Discrete building-blocks analysis
│   ├── unconventional/          # Search/audit tools
│   ├── tfpt-test-results.pdf    # Engineering report
│   ├── tfpt-test-results-physics.pdf  # Physics report
│   └── README.md                # Suite documentation
│
├── E8ChainSolver/               # E8 nilpotent orbit analysis
│   └── e8_orbit_engine/         # γ(n) damping validation
│       ├── src/                 # Python source
│       ├── data/                # Orbit data
│       └── README.md            # Usage
│
├── Pyrate3/                     # PyR@TE3 RGE integration
│   ├── models/                  # YAML model definitions
│   └── pyrate/                  # PyR@TE library
│
└── old-latex-and-docs/          # Historical documents
    └── README.txt
```

---

## Getting Started

### Prerequisites

- Python 3.9+
- TeX distribution (TeX Live / MacTeX)
- Optional: CAMB for CMB calculations

### Installation

```bash
# Clone the repository
git clone https://github.com/sthamann/tfpt-theoryv2
cd tfpt-theoryv2

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r tfpt-suite/requirements.txt
```

### Running the Suite

```bash
# List all modules
python3 tfpt-suite/run_suite.py list-modules

# Run all modules (engineering mode)
python3 tfpt-suite/run_suite.py run-all

# Run in physics mode (stricter)
python3 tfpt-suite/run_suite.py --mode physics run-all

# Build PDF report
python3 tfpt-suite/run_suite.py build-report
```

### Compiling the Paper

```bash
cd latex
pdflatex -interaction=nonstopmode tfpt-theory-fullv25.tex
pdflatex -interaction=nonstopmode tfpt-theory-fullv25.tex
```

---

## Components

### tfpt-suite (Verification Suite)

Modular, reproducible Python framework turning TFPT claims into testable computations.

| Category | Purpose | Key Modules |
|----------|---------|-------------|
| **Core Invariants** | c₃, φ₀, β, CFE | `core_invariants`, `discrete_consistency_uniqueness` |
| **α Sector** | Precision audit | `alpha_precision_audit`, `alpha_on_shell_bridge` |
| **RG & Gauge** | Two-loop running | `two_loop_rg_fingerprints`, `unification_gate` |
| **Cosmology** | Bounce, inflation | `bounce_perturbations`, `k_calibration` |
| **Flavor** | CKM/PMNS | `ckm_full_pipeline`, `pmns_full_pipeline` |
| **Dashboard** | Global tests | `global_consistency_test`, `predictions_dashboard` |

### E8 Orbit Engine

Validates γ(n) ≈ 0.834 + 0.108n + 0.0106n² from E8 nilpotent orbits:

```bash
e8-orbit fit data/nilpotent_orbits.csv
```

### theoryv3 (Analysis Branch)

Discrete building-blocks analysis: π → invariants → g=5 → α.

```bash
python3 tfpt-suite/theoryv3/run_theoryv3.py run-all
```

---

## Predictions & Falsification

### Direct Falsification Paths

| Claim | Prediction | Falsification Test | Priority |
|-------|------------|-------------------|----------|
| Fine-structure constant | α⁻¹ = 137.035999216 | \|z\| > 5 vs CODATA | Critical |
| Cosmic birefringence | β = 0.2424° | \|z\| > 5 vs CMB | Critical |
| RG fingerprint | α₃(1 PeV) ≈ φ₀ | \|z\| > 3 at colliders | High |
| Axion signal | ν = 15.76 GHz | No signal at ADMX-G2 | High |
| Inflation r | r ≈ 0.0038 | \|z\| > 5 from CMB-S4 | Medium |
| PMNS mixing | sin²θ₁₃ ≈ 0.0231 | \|z\| > 5 vs oscillation | Medium |

### Upcoming Decisive Tests

- **LiteBIRD / CMB-S4**: Cosmic birefringence, inflation r
- **ADMX-G2 / MADMAX**: Axion at 15.76 GHz
- **Hyper-Kamiokande / DUNE**: PMNS precision

---

## Reproduce It Yourself

### Two-Defect (g=5) — Primary Result

```python
import math

c3 = 1 / (8 * math.pi)
phi_tree = 1 / (6 * math.pi)
delta_top = 3 / (256 * math.pi**4)
b1 = 41 / 10
k = 2   # double cover
g = 5   # SU(5) holonomy degeneracy
delta2 = (g / 4) * delta_top**2

alpha = 1 / 137  # Initial guess
for _ in range(10):  # Self-consistency loop
    phi0 = phi_tree + delta_top*math.exp(-k*alpha) + delta2*math.exp(-2*k*alpha)
    for _ in range(50):  # Newton iteration
        const = 8 * b1 * c3**6 * math.log(1/phi0)
        f = alpha**3 - 2*c3**3*alpha**2 - const
        df = 3*alpha**2 - 4*c3**3*alpha
        alpha -= f / df

print(f"alpha_inv = {1/alpha:.12f}")  # 137.035999216158
```

**Result:** α⁻¹ = 137.035999216 (z = 1.87 vs CODATA 2022)

### Falsify Without Believing

1. Run the script above
2. Try changing any constant — the result breaks
3. Try adding a fit parameter — there's nowhere to put it

---

## Publications

### Main Paper

**"Topological Fixed Point Theory (TFPT): A Parameter-Free Derivation of the Fine-Structure Constant from Topology, Geometry, and Quantum Consistency"**

- **Version:** 2.5 (January 27, 2026)
- **Authors:** Stefan Hamann & Alessandro Rizzo
- **Download:** [Zenodo](https://zenodo.org) | [GitHub](https://github.com/sthamann/tfpt-theoryv2)

### Resources

- **Website:** [www.fixpoint-theory.com](https://www.fixpoint-theory.com)
- **GitHub:** [github.com/sthamann/tfpt-theoryv2](https://github.com/sthamann/tfpt-theoryv2)
- **Full Suite Report:** `tfpt-suite/tfpt-test-results-physics.pdf`

---

## Why It Matters

- **Constants are no longer arbitrary inputs** — They emerge from pure mathematical constraints
- **α emerges as a forced solution** — Uniquely determined by topological and geometric fixed points
- **ppm-precision without parameters** — No adjustable quantities
- **Testable by upcoming experiments** — CMB-S4, ADMX-G2, Hyper-K, LiteBIRD

---

## License

MIT License

---

## Citation

```bibtex
@article{hamann2026tfpt,
  title={Topological Fixed Point Theory: A Parameter-Free Derivation of the 
         Fine-Structure Constant from Topology, Geometry, and Quantum Consistency},
  author={Hamann, Stefan and Rizzo, Alessandro},
  year={2026},
  version={2.5}
}
```

---

*Last updated: 2026-01-28*
