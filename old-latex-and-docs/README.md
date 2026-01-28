# Old LaTeX Documents and PDFs

This folder contains archived LaTeX source files and PDF documents from earlier stages of TFPT (Topological Fixed Point Theory) development. These files were extracted from original PDF documents and serve as historical reference material.

## Compilation

All `.tex` files require **LuaLaTeX** or **XeLaTeX** due to Unicode math symbols in the extracted text.

```bash
lualatex <filename>.tex
```

## File Index

### Main Paper

| File | Source PDF | Description |
|------|-----------|-------------|
| `paper_v1_06_01_09_2025.tex` | Paper V1.06 - 01.09.2025.pdf | Comprehensive 81-page paper covering the complete TFPT framework (German: "Von Topologie zu Dynamik: Die Ordnung hinter α und den Naturkonstanten"). Includes genetic algorithm development, E₈ cascade, inflation, and Standard Model derivations. |

### Theory Notes

| File | Source PDF | Description |
|------|-----------|-------------|
| `theory_root.tex` | Theory-Root.pdf | **Structural Equivalence of Fields and a Unified Fixed-Point Equation** — Self-contained presentation showing how c₃ = 1/(8π) and φ₀ yield the cubic fixed point for α, the UFE with torsion, and cosmic birefringence prediction. Includes TikZ visualizations. |
| `unified_field_equation.tex` | Unified_Field_Equation.pdf | **On the Unified Field Equation as a Consequence of a Quantum Fixed-Point Condition** — Derivation of the UFE from the cubic condition, including the birefringence law β = Δa/(4π) and numerical comparison with Planck PR4 data. |
| `update_tfptv1_07sm.tex` | Update TFPTV1.07SM.pdf | **TFPT Complete Standard Model Derivations (v1.0.8)** — Integration of CFE with geometric self-consistency, E₈ log-exact cascade, Z₃ flavor architecture, and mixing matrix derivations (CKM/PMNS). |

### Technical Notes

| File | Note ID | Source PDF | Description |
|------|---------|-----------|-------------|
| `eliminating_k.tex` | H1 | Eliminating K.pdf | **κ² from φ₀ and c₃** — Shows that gravitational coupling κ in the UFE is determined by φ₀/c₃², removing G as an independent parameter. Derives ξ ≈ 0.748328 (close to 3/4 at tree level). |
| `improvedcubic.tex` | — | ImprovedCubic.pdf | **Geometric Self-Consistency on the Orientable Double Cover** — Derivation of the backreaction formula φ₀(α) = φ_tree + δ_top(1−2α), improving α⁻¹ agreement with CODATA from +3.67 ppm to −0.064 ppm. |
| `five_problems.tex` | S1 | Five Problems.pdf | **Solutions to Foundational Problems in Physics** — Addresses UV completion via Asymptotic Safety, singularity resolution via torsion bounce, black holes as spacetime vortices, dark matter (TFPT axion), and dark energy origins. |

### Supplementary

| File | Description |
|------|-------------|
| `TFPT Complete Proof - Fine Structure Constant Derivation.pdf` | PDF document (original, not LaTeX) |
| `README.txt` | Original plain-text readme with file listing |

## Key Invariants (Reference)

The documents consistently use:

- **Topological invariant:** c₃ = 1/(8π) ≈ 0.0397887
- **Geometric scale (tree):** φ_tree = 1/(6π) ≈ 0.0530516
- **Topological surcharge:** δ_top = 3/(256π⁴) ≈ 0.0001204
- **Combined baseline:** φ₀ = φ_tree + δ_top ≈ 0.0531720
- **Self-consistent:** φ₀(α) = φ_tree + δ_top(1−2α) ≈ 0.0531702

## Relationship to Current Codebase

These documents represent the theoretical foundation that the `tfpt-suite` numerical implementation is based on. The suite validates these derivations computationally—see `/tfpt-suite/README.md` for details on the test infrastructure.
