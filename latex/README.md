# LaTeX sources

- `tfpt-theory-fullv25.tex`: TFPT paper snapshot (**Version 2.5, 27 Jan 2026**). **Standalone**: suite snippets/manifests/reference-ledger are inlined (no `\input{...}` required) and **all plots are LaTeX-native** (TikZ/pgfplots; no PNG `\includegraphics`).
- `tfpt-theory-fullv24.tex`: TFPT paper snapshot (**Version 2.4, 22 Jan 2026**).
- `tfpt-theory-fullv21.tex`: TFPT paper snapshot (Version 2.2, 21 Jan 2026).
- `tfpt-v25-snippets.tex`: v2.5 numeric/tables/plots pack (auto-extracted from `tfpt-suite` reports; intended to be `\input{}`-ed by the paper source used for v2.5).

## Changes applied (v2.1 hardening pass)

- Marked the Möbius **seam-as-boundary** step as an explicit **Postulate A3** (instead of “it follows”), added a short plausibility note, and added **Appendix F** with a clear “what remains to be formalized” ledger.
- Removed “magic number” framing for `48` by rewriting `\delta_{\mathrm{top}}` as an **orbifold deficit sum rule** (and updated the coefficient appendix accordingly).
- Tightened the `b_1` wording: **hypercharge-trace UV log coefficient**, and the CFE as a **stationarity/closure condition**, not an RG differential equation; added a reviewer-safe disclaimer.
- Added a compact **`\beta_{\text{rad}}` dictionary** and fixed the degrees conversion in the `\beta_{\text{rad}}` identity summary.
- Added an **Identity Catalogue** appendix with status tags (**D/K/L/P**) so “definitions vs consequences vs predictions” are explicit.
- Unified the **Cabibbo** formula and added the exactly equivalent `\beta_{\text{rad}}` form.
- Updated the **high-precision constants table** to match the recomputation from the canonical definitions.
- Closed the “single anchor” presentation of `\gamma(0)` by setting **`\gamma(0)=5/6`** from already-used discrete normalizations (SU(5) and Möbius boundary counting), and updated the E$_8$/PMNS text accordingly.
- Replaced the birefringence **`\Delta a=\varphi_0`** phrasing as a bare postulate by a minimal **late-time shifted axion attractor** dynamical completion (with a new tomography test `\beta(z)`).
- Reworked the inflation section to an **$R^2$ (Starobinsky) completion** with **`M/\bar M_{\rm Pl}=\sqrt{8\pi}\,c_3^4`**, making **$A_s,n_s,r$** explicit and testable.
- Removed the undefined per-block prefactor `G_B` by adopting the minimal choice **`G_B = 1`** in the block constant formula (no hidden knobs).

## Changes applied (v2.2 hardening pass)

- Replaced “varying $\alpha$” with **modulus closure**: a coupling modulus $\sigma$ and $\delta\Gamma_{\text{eff}}/\delta\sigma=0$.
- Reframed **$b_1$** as a **spectral trace coefficient** from APS boundary determinants (Appendix boundary).
- Replaced Postulate A3 with **$\eta$-gluing seam lemma** (explicit spectral flow) and rewrote Appendix F.
- Reframed **$\delta_{\text{top}}$** via the **spin-lifted $\mathbb{Z}_2$ deficit factor** and four-sector seam symmetry.
- Replaced the linear backreaction with the **exact exponential closure** $\phiz(\alpha)=\phitree+\deltatop e^{-2\alpha}$.
- Added a **geometric renormalization lemma**: $\ln(1/\phiz)=\ln(\Lambda/\mu_{\text{geo}})$.
- Added **cover-degree scaling** and **$c_2$ next-order response** to treat the ppm residual systematically.
- Split the axion sector into **$a_{\text{top}}$ vs $a_{\text{QCD}}$** and proved $\theta_i=\pi$ as a topological sector.
- Added **$R^2$ derivation (Sec. 8.2–8.3, Appendix K)** and **bounce perturbations (Sec. 8.4, Appendix L)**.
- Updated claim map, assumption ledger, status table, verification suite modules, and references (APS/$\eta$-gluing).

## Changes applied (v2.3 hardening pass)

- Replaced backreaction “modeling” with a defect-sector partition-function lemma; exponent $2$ now fixed by deck-invariant normalization.
- Clarified $\alpha$ wording as a **closed-form defect-sector reweighting** consequence in claim map/status tables.
- Counted $\delta_{\text{top}}$ sectors from cut-geometry patches $(C_{1,2},\Gamma_\pm)$; labeled the “48” mnemonic as such.
- Quantized $\Delta a_{\text{top}} = n\varphi_0$ (minimal $n=1$) and marked $s(t)$ as an optional dynamical realization.
- Updated axion DM relic density with $C_{\mathrm{str}}$ and explicit pre/post-inflation branches ($N_{\mathrm{DW}}=1$), including the two $\theta_i$ cases.
- Labeled inflation/bounce claims as conditional (K1–K4 list + pending modules) and added a bounce status note.
- Added a K-theory note linking spectral flow to $\mathrm{wind}(\det U_\Gamma)$ in Appendix F.

## Changes applied (v2.4 hardening pass) (Jan 2026)

- Fixed the **Maxwell convention** explicitly as $F=dA$ (independent of the torsionful connection), added in Conventions + UFE section and recorded in the Assumption Ledger.
- Clarified that the Maxwell equations use the **Levi-Civita covariant derivative** ($\nabla^{\mathrm{LC}}$) in the UFE theorem statement.
- Reframed the UFE theorem as **full-connection form first**, then the **Levi-Civita + contorsion rewrite** for computations; added explicit torsion/contorsion definitions.
- UFE notation cleanup: made the displayed gravitational equation explicitly $\hat{G}_{AB}(\hat{\Gamma})$ (no $\Gamma/\hat{\Gamma}$ mismatch), stated $\hat{R}(\hat{\Gamma})$ is the Ricci scalar of the full connection, and removed ambiguity about which $\nabla$ is used when expanding in Levi-Civita + contorsion variables.
- RG fingerprints wording: replaced “measured” by **two-loop extrapolated** and stated the PDG $M_Z$ boundary-condition origin.
- Simplified the $\xi$ definition by removing a redundant identity.
- Clarified the CFE root figure caption as **baseline with fixed $\varphi_0$**.
- Added the **$g_{a\gamma\gamma}$ vs.\ $c_3$ glue** in the birefringence transport law and recorded it in the Identity Catalogue.
- Updated frontmatter to **Version 2.4 (22 January 2026)**.
- Updated `tfpt-theory-fullv24.tex` status dashboard to reflect current suite closure: **unification gate now passes** under a discrete policy scan (mismatch \(\sim 0.38\%\)), and a new “suite-closure (physics-mode)” block summarizes CKM/PMNS/global χ², DM/Λ/BBN/GW/QED/arrow proxy closures.
- Updated the axion DM relic-density narrative to the shipped default: post-inflation RMS branch with a discrete strings/domain-walls factor \(C_{\mathrm{str}}=7/3\) (no longer the “0.02\% axion fraction” placeholder).

## Changes applied (v2.5 reproducibility upgrade) (Jan 2026)

- Centralized numerics/tables/plots via `\input{tfpt-v25-snippets.tex}` to avoid number drift across the paper.
- Added a **Kernel + Compiler** frontmatter section (“Kernel and Discrete Grammar”), including a **constant-factory pipeline** diagram and a **derived constants** table with per-module evidence pointers.
- Added a **crosslink map** making **holonomy degeneracy \(g=5\)** an explicit global node and propagated it into the paper where relevant:
  - \(\delta_2=(g/4)\delta_{\text{top}}^2\) with \(g=5\) (two-defect sector)
  - \(\gamma(0)=g/(g+1)=5/6\) (E\(_8\) damping anchor)
- Added an explicit **Comparison and Scheme Policy** section (primary reference: \(\bar\alpha^{(5)}(M_Z)\) under an MS\(\overline{\mathrm{MS}}\)-at-\(M_Z\) policy; \(\alpha(0)\) is treated as diagnostic IR/on-shell quantity once scheme/matching is explicit).
- Fixed a TeX control-sequence footgun: renamed the `A_s` snippet macros to avoid digits in macro names (prevents stray `55...` tokens leaking into page 1 output).
- Fixed Abstract/claim-map wording for $\alpha$: one-defect truncation vs **derived two-defect** \(g=5\) value; explicit separation of **ppm** vs **strict z-score**.
- Added **alpha precision audit** content (second-order defect term $\delta_2$, sensitivity plot, and benchmark table).
- Updated the $\beta_{\rm rad}$ identity table: added a numeric **$\Omega_b$ benchmark** plus comparison table and updated the status dashboard accordingly.
- Recorded **$R^2$ effective-action benchmarks** and Starobinsky $A_s$ values in Appendix K (conditional on K1–K4, per module notes).
- Extended Appendix L with the **bounce transfer-function plot** (using direct `coordinates` for robustness), Wronskian diagnostic note, and the **$k\to\ell$ calibration** table.
- Added/updated compact benchmarks for **APS $\eta$-gluing** (worked example) and the **M\"obius cusp classification** statement.
- Added v2.5 benchmark inserts for **CKM** (baseline pipeline matrix + UV Yukawa texture block + Yukawa heatmap) and **PMNS Z$_3$ breaking** (permutation + variant table).
- Added a **present-day torsion regimes** subsection and a **torsion bounds** table (explicitly separating local-vacuum vs cosmological-torsion assumptions).
- Added a **dark energy (Λ) candidate** subsection (defect-suppressed torsion condensate / discrete normalization scan; explicit “mechanism interface” note).
- Updated the **verification suite summary** to the current `tfpt-suite` module structure and added the **global χ² dashboard** note including the “excluding α” total (alpha-dominance transparency).
- Eliminated remaining reviewer-visible consistency traps:
  - α is now single-source across theorems/appendices: **two-defect \(g=5\)** is the headline prediction; **one-defect** is explicitly a diagnostic truncation baseline.
  - Upgraded the α sensitivity figure legend to show **\(\delta_2=0\)** vs **\(g=4\)** vs **\(g=5\)** cases.
  - Added an auto-generated **Prediction Ledger and Falsification** appendix (exported by `tfpt-suite/run_suite.py export-prediction-ledger`, and inlined into the standalone snapshot).
  - Made the CKM benchmark hygiene explicit: the suite currently uses an **external \(\delta_{\rm CKM}\)** override for benchmarking, and this is stated in the paper.
  - Marked the flagged low-\(\hat{k}\) bounce-transfer points visually (Wronskian diagnostic), to avoid “numerics shaky” reviewer impressions.

## Build

Use your local TeX distribution (e.g. TeX Live / MacTeX):

```bash
# Standalone v2.5 snapshot: compile directly (no external TeX inputs needed)
pdflatex -interaction=nonstopmode -halt-on-error tfpt-theory-fullv25.tex
pdflatex -interaction=nonstopmode -halt-on-error tfpt-theory-fullv25.tex
pdflatex -interaction=nonstopmode -halt-on-error tfpt-theory-fullv24.tex
pdflatex -interaction=nonstopmode -halt-on-error tfpt-theory-fullv24.tex
```

