# Progress (TFPT paper edits)

## v2.1 hardening pass (Jan 2026)

- [x] **Postulate A3 clarity**: label the seam/effective-boundary step explicitly; add short plausibility note; add Appendix F “open formalization” ledger.
- [x] **`48` coefficient**: remove numerology vibe by presenting `\delta_{\mathrm{top}}` as an orbifold deficit sum rule; update the coefficient appendix accordingly.
- [x] **`b_1` precision**: remove “topological protection” phrasing for `b_1`; emphasize hypercharge-trace UV log coefficient and CFE as stationarity (not RG running).
- [x] **`\beta_{\text{rad}}` dictionary**: add a small dictionary box; fix the degrees conversion in the identity summary.
- [x] **Identity status**: add appendix legend D/K/L/P and tag representative identities.
- [x] **Cabibbo consistency**: ensure the factor sits outside the square root; add equivalent `\beta_{\text{rad}}` form.
- [x] **Numerical constants**: align the high-precision constants table with the canonical recomputation.

## v2.1 closure of open anchors (Jan 2026)

- [x] **`γ(0)` without physical anchor**: set `γ(0)=5/6` from existing SU(5) + Möbius normalizations; update E$_8$ cascade wording and downstream PMNS usage.
- [x] **Dynamic `Δa` (no bare postulate)**: replace “Δa = φ0” with a late-time shifted cosine/attractor completion and add the `β(z)` tomography signature.
- [x] **Inflation completion**: switch to $R^2$ (Starobinsky) with `M/\bar M_{\rm Pl}=\sqrt{8\pi}\,c_3^4` and make $A_s,n_s,r$ explicit; mark the remaining derivation as an open calculation task.

## v2.2 hardening pass (Jan 2026)

- [x] **Modulus closure (A6)**: introduce coupling modulus $\sigma$ and replace “vary $\alpha$” with $\delta\Gamma_{\text{eff}}/\delta\sigma=0$, equivalent to $\partial U/\partial\alpha=0$.
- [x] **$b_1$ spectral index**: tie the hypercharge trace to the APS boundary determinant coefficient; update Appendix boundary derivation and reviewer-safe statement.
- [x] **Seam term via $\eta$-gluing**: replace Postulate A3 with a gluing/spectral-flow lemma; explicit $D_\Gamma$, $U_\Gamma$, and $\mathrm{SF}$ computation in Appendix F.
- [x] **$\delta_{\text{top}}$ via spin lift**: replace orbifold folklore with the spin-lifted $\mathbb{Z}_2$ deficit factor and four-sector seam symmetry; update Appendix 48 derivation.
- [x] **Backreaction closure**: replace linear form by exact $\phiz(\alpha)=\phitree+\deltatop e^{-2\alpha}$ and fix the full series (no free coefficients).
- [x] **Two-field axion sector**: separate $a_{\text{top}}$ (birefringence) and $a_{\text{QCD}}$ (DM); add topological $\theta_i=\pi$ sector argument.
- [x] **Inflation derivation + bounce**: add Sections 8.2–8.4 with derivation statements; add Appendix K (effective action) and Appendix L (bounce perturbations).
- [x] **Frontmatter/status updates**: claim map, assumption ledger, status table, verification suite modules, geometric-log lemma, and $\Omega_b$ conjecture wording.
- [x] **References**: add APS and $\eta$-gluing citations.

## v2.3 hardening pass (Jan 2026)

- [x] **Defect-sector reweighting lemma**: derive $\delta_{\text{top}}(\alpha)$ from partition logic.
- [x] **$\alpha$ status wording**: claim map/status use closed-form defect reweighting.
- [x] **Sector count $F=4$**: cut-geometry patch counting in Appendix G + A4.
- [x] **Axion DM accounting**: add $C_{\mathrm{str}}$ and two $\theta_i$ scenarios.
- [x] **$\Delta a_{\text{top}}$ quantization**: define $n\phiz$; $s(t)$ as optional realization.
- [x] **``48'' origin consistency**: label mnemonic and point to Appendix G.
- [x] **Inflation/bounce conditioning**: K1--K4 list + pending modules.
- [x] **K-theory note**: $\mathrm{SF}(U_\Gamma)$ as winding of $\det U_\Gamma$.
- [x] **PQ scenario branches**: pre/post inflation + $N_{\mathrm{DW}}=1$ note.
- [x] **Axiom appendix links**: explicit Appendix F/G references.

## v2.4 hardening pass (Jan 2026)

- [x] **Maxwell convention**: state explicitly $F=dA$ (independent of torsionful connection); add in Conventions + UFE section; record in Assumption Ledger.
- [x] **Maxwell $\nabla$ clarity**: denote Levi-Civita derivative explicitly ($\nabla^{\mathrm{LC}}$) in the Maxwell equations.
- [x] **UFE theorem framing/notation**: define $\hat{G}_{AB}(\hat{\Gamma})$ explicitly (full torsionful connection), stabilize $\Gamma/\hat{\Gamma}$ usage, and clarify Levi-Civita derivatives when expanding in $\Gamma(\hat{g})+K$ variables; define torsion and contorsion.
- [x] **UFE action clarity**: state explicitly that $\hat{R}(\hat{\Gamma})$ is the Ricci scalar of the full connection.
- [x] **Birefringence coupling glue**: add explicit $d\beta/d\eta = -(1/2)g_{a\gamma\gamma}\,da_{\text{top}}/d\eta = 2c_3\,da_{\text{top}}/d\eta$ line; record in Identity Catalogue.
- [x] **RG fingerprints wording**: replace “measured” with “two-loop extrapolated” and note the PDG $M_Z$ boundary-condition origin.
- [x] **$\xi$ cleanup**: remove redundant identity in the $\xi$ definition.
- [x] **CFE figure caption clarity**: mark $\alpha^{-1}\approx 137.0365$ as baseline with fixed $\varphi_0$.
- [x] **Version bump**: update frontmatter to Version 2.4 (22 January 2026).
- [x] **Suite sync (Jan 2026)**: update reproducibility/status wording to reflect that `effective_action_r2` and `bounce_perturbations` are implemented in `tfpt-suite/`, and clarify the CKM disposition line accordingly.

## v2.5 reproducibility upgrade (Jan 2026)

- [x] **Version bump**: update frontmatter to Version 2.5 (27 January 2026).
- [x] **Numerics centralization**: include `tfpt-v25-snippets.tex` and replace hard-coded core numbers with macros.
- [x] **Kernel / compiler structure**: add “Kernel and Discrete Grammar” + constant-factory table + crosslink map (holonomy degeneracy \(g=5\)).
- [x] **Policy layer**: add explicit “Comparison and Scheme Policy” section (primary reference \(\bar\alpha^{(5)}(M_Z)\); \(\alpha(0)\) diagnostic once matching is explicit).
- [x] **TeX macro name fix**: rename the `A_s` macros to avoid digits in control sequence names (prevents “552.016…” style artifacts).
- [x] **Abstract + claim-map hygiene**: fix $\alpha$ wording; state baseline + refined value and separate ppm vs strict z-score.
- [x] **Alpha precision audit**: add $\delta_2$ benchmark table + sensitivity plot (and record the effective $k$).
- [x] **$\beta_{\text{rad}}$ table / $\Omega_b$**: insert numeric benchmark + comparison table; update status dashboard.
- [x] **$R^2$ appendix benchmarks**: record effective-action coefficients and $A_s$ values; add conditional status note (K1–K4).
- [x] **Bounce appendix**: add transfer-function plot + Wronskian warning and the $k\to\ell$ calibration table (and make the plot robust by using direct `coordinates`).
- [x] **APS $\eta$-gluing / cusps**: add worked example note for $\eta(0,a)$ and update cusp classification statement + scan note.
- [x] **CKM / PMNS**: insert CKM baseline pipeline matrix benchmark + UV Yukawa texture block (+ Yukawa heatmap); add PMNS Z$_3$ breaking permutation + variant table; update wording to conditional.
- [x] **Global consistency**: add $\chi^2$ dashboard note explaining alpha-dominance and the “excluding alpha” total.
- [x] **Torsion bounds**: add bounds table and clarify local-vacuum vs cosmological-torsion regimes.
- [x] **Dark energy (Λ) candidate**: add discrete scan / torsion-condensate interface note (explicit mechanism gap).
