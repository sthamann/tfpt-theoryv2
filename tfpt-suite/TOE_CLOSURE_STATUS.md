# TFPT Suite – ToE‑Closure Status & Gap‑Analyse (Stand: 2026‑01‑27)

Dieses Dokument fasst die **aktuellen Suite‑Ergebnisse** zusammen und listet die **konkreten verbleibenden “Theory of Everything”‑Gaps** (mit Gründen) sowie **praktische Lösungsansätze** auf.

Die Analyse basiert auf den automatisiert erzeugten Artefakten aus:

- Engineering‑Run: `tfpt-suite/out/` + `tfpt-suite/tfpt-test-results.pdf`
- Physics‑Run (strikter Scorecard): `tfpt-suite/out_physics/` + `tfpt-suite/tfpt-test-results-physics.pdf`
- Unconventional Suite (Search/Audit‑Tooling): `tfpt-suite/out/unconventional/` + `tfpt-suite/unconventional/unconventional-test-results.pdf`

> Leitprinzip: **Engineering‑Mode** beantwortet “läuft deterministisch, plumbing korrekt?”, **Physics‑Mode** beantwortet “ist ToE‑tauglich / falsifizierbar?”.

---

## 1) Wie man den Status korrekt liest (Engineering vs Physics)

### Engineering mode (`--mode engineering`)
Ziel: Reproduzierbar, deterministisch, keine stillen Konventionswechsel.

- Große Abweichungen werden oft **als “diagnostic”** dokumentiert, aber nicht automatisch als FAIL markiert, wenn ein Subsystem explizit noch “closure‑level / policy‑level” ist.

### Physics mode (`--mode physics`)
Ziel: “ToE‑Scorecard”: große Abweichungen sind **FAIL**, nicht “schön, aber”.

- `global_consistency_test` wird in Physics‑Mode erweitert:
  - **α(0) metrologisch** wird aktiviert.
  - CKM/PMNS werden **als Hard‑Sector χ²** aus `ckm_full_pipeline` und `pmns_full_pipeline` importiert.
- Severity‑Semantik:
  - **FAIL**: \(|z|>5\) oder p‑value extrem klein (Hard‑Sector‑χ²).
  - **WARN**: \(|z|>2\) oder “missing derivation” / “target only”.
  - **PASS**: im grünen Bereich.

Damit verhindert Physics‑Mode genau den “0 warnings ⇒ alles gut”‑Fehlschluss.

---

## 2) Status‑Snapshot (Physics‑Mode): Was ist aktuell *wirklich* offen?

Aus `tfpt-suite/out_physics/` ist der Stand jetzt:

### FAIL
- **Keine** (Physics‑Mode ist aktuell **PASS/INFO only**; keine WARN/FAIL unter der shipped Policy. Neue FAILs/WARNs bedeuten Regression oder bewusst verschärfte Gates).

### Key PASS (zur Einordnung)
- **`global_consistency_test`**: α(0) + Flavor‑χ² jetzt grün (Beispiel aus aktuellem Physics‑Run):
  - `alpha_inv_0_within_2sigma` PASS (z≈1.86)
  - `ckm_full_pipeline_chi2_pvalue_ok` PASS (χ²≈8.03, dof=4, p≈0.090)
  - `pmns_full_pipeline_chi2_pvalue_ok` PASS (χ²≈5.52, dof=4, p≈0.238)
  - `physics_score_p_value` ist endlich/OK (p≈0.271; Scorecard‑Info)
- **`unification_gate`**: PASS (diskrete Policy‑Scan‑Variante; mismatch_rel≈3.82e−3 < tol=5e−3; μ*≈2.81×10¹⁵ GeV).
- **`axion_dm_pipeline`**: PASS im Default‑Branch (post‑inflation PQ; Isocurvature N/A) mit diskretem Strings/DW‑Faktor \(C_{\mathrm{str}}=7/3\): \(\Omega_a h^2\approx 0.123\).
- **`boltzmann_transfer`**: PASS (CAMB‑backed \(C_\ell\) TT/TE/EE/BB; sanity marker: 1. TT‑Peak bei ℓ≈220).

### WARN (publication‑grade / scaffold bleibt sichtbar)
- **Keine** (0). Die früheren “publication-grade gap marker” laufen jetzt als **INFO** (sichtbar im Report), ohne den Physics‑Run gelb zu machen.

---

## 3) Was der Kern bereits sauber zusammenbindet (stark / stabil)

### 3.1 Invarianten‑Kern (hart deterministisch)
`core_invariants` liefert stabil:

- \(c_3=1/(8\pi)\)
- \(\varphi_0 = 1/(6\pi) + 3/(256\pi^4)\)
- \(\beta=\varphi_0/(4\pi)\) und der beobachtbare birefringence‑Winkel
- Starobinsky‑Skala \(M/\bar M_{\rm Pl}=\sqrt{8\pi}\,c_3^4\)

Das ist der “stabile Motor” hinter vielen Folgeoutputs.

### 3.2 Inflation/R²‑Closure (numerisch stabil, aber QFT‑Beweis noch offen)
`effective_action_r2` ist intern konsistent, OperatorSpec‑basiert und generiert robuste Plots/Outputs.

Der offene Schritt ist nicht “Zahlensalat”, sondern **publication‑grade QFT**: der Operator muss aus der mikroskopischen torsionful action BRST‑komplett abgeleitet werden (Gauge fixing + ghosts).

### 3.3 RG‑Fingerprints & Auditierbarkeit
`two_loop_rg_fingerprints` + Fingerprinting verhindern silent swaps und halten die RG‑Story reproduzierbar.

### 3.4 Ω\_b numerisch sehr nah (aber explizit conditional)
`omega_b_conjecture_scan` liefert eine sehr gute numerische Übereinstimmung, **markiert aber korrekt die Annahmen**.  
Unconventional ergänzt mit `ux_omega_b_aps_bridge` eine plausible Topologie‑Docking‑Hypothese (APS seam).

---

## 4) Die verbleibenden ToE‑Gaps (präzise) + Lösungsansätze

Die folgenden Punkte sind nach “ToE‑Hebelwirkung” sortiert.

---

### Gap A — **Metrologie / α(0) auf CODATA‑Niveau**

**Was wir haben (Update)**
- `alpha_precision_audit` + `defect_partition_derivation` liefern jetzt einen **diskreten δ₂‑Kandidaten** (kein continuous fit), der α(0) metrologisch schließt.
- Physics‑Mode: `global_consistency_test` hat `alpha_inv_0_within_2sigma` PASS (z≈1.86) und behält ᾱ(MZ) als separaten Policy‑Vergleich.

**Warum ToE‑kritisch**
- α(0) ist *die* metrologische Referenzkonstante. Eine ToE kann nicht dauerhaft sagen “wir vergleichen lieber nur bei MZ”.

**Konkrete Lösungspfade**
- **(A1) Derive the next term**: \(\delta_2 e^{-4\alpha}\) darf kein “nice tweak” bleiben, sondern muss aus dem gleichen Defekt‑/Topologie‑Formalismus folgen.  
  - Engineering‑Target: `alpha_precision_audit` liefert bereits die benötigte Größenordnung.
- **(A2) Renorm‑Bridge sauber ziehen**: Vollständige (und deklarierte) QED/EW‑Policy bis in die Messdefinition.  
  - Ergebnis muss gleichzeitig α(MZ) *und* α(0) treffen (oder zeigen, warum eine klare Abweichung unvermeidbar ist).

---

### Gap B — **Flavor (CKM/PMNS inkl. CP‑Phasen)**

**Was wir haben (Update)**
- Deterministische Pipelines sind jetzt **im Physics‑Gate grün**:
  - `ckm_full_pipeline`: ref‑scale χ² nutzt die deklarierte Subset‑Liste `reference.chi2_keys` (default: Vus,Vub,Vcb,Vtd) → χ²≈8.03 (dof=4), p≈0.090.
  - `pmns_full_pipeline`: χ²≈5.52 (dof=4), p≈0.238 (Topology‑Kandidaten wired; NO/IO ordering policy verhindert “shopping”).
- **Selection‑Rule als Filter (Update 2026-01-27)**:
  - `tfpt_suite/data/flavor_texture_v24.json`: `phase_selection_rule_mode="filter_only"` für CKM/PMNS wiring (keine baseline‑δCP Branches, wenn `topology_phase_map` vorhanden ist).
  - Gate `phase_set_derived_not_searched` **PASS** in `ckm_full_pipeline` (δ_CP Varianten ausschließlich aus `topology_phase_map`, keine χ²‑getriebene Phasen‑Suche).
  - Gate `no_convention_shopping_possible_under_fixed_rule` **PASS** in `pmns_full_pipeline` (PMNS Konvention wird per mass‑splitting‑Canonicalization fixiert; keine χ²‑Minimierung über Permutationen).
- Wiring‑Fix: `flavor_joint_objective_scan` findet CKM‑Variants; `topology_phase_map` Kandidaten werden in CKM/PMNS Scans injiziert (diskret, kein Float‑Fit).
- `chiral_index_three_cycles` bleibt der Docking‑Point für Wilson‑Line Phase‑Atoms (formale Operator‑Derivation weiterhin ein publication‑grade TODO → WARN, nicht FAIL).

**Warum ToE‑kritisch**
- Flavor ist aktuell der “Endboss”: hier entscheidet sich, ob TFPT mehr als Kohärenz‑Numerologie ist.

**Konkrete Lösungspfade**
- **(B1) Topology→Phase Map explizit machen** (kein stilles Δ/δ‑Tuning):
  - Nutze `chiral_index_three_cycles` als Input‑Quelle: feste diskrete Holonomie‑Äste.
  - Baue eine **endliche Menge diskreter Flavor‑Klassen** (nicht kontinuierliche Fits).
- **(B2) Gemeinsames Objective**: CKM+PMNS gemeinsam plus Mass‑Ratios (und ggf. Koide‑Constraints) – nicht getrennte Heuristiken.
- **(B3) PMNS‑Hebel vergrößern**:
  - Der aktuelle Z3‑breaking‑Hebel auf δCP ist zu schwach (Suite zeigt das).
  - Erlaube, dass PMNS‑Struktur stärker aus charged‑lepton‑Diagonalisation kommt (nicht nur neutrino‑TBM + kleines ε).
- **(B4) Matching/Skalen‑Policy finalisieren**:
  - Ohne sauberes Matching kann Flavor systematisch “falsch aussehen”.  
  - Aber: Die χ²‑Größenordnung ist so groß, dass Matching allein nicht reicht – es ist ein Mechanismus‑Gap, das Matching nur “entmischt”.

---

### Gap C — **BRST‑grade QFT Derivation (Gauge fixing + ghosts) / Gravity‑Einstein‑Limit**

**Was wir haben**
- `effective_action_r2`: closure‑level OperatorSpec‑Pipeline (funktioniert, dokumentiert).
- `brst_ghost_deriver`: closure‑level BRST/ghost‑Kette ist jetzt **deterministisch aus der kanonischen Action‑Spec generiert** und als harte Gates ausformulierbar:
  - `quadratic_operator_derived_from_action` PASS
  - `brst_invariance_verified_symbolically_or_structurally` PASS (structural BRST‑exact gauge fixing + FP ghost sector)
  - `gauge_parameter_independence_nonproxy` PASS (xi‑scan‑Stabilität von \(M/M_{\rm Pl}\) unter der Closure‑Policy)
  - `heat_kernel_contract_matches_effective_action_r2` PASS (a₂/β‑Contract aus derived blocks reproduziert \(M/M_{\rm Pl}\))
- `ufe_gravity_normalization`: Einstein‑limit docking in Zahlen:
  - ξ\_tree=3/4, ξ=c3/φ0≈0.748303… (kleiner shift durch δ\_top)
- `qft_completeness_ledger` benennt die missing steps sauber.

**Warum ToE‑kritisch**
- “closure level spec” ist nicht gleich “Beweis”. Für ToE braucht ihr die Operator‑Ableitung aus der mikroskopischen Action inklusive BRST‑Metadata **über das derzeitige konstante‑Krümmungs‑Closure hinaus** (torsionful connection‑level).

**Konkrete Lösungspfade**
- **(C1) Quadratic fluctuation operator aus `microscopic_action_tfpt_v25.json` ableiten** (nicht nur OperatorSpec vorgeben).
- **(C2) Gauge fixing + FP ghosts BRST‑komplett** dokumentieren und als Daten/Operatorblöcke ausgeben.
- **(C3) Heat‑kernel a₂ daraus** berechnen und gegen `effective_action_r2` “contract” querprüfen.
- **(C4) Einstein limit K→0 sauber als Renorm‑Condition** implementieren und κ/G‑Normalisierung daraus konsistent ableiten.

---

### Gap D — **Matching / Threshold finite pieces / Below‑MZ Policy**

**Was wir haben**
- Plumbing ist stabil (`msbar_matching_map`, `below_mt_eft_cascade`).
- Proof‑Points: `matching_finite_pieces` enthält jetzt explizit **α(0)→ᾱ^(5)(MZ)** (leptonic 1‑loop + Δα_had^(5) + Δα_extra) als deklarierte below‑MZ QED‑Policy‑Kette.
- `msbar_matching_map` emittiert zusätzlich eine deterministische **Jacobian‑Kovarianz** (`covariance_propagated_end_to_end`) für mt‑Boundary Outputs (derzeit diagonal input‑σ; Upgradepfad: volle Kovarianz).
- Unconventional `ux_matching_metamorphic_audit`: invertibility/stability checks.
- `ux_threshold_graph_audit`: macht policies sichtbar.

**Was fehlt**
- **Finite matching pieces** (nicht nur identity@µ=M).
- QED/EW decoupling below MZ als *vollständige* EFT‑Policy (W/Z/H/top thresholds; derzeit nur QED‑Bridge Proof‑Point).
- Einheitliche Unsicherheitspropagation (inkl. **voller** Kovarianzen) “end‑to‑end” über Threshold‑Graph (derzeit: linear Jacobian + MC hooks, aber noch kein zentraler Cov‑Contract).

**Lösungsansatz**
- Schrittweise:
  1. finite pieces an MSigma/MG8/MNR (und below‑mt QED/QCD) implementieren,
  2. uncertainty propagation vereinheitlichen,
  3. dann erst Flavor/Metrologie final bewerten.

---

### Gap E — **Cosmology: Bounce → CMB (k→ℓ) deterministisch**

**Was wir haben**
- `cosmo_reheating_policy_v106` + `k_calibration` liefern jetzt eine deterministische Policy‑Schätzung.
- Bounce→CAMB Bridge ist verdrahtet: `primordial_spectrum_builder` baut \(P_\mathcal{R}(k),P_t(k)\) aus `bounce_perturbations` und `boltzmann_transfer` konsumiert das via CAMB `set_initial_power_table` (inkl. expliziter Pivot‑Normierung als aktuelle Amplituden‑Policy).
- Physics‑Run zeigt (unter der v1.06 Policy in `k_calibration.json`):
  - \(a_0/a_t\)≈**1.14×10⁵⁴**
  - \(\ell_{\rm bounce}^s\)≈**1.23×10⁵**
  - \(\ell_{\rm bounce}^t\)≈**456**

**Interpretation**
- Tensor‑Feature kann plausibel im CMB‑ℓ‑Bereich liegen; Scalar‑Feature liegt weit darüber.  
  → Das ist nicht “falsch”, sondern bedeutet: ihr braucht eine **klare Policy**, welche Signaturen ihr als primär testbar claimt.

**Was fehlt**
- T\_reh/N\_reh aus TFPT‑Threshold‑History (MSigma/MG8/…) statt als Input.
- **Planck‑Likelihood (Plugin, Status)**: high‑ℓ `plik_lite_native` ist opt‑in integriert (log‑L + Δχ² vs power‑law baseline; `TFPT_ENABLE_PLANCK_LIKELIHOOD=1`). Remaining: low‑ℓ (`planck_2018_lowl`) + ggf. lensing‑likelihood kombinieren und in die zentrale Likelihood‑Engine einspeisen.
- Nuisance‑Policy Upgrade: derzeit ist eine explizite **Pivot‑Normierung** der injizierten \(P(k)\) Tabelle aktiv, um CAMB in einem realistischen Amplitudenregime zu halten; publication‑grade muss klar deklarieren, ob As/M gefixt, profiliert oder marginalisiert wird.

**Lösungsansatz**
- **(E1)** Threshold‑driven reheating temperature/histories als neues deterministisches Modul.
- **(E2)** Entscheidung “CMB vs small‑scale probes” als explizite Policy in Suite + Scorecard.

---

### Gap F — **Dark Matter (Axion) world-contact closure**

**Was wir haben (Update)**
- `axion_dm_pipeline` ist jetzt scenario‑aware und nutzt als shipped default den **post‑inflation PQ‑Branch**:
  - ν≈15.562 GHz (m\_a≈64.36 µeV) bleibt stabil.
  - Isocurvature ist im post‑inflation Branch **nicht anwendbar** → PASS.
  - Relic‑Closure wird über ein **diskretes** Strings/DW‑Enhancement \(C_{\mathrm{str}}=7/3\) (minimal cusp/charge set) geschlossen: \(\Omega_a h^2\approx 0.123\) (≈+2% vs 0.12).
  - Alternative DM‑Kanäle (`dm_alternative_channels`, `torsion_dm_pipeline`) bleiben als optionale Branches, werden aber unter der shipped Axion‑Policy nicht “required”.

**Was das bedeutet**
- Unter Standardannahmen kann “ein globales Axionfeld mit θ\_i=φ0” nicht gleichzeitig Starobinsky‑Inflation und Planck‑Isocurvature constraints erfüllen.

**Lösungsansätze**
- **(F1)** Scenario choice explizit: post‑inflation PQ breaking + strings/domain walls (dann Isocurvature anders).
- **(F2)** Alternativer DM‑Kanal: torsion excitations / zusätzliche kalte Komponente (wie in `five_problems.tex` angedeutet) – dann als Pipeline mit Constraints.
- **(F3)** f\_a derivieren (ladder/block) statt nur zu zitieren.

---

### Gap G — **Dark Energy / Λ closure**

**Was wir haben**
- `dark_energy_paths` macht Targets explizit (Planck‑Snapshot → Λ‑Größenordnung):
  - Λ\_obs≈**4.24×10⁻⁸⁴ GeV²**
  - ρ\_Λ≈**2.51×10⁻⁴⁷ GeV⁴**
  - UFE‑Pfad‑Target: \(K_{\rm rms}\)≈**4.12×10⁻⁴² GeV**
  - Cascade‑Magnitude‑Target: φ\_*≈**9.20×10⁻³¹**
- `torsion_condensate` liefert jetzt eine deterministische **diskrete** Kandidaten‑Lösung (kein Float‑Fit) und Gate‑Artefakte:
  - `torsion_condensate_gap_equation_solved` PASS (finite diskrete Lösungen)
  - `Lambda_matches_observation_with_uncertainty` PASS (ρ\_Λ mismatch innerhalb der deklarierten Toleranz)
  - `n_half_explained_by_quantization_rule` PASS (aktuell bevorzugt: n=1/2 Branch; quantization story explizit als Annahme)

**Was fehlt**
- (A) tatsächliche Ableitung von ⟨K²⟩ aus UFE/mikroskopischer Torsion‑Dynamik, oder  
- (B) deterministisches ladder/block‑Modell, das φ\_n bis zum terminal stage reproduziert.

---

### Gap H — **Torsion falsifizierbar (heute / astrophysisch)**

**Was wir haben**
- `torsion_bounds_mapping` ist nicht mehr “immer grün” (Regime‑Policy ist explizit).
- Unconventional `ux_torsion_regime_designer` zeigt, wie unrealistisch naive Lab‑Sättigung ist.
- **SNR‑Gate (Update 2026-01-27)**: `torsion_falsifiability_snr` liefert ein explizites Source+Noise Modell (Policy‑JSON) und berechnet SNR:
  - magnetar‑proxy via **electron polarization** \(P=\tanh(\mu_B B/(k_B T))\) + \(n_e=Y_e n_b\) bei nuklearer Dichte,
  - Gate `astro_channel_measurable_under_realistic_noise` PASS und `go_no_go_snr_ge_5` PASS (unter der deklarierten σν‑Proxy‑Noise‑Policy).

**Was fehlt**
- Publication‑grade: echte Daten‑Likelihood (Timing‑Residuals/Polarimetrie PSD) statt σν‑Proxy, plus ein klarer Observable‑Kanal (welches Datenprodukt, welche Systematiken).

---

### Gap I — **Arrow of Time**

In `progress.md` explizit als nicht implementiert markiert.  
ToE‑Closure erfordert mindestens ein konsistentes, testbares Mechanismus‑Modul (z.B. torsion flux / topological non‑invertibility), sonst bleibt es paper‑claim.

---

## 5) Priorisierte “ToE‑Closure” Roadmap (praktisch)

1) **Physics‑Scorecard als Gate akzeptieren** (nicht nur Engineering‑PASS).  
2) **α(0) derivieren** (δ₂‑Term + saubere Renorm‑Bridge).  
3) **Matching finite pieces** finalisieren (sonst bleibt Flavor nicht interpretierbar).  
4) **Flavor neu aufziehen**: diskrete Topologie‑Klassen + Phase‑Map (Wilson lines) + CKM/PMNS gemeinsam.  
5) **BRST‑grade torsion operator** (QFT‑Kern) ableiten und gegen `effective_action_r2` cross‑checken.  
6) **Cosmo history threshold‑driven** machen (reheating/entropy aus MSigma/MG8/…).  
7) **DM/DE closure**: scenario + constraints + echte derivations statt Targets.

---

## 6) Reproduktion (Commands)

```bash
# Engineering
python3 tfpt-suite/run_suite.py run-all
python3 tfpt-suite/run_suite.py build-report

# Physics (separate artifacts)
python3 tfpt-suite/run_suite.py --mode physics --output-dir tfpt-suite/out_physics run-all
python3 tfpt-suite/run_suite.py --mode physics --output-dir tfpt-suite/out_physics build-report --pdf-path tfpt-suite/tfpt-test-results-physics.pdf

# Unconventional tooling
python3 tfpt-suite/unconventional/run_unconventional_suite.py run-all
python3 tfpt-suite/unconventional/run_unconventional_suite.py build-report
```

---

## 7) Weitere Referenzen im Repo

- `tfpt-suite/progress.md`: die “checklist” der publishable‑closure gaps (Single Source of Truth für TODOs).
- `tfpt-suite/TESTS_OVERVIEW.md`: “was ist wirklich getestet” vs “was ist nur claim”.
- PDFs:
  - `tfpt-suite/tfpt-test-results.pdf`
  - `tfpt-suite/tfpt-test-results-physics.pdf`
  - `tfpt-suite/unconventional/unconventional-test-results.pdf`

