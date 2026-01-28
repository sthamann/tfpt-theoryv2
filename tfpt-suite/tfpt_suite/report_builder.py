from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# =============================================================================
# MODULE GROUPING BY THEORETICAL QUESTION
# =============================================================================

# Conventional suite groups: modules organized by the theoretical question they answer
CONVENTIONAL_MODULE_GROUPS: list[tuple[str, str, list[str]]] = [
    (
        "I. Fundamental Invariants",
        "What algebraic constants does TFPT define and are they consistent?",
        [
            "core_invariants",
            "discrete_consistency_uniqueness",
            "discrete_complexity_minimizer",
        ],
    ),
    (
        "II. Fine-Structure Constant α",
        "Can TFPT predict the fine-structure constant from first principles?",
        [
            "alpha_precision_audit",
            "alpha_on_shell_bridge",
        ],
    ),
    (
        "III. Renormalization Group Running",
        "How do gauge couplings run and where do they cross TFPT invariants?",
        [
            "two_loop_rg_fingerprints",
            "unification_gate",
            "msbar_matching_map",
            "below_mt_eft_cascade",
            "stability_unitarity_audit",
        ],
    ),
    (
        "IV. QFT Consistency",
        "Is the theory anomaly-free and QFT-complete?",
        [
            "anomaly_cancellation_audit",
            "qft_completeness_ledger",
        ],
    ),
    (
        "V. Gravity & R² Sector",
        "Does TFPT derive the Starobinsky R² coefficient and Planck scale?",
        [
            "effective_action_r2",
            "ufe_gravity_normalization",
            "aps_eta_gluing",
        ],
    ),
    (
        "VI. Cosmology & Bounce",
        "How do perturbations propagate through the bounce and map to CMB multipoles?",
        [
            "bounce_perturbations",
            "k_calibration",
            "cosmo_reheating_policy_v106",
        ],
    ),
    (
        "VII. Flavor: CKM Matrix",
        "Can the Z₃ Möbius texture reproduce the CKM matrix?",
        [
            "mobius_cusp_classification",
            "mobius_delta_calibration",
            "mobius_z3_yukawa_generator",
            "ckm_full_pipeline",
        ],
    ),
    (
        "VIII. Flavor: PMNS Matrix & Neutrinos",
        "Does TFPT correctly predict neutrino mixing angles and masses?",
        [
            "pmns_z3_breaking",
            "pmns_mechanism_bridge",
            "pmns_full_pipeline",
            "seesaw_block",
        ],
    ),
    (
        "IX. Torsion & Birefringence",
        "Is TFPT torsion compatible with experimental bounds and testable via birefringence?",
        [
            "torsion_bounds_mapping",
            "birefringence_tomography",
        ],
    ),
    (
        "X. Baryons & Dark Matter",
        "Can TFPT explain the baryonic density Ωb and axion DM?",
        [
            "omega_b_conjecture_scan",
            "axion_dm_pipeline",
            "axion_scenario_matrix",
            "baryogenesis_placeholder",
        ],
    ),
    (
        "XI. Dark Energy",
        "Which mechanisms for Λ are TFPT-compatible?",
        [
            "dark_energy_paths",
        ],
    ),
    (
        "XII. Global Consistency & Dashboard",
        "How well does TFPT overall fit the experimental data?",
        [
            "predictions_dashboard",
            "global_consistency_test",
        ],
    ),
    (
        "XIII. Topology & Chirality",
        "Does TFPT explain three families via the chiral index?",
        [
            "chiral_index_three_cycles",
            "defect_partition_derivation",
        ],
    ),
    (
        "XIV. Additional Checks",
        "Additional consistency checks (masses, BBN, GW)",
        [
            "mass_spectrum_minimal",
            "bbn_neff_sanity",
            "gw_background_bounds",
            "g2_and_lamb_shift_proxy",
        ],
    ),
]

# Unconventional suite groups
UNCONVENTIONAL_MODULE_GROUPS: list[tuple[str, str, list[str]]] = [
    (
        "A. Matching & Threshold Robustness",
        "Are the matching primitives and threshold transitions correctly implemented?",
        [
            "ux_matching_metamorphic_audit",
            "ux_threshold_graph_audit",
        ],
    ),
    (
        "B. k→ℓ Bridge Feasibility",
        "Can the bounce reach CMB multipoles under plausible priors?",
        [
            "ux_cosmo_history_sampler",
        ],
    ),
    (
        "C. Flavor Convention Search",
        "Which discrete flavor conventions perform best (with holdout)?",
        [
            "ux_flavor_holdout_search",
        ],
    ),
    (
        "D. Gravity Gauge-Fixing",
        "Which gauge-fixing choice minimizes non-minimal operator structure?",
        [
            "ux_gravity_gaugefix_ga",
        ],
    ),
    (
        "E. Ωb-APS Bridge",
        "Can the APS seam term explain the Ωb coefficient?",
        [
            "ux_omega_b_aps_bridge",
        ],
    ),
    (
        "F. Torsion Regime Design",
        "Which physical regimes could make TFPT torsion falsifiable?",
        [
            "ux_torsion_regime_designer",
        ],
    ),
]

# Modules that use RG running and should display PyR@TE config info
RG_MODULES = {
    "two_loop_rg_fingerprints",
    "ckm_full_pipeline",
    "pmns_full_pipeline",
    "pmns_mechanism_bridge",
    "stability_unitarity_audit",
    "unification_gate",
    "ux_threshold_graph_audit",
    "ux_flavor_holdout_search",
}


def _workspace_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _relpath_display(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _sanitize_paths(text: str) -> str:
    root = str(_workspace_root().resolve())
    if not root.endswith("/"):
        root_slash = root + "/"
    else:
        root_slash = root
        root = root.rstrip("/")
    return text.replace(root_slash, "").replace(root, "")


def _zwsp_break(text: str) -> str:
    """Previously added zero-width spaces for line breaking, but these render as black boxes
    in PDFs with WinAnsiEncoding. Now disabled - just returns the text unchanged."""
    # Zero-width space (U+200B) is not supported by WinAnsiEncoding fonts
    # Simply return text unchanged to avoid black boxes
    return text


def _truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    if value is None:
        return False
    return bool(value)


def _wrap_lines(text: str, max_len: int) -> str:
    out: list[str] = []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if len(line) <= max_len:
            out.append(line)
            continue
        start = 0
        while start < len(line):
            out.append(line[start : start + max_len])
            start += max_len
    return "\n".join(out)


def _replace_unicode_symbols(text: str) -> str:
    """Replace Unicode symbols and problematic characters that may not render correctly in PDF with ASCII equivalents.
    
    This function uses a comprehensive replacement map for all Unicode characters found in TFPT data files.
    The replacements are designed to be readable in a scientific/physics context.
    """
    # Comprehensive replacement map for all Unicode characters found in TFPT data
    # Organized by category for maintainability
    replacements = {
        # Greek lowercase letters
        "α": "alpha",
        "β": "beta",
        "γ": "gamma",
        "δ": "delta",
        "ε": "epsilon",
        "ζ": "zeta",
        "η": "eta",
        "θ": "theta",
        "κ": "kappa",
        "λ": "lambda",
        "μ": "mu",
        "ν": "nu",
        "ξ": "xi",
        "π": "pi",
        "ρ": "rho",
        "σ": "sigma",
        "τ": "tau",
        "φ": "phi",
        "χ": "chi",
        "ψ": "psi",
        "ω": "omega",
        # Greek uppercase letters
        "Γ": "Gamma",
        "Δ": "Delta",
        "Λ": "Lambda",
        "Σ": "Sigma",
        "Ω": "Omega",
        # Mathematical operators and relations
        "≈": "~",
        "≃": "~",
        "≡": "===",
        "≠": "!=",
        "≤": "<=",
        "≥": ">=",
        "≲": "<=~",
        "∝": "~",
        "±": "+/-",
        "−": "-",  # Unicode minus sign
        "×": "*",
        "·": "*",
        "∞": "inf",
        # Arrows
        "→": "->",
        "←": "<-",
        "↔": "<->",
        "⇒": "=>",
        "⇐": "<=",
        # Calculus and set theory
        "∫": "Int",
        "∮": "Int",
        "∂": "d",
        "∇": "nabla",
        "√": "sqrt",
        "∈": "in",
        "⊃": "supset",
        # Superscripts
        "⁻": "^-",
        "⁺": "^+",
        "⁰": "^0",
        "¹": "^1",
        "²": "^2",
        "³": "^3",
        "⁴": "^4",
        "⁵": "^5",
        "⁶": "^6",
        "⁷": "^7",
        "⁸": "^8",
        "⁹": "^9",
        # Subscripts
        "₀": "_0",
        "₁": "_1",
        "₂": "_2",
        "₃": "_3",
        "₄": "_4",
        "₅": "_5",
        "₆": "_6",
        "₇": "_7",
        "₈": "_8",
        "₉": "_9",
        # Special mathematical symbols
        "ℓ": "l",  # script l
        "′": "'",  # prime
        "†": "+",  # dagger (used for adjoint)
        "⋆": "*",  # star operator
        "□": "[]",  # d'Alembertian
        "°": "deg",
        # Brackets
        "⟨": "<",
        "⟩": ">",
        # Combining characters (often problematic)
        "̄": "",  # combining macron - remove
        "̂": "",  # combining circumflex - remove
        # Dashes and punctuation
        "–": "-",  # en dash
        "—": "--",  # em dash
        "'": "'",  # right single quote
        "'": "'",  # left single quote
        """: '"',  # left double quote
        """: '"',  # right double quote
        # Accented characters (common in names/references)
        "ö": "o",
        "ä": "a",
        "ü": "u",
        "é": "e",
        "è": "e",
        "ê": "e",
        "à": "a",
        "â": "a",
        "ç": "c",
        "ñ": "n",
        "ý": "y",
        "ĝ": "g",
        # Micro sign (distinct from Greek mu)
        "µ": "mu",
    }
    
    result = text
    for unicode_char, ascii_replacement in replacements.items():
        result = result.replace(unicode_char, ascii_replacement)
    
    # Final pass: replace any remaining non-ASCII characters with '?'
    # This ensures we never have encoding issues in the PDF
    final = []
    for char in result:
        if ord(char) <= 127:
            final.append(char)
        else:
            # Unknown non-ASCII character - use '?' as fallback
            final.append('?')
    
    return ''.join(final)


def _format_mpf_string(text: str) -> str:
    """Format mpf('...') strings to be more readable."""
    import re
    
    # Replace mpf('...') with just the number, formatted nicely
    def format_mpf(match):
        num_str = match.group(1)
        try:
            # Try to parse as float
            num = float(num_str)
            # Format with reasonable precision
            if abs(num) < 1e-3 or abs(num) > 1e6:
                return f"{num:.6e}"
            else:
                return f"{num:.10f}".rstrip("0").rstrip(".")
        except ValueError:
            return num_str
    
    # Replace mpf('...') patterns
    text = re.sub(r"mpf\(['\"]([^'\"]+)['\"]\)", format_mpf, text)
    
    # Also handle mpf without quotes
    text = re.sub(r"mpf\(([^)]+)\)", format_mpf, text)
    
    return text


def _format_numerical_data(text: str) -> str:
    """Format numerical data in text for better readability."""
    import re
    
    # First format mpf strings
    text = _format_mpf_string(text)
    
    # Format very long decimal numbers (more than 15 digits after decimal)
    def format_long_decimal(match):
        num_str = match.group(0)
        try:
            num = float(num_str)
            if abs(num) < 1e-3 or abs(num) > 1e6:
                return f"{num:.6e}"
            else:
                # Round to 12 significant digits
                return f"{num:.12g}"
        except ValueError:
            return num_str
    
    # Match numbers with many decimal places
    text = re.sub(r"\d+\.\d{15,}", format_long_decimal, text)
    
    return text


def _format_report_text(text: str) -> str:
    """Format report.txt content for better PDF readability."""
    # Replace Unicode symbols FIRST - this is critical for PDF rendering
    text = _replace_unicode_symbols(text)
    
    # Format numerical data
    text = _format_numerical_data(text)
    
    # Apply Unicode replacement again after numerical formatting (in case formatting introduced new issues)
    text = _replace_unicode_symbols(text)
    
    # Break up very long lines (e.g., long dictionary outputs)
    lines = text.split("\n")
    formatted_lines = []
    for line in lines:
        # If line is very long and contains key=value pairs, try to break it
        if len(line) > 150 and "=" in line:
            # Try to break at commas or spaces after key=value pairs
            parts = []
            current = ""
            for char in line:
                current += char
                if len(current) > 100 and (char == "," or char == " "):
                    parts.append(current.rstrip())
                    current = ""
            if current:
                parts.append(current)
            if len(parts) > 1:
                formatted_lines.extend(parts)
            else:
                formatted_lines.append(line)
        else:
            formatted_lines.append(line)
    
    return "\n".join(formatted_lines)


def _read_text_if_exists(path: Path) -> str | None:
    try:
        if path.is_file():
            return path.read_text(encoding="utf-8")
    except Exception:
        return None
    return None


def _normalize_heading(text: str) -> str:
    return " ".join([part for part in "".join([c.lower() if c.isalnum() else " " for c in str(text)]).split() if part]).strip()


def _strip_blank_lines(lines: list[str]) -> list[str]:
    i = 0
    j = len(lines)
    while i < j and not lines[i].strip():
        i += 1
    while j > i and not lines[j - 1].strip():
        j -= 1
    return lines[i:j]


def _parse_markdown_h2_sections(md: str) -> list[tuple[str, str]]:
    sections: list[tuple[str, str]] = []
    cur_title: str | None = None
    cur_lines: list[str] = []

    for raw in (md or "").splitlines():
        line = raw.rstrip("\n")
        if line.startswith("## "):
            if cur_title is not None:
                body = "\n".join(_strip_blank_lines(cur_lines))
                sections.append((cur_title, body))
            cur_title = line[3:].strip()
            cur_lines = []
            continue
        if line.startswith("# "):
            continue
        if cur_title is None:
            continue
        cur_lines.append(line)

    if cur_title is not None:
        body = "\n".join(_strip_blank_lines(cur_lines))
        sections.append((cur_title, body))
    return sections


def _extract_doc_sections(md: str, *, wanted_prefixes: list[str]) -> list[tuple[str, str]]:
    wanted = [_normalize_heading(x) for x in wanted_prefixes]
    out: list[tuple[str, str]] = []
    for title, body in _parse_markdown_h2_sections(md):
        h = _normalize_heading(title)
        for pref in wanted:
            if h.startswith(pref):
                out.append((title, body))
                break
    return out


def _extract_pyrate_info(results_payload: dict[str, Any]) -> dict[str, str] | None:
    """
    Extract PyR@TE model configuration info from results.json if available.
    Returns dict with model_name, yaml_source, yaml_sha256, pythonoutput_dir, etc.
    """
    # Check various locations where PyR@TE info might be stored
    info: dict[str, str] = {}

    # Check model_fingerprint (used by two_loop_rg_fingerprints)
    fingerprint = results_payload.get("results", {}).get("model_fingerprint", {})
    if fingerprint:
        if fingerprint.get("model_name_expected"):
            info["model_name"] = fingerprint["model_name_expected"]
        if fingerprint.get("yaml_source"):
            info["yaml_source"] = fingerprint["yaml_source"]
        if fingerprint.get("yaml_source_sha256"):
            info["yaml_sha256"] = fingerprint["yaml_source_sha256"][:16] + "..."
        if fingerprint.get("pythonoutput_module_file"):
            info["pythonoutput_file"] = fingerprint["pythonoutput_module_file"]

    # Check generation info
    gen = results_payload.get("results", {}).get("generation", {})
    if gen:
        if gen.get("model_module") and "model_name" not in info:
            info["model_name"] = gen["model_module"]
        if gen.get("pythonoutput_dir") and "pythonoutput_file" not in info:
            info["pythonoutput_dir"] = gen["pythonoutput_dir"]
        config = gen.get("config", {})
        if config.get("yaml_source") and "yaml_source" not in info:
            info["yaml_source"] = config["yaml_source"]
        if config.get("yaml_source_sha256") and "yaml_sha256" not in info:
            info["yaml_sha256"] = config["yaml_source_sha256"][:16] + "..."

    # Check rg_config (used by some modules)
    rg_config = results_payload.get("results", {}).get("rg_config", {})
    if rg_config:
        if rg_config.get("model_module") and "model_name" not in info:
            info["model_name"] = rg_config["model_module"]
        if rg_config.get("yaml_source") and "yaml_source" not in info:
            info["yaml_source"] = rg_config["yaml_source"]

    # Check pyrate_config
    pyrate_config = results_payload.get("results", {}).get("pyrate_config", {})
    if pyrate_config:
        if pyrate_config.get("model") and "model_name" not in info:
            info["model_name"] = pyrate_config["model"]
        if pyrate_config.get("yaml_file") and "yaml_source" not in info:
            info["yaml_source"] = pyrate_config["yaml_file"]

    return info if info else None


@dataclass(frozen=True)
class ModuleArtifacts:
    module_id: str
    title: str
    out_dir: Path
    meta: dict[str, Any]
    results_payload: dict[str, Any]
    report_text: str
    plots: list[Path]

    @property
    def checks(self) -> list[dict[str, Any]]:
        return list(self.results_payload.get("checks", []))

    @property
    def warnings(self) -> list[str]:
        return list(self.results_payload.get("warnings", []))

    @property
    def spec(self) -> dict[str, Any]:
        return dict(self.results_payload.get("spec", {}))

    def checks_summary(self) -> tuple[int, int, int, int]:
        total = len(self.checks)
        passed_strict = 0
        warned = 0
        failed = 0
        for c in self.checks:
            sev = str(c.get("severity") or "").strip().upper()
            if not sev:
                sev = "PASS" if _truthy(c.get("passed")) else "FAIL"

            if sev == "PASS":
                passed_strict += 1
            elif sev == "WARN":
                warned += 1
            elif sev == "FAIL":
                failed += 1
            else:
                passed_strict += 1
        return passed_strict, total, warned, failed

    def get_pyrate_info(self) -> dict[str, str] | None:
        return _extract_pyrate_info(self.results_payload)


def discover_module_artifacts(out_dir: Path) -> list[ModuleArtifacts]:
    if not out_dir.exists():
        raise FileNotFoundError(f"Output directory not found: {out_dir}")

    modules: list[ModuleArtifacts] = []
    for child in sorted(out_dir.iterdir()):
        if not child.is_dir():
            continue
        meta_path = child / "meta.json"
        results_path = child / "results.json"
        report_path = child / "report.txt"
        if not results_path.exists():
            continue

        meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        report_text = report_path.read_text(encoding="utf-8") if report_path.exists() else ""

        # Find plot files - check both files and ensure they're actually files, not directories
        plots = sorted([
            p for p in child.iterdir() 
            if p.is_file() and p.suffix.lower() in {".png", ".jpg", ".jpeg"}
        ])

        module_id = str(meta.get("module", {}).get("id") or child.name)
        title = str(meta.get("module", {}).get("title") or payload.get("spec", {}).get("name") or module_id)

        modules.append(
            ModuleArtifacts(
                module_id=module_id,
                title=title,
                out_dir=child,
                meta=meta,
                results_payload=payload,
                report_text=report_text,
                plots=plots,
            )
        )

    return modules


def _group_modules(
    mods: list[ModuleArtifacts],
    groups: list[tuple[str, str, list[str]]],
) -> list[tuple[str, str, list[ModuleArtifacts]]]:
    """
    Group modules by their theoretical question.
    Returns list of (group_title, group_question, [modules]).
    Modules not in any group go to "Other".
    """
    mod_by_id = {m.module_id: m for m in mods}
    used_ids: set[str] = set()
    result: list[tuple[str, str, list[ModuleArtifacts]]] = []

    for group_title, group_question, module_ids in groups:
        group_mods: list[ModuleArtifacts] = []
        for mid in module_ids:
            if mid in mod_by_id:
                group_mods.append(mod_by_id[mid])
                used_ids.add(mid)
        if group_mods:
            result.append((group_title, group_question, group_mods))

    # Remaining modules
    remaining = [m for m in mods if m.module_id not in used_ids]
    if remaining:
        result.append(("XV. Other Modules", "Additional modules without specific grouping", remaining))

    return result


def build_tfpt_test_results_pdf(*, out_dir: Path, pdf_path: Path, verification_mode: str | None = None) -> Path:
    """
    Build a grouped, nicely formatted PDF report.
    
    Args:
        out_dir: Directory containing module outputs
        pdf_path: Output PDF path
        verification_mode: Optional verification mode ("engineering" or "physics")
    """
    from datetime import datetime

    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import (
        Image,
        KeepTogether,
        PageBreak,
        Paragraph,
        Preformatted,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    mods = discover_module_artifacts(out_dir)
    if not mods:
        raise RuntimeError(f"No module outputs found under: {out_dir}")

    pdf_path = pdf_path.resolve()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine if unconventional
    out_rel = _relpath_display(out_dir).replace("\\", "/")
    is_unconventional = ("out/unconventional" in out_rel) or out_dir.name == "unconventional"
    
    # Determine verification mode (default to engineering if not specified)
    if verification_mode is None:
        # Try to infer from pdf_path or out_dir
        pdf_str = str(pdf_path).lower()
        out_str = str(out_dir).lower()
        if "physics" in pdf_str or "physics" in out_str or "out_physics" in out_str:
            verification_mode = "physics"
        else:
            verification_mode = "engineering"
    is_physics_mode = verification_mode == "physics"

    # Select appropriate grouping
    groups = UNCONVENTIONAL_MODULE_GROUPS if is_unconventional else CONVENTIONAL_MODULE_GROUPS
    grouped_modules = _group_modules(mods, groups)

    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.8 * cm,
        title="TFPT Test Results",
    )

    styles = getSampleStyleSheet()

    # Custom styles with better typography
    styles.add(ParagraphStyle(
        name="TitleCustom",
        parent=styles["Title"],
        fontSize=22,
        spaceAfter=12,
        textColor=colors.HexColor("#1a365d"),
    ))
    styles.add(ParagraphStyle(
        name="GroupTitle",
        parent=styles["Heading1"],
        fontSize=16,
        spaceBefore=18,
        spaceAfter=8,
        textColor=colors.HexColor("#2c5282"),
        borderWidth=1,
        borderColor=colors.HexColor("#bee3f8"),
        borderPadding=6,
        backColor=colors.HexColor("#ebf8ff"),
    ))
    styles.add(ParagraphStyle(
        name="GroupQuestion",
        parent=styles["BodyText"],
        fontSize=10,
        leading=14,
        spaceAfter=12,
        textColor=colors.HexColor("#4a5568"),
        leftIndent=10,
        fontName="Helvetica-Oblique",
    ))
    styles.add(ParagraphStyle(
        name="ModuleTitle",
        parent=styles["Heading2"],
        fontSize=12,
        spaceBefore=14,
        spaceAfter=6,
        textColor=colors.HexColor("#2d3748"),
    ))
    styles.add(ParagraphStyle(
        name="SectionHead",
        parent=styles["Heading3"],
        fontSize=10,
        spaceBefore=8,
        spaceAfter=4,
        textColor=colors.HexColor("#4a5568"),
    ))
    styles.add(ParagraphStyle(
        name="BodyCustom",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
    ))
    styles.add(ParagraphStyle(
        name="PyRateInfo",
        parent=styles["BodyText"],
        fontSize=8,
        leading=10,
        textColor=colors.HexColor("#553c9a"),
        backColor=colors.HexColor("#faf5ff"),
        borderWidth=1,
        borderColor=colors.HexColor("#d6bcfa"),
        borderPadding=4,
        leftIndent=5,
        rightIndent=5,
    ))
    styles.add(ParagraphStyle(
        name="CodeSmall",
        parent=styles["BodyText"],
        fontName="Courier",
        fontSize=7,
        leading=8.5,
    ))

    available_w = float(A4[0] - doc.leftMargin - doc.rightMargin)

    # Table cell styles
    cell = ParagraphStyle(name="TblCell", parent=styles["BodyText"], fontSize=7, leading=9, wordWrap="CJK")
    cell_bold = ParagraphStyle(name="TblCellBold", parent=cell, fontName="Helvetica-Bold")
    cell_code = ParagraphStyle(name="TblCellCode", parent=cell, fontName="Courier")

    def P(txt: str, *, style: ParagraphStyle = cell) -> Paragraph:
        from xml.sax.saxutils import escape
        # Replace Unicode symbols BEFORE escaping
        txt = _replace_unicode_symbols(str(txt))
        return Paragraph(escape(txt), style)

    story: list[Any] = []

    # === TITLE PAGE ===
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if is_unconventional:
        title = "TFPT Unconventional Test Results"
    elif is_physics_mode:
        title = "TFPT Test Results (Physics Mode)"
    else:
        title = "TFPT Test Results (Engineering Mode)"
    
    story.append(Paragraph(title, styles["TitleCustom"]))
    story.append(Spacer(1, 0.3 * cm))

    module_count = len(mods)

    if is_unconventional:
        intro_text = """
        <b>Unconventional Suite</b>: Search/Audit/Conjecture tooling for TFPT.<br/><br/>
        
        This PDF contains the results of the <b>Unconventional Suite</b>, which consists of {module_count} modules.
        These modules are <b>not</b> publication-grade physics derivations, but rather 
        tools to accelerate closing ToE gaps (GA search, enumeration, metamorphic tests).<br/><br/>
        
        <b>Content:</b> Each module is grouped by its theoretical question and 
        provides search/audit/design tools for specific TFPT aspects. The results serve 
        exploration and finding solution paths, not final verification.<br/><br/>
        
        <b>Difference from other variants:</b> In contrast to the Conventional Suite (Engineering/Physics), 
        this suite focuses on experimental search procedures and design tools, not on 
        deterministic verification.
        """.format(module_count=module_count)
    elif is_physics_mode:
        intro_text = """
        <b>Conventional Suite - Physics Mode</b>: Strict physics validation of TFPT theory.<br/><br/>
        
        This PDF contains the results of the <b>Conventional Suite in Physics Mode</b> with all 
        {module_count} modules. Physics Mode uses <b>stricter interpretation criteria</b> than 
        Engineering Mode.<br/><br/>
        
        <b>Content:</b> All modules implement deterministic, reproducible checks for the 
        central TFPT predictions and consistency tests. Modules are grouped by theoretical questions 
        to clarify the connection between tests and theory.<br/><br/>
        
        <b>Differences from Engineering Mode:</b><br/>
        • <b>Stricter Checks:</b> Large deviations from experimental references are upgraded to 
          WARN/FAIL (in Engineering Mode only numerical consistency)<br/>
        • <b>Extended Scorecards:</b> Additional metrics for alpha(0) metrology and flavor chi^2 imports<br/>
        • <b>Missing-Derivation Flags:</b> Missing derivations are explicitly marked as warnings<br/>
        • <b>Usage:</b> Pre-submission audit, identification of real physics gaps<br/><br/>
        
        <b>Difference from Unconventional Suite:</b> This suite contains publication-grade 
        verification modules, while the Unconventional Suite contains search/design tools.
        """.format(module_count=module_count)
    else:
        intro_text = """
        <b>Conventional Suite - Engineering Mode</b>: Deterministic verification of TFPT theory.<br/><br/>
        
        This PDF contains the results of the <b>Conventional Suite in Engineering Mode</b> with all 
        {module_count} modules. Engineering Mode focuses on deterministic execution and explicit 
        assumptions.<br/><br/>
        
        <b>Content:</b> All modules implement deterministic, reproducible checks for the 
        central TFPT predictions and consistency tests. Modules are grouped by theoretical questions 
        to clarify the connection between tests and theory.<br/><br/>
        
        <b>Check Semantics:</b> PASS/FAIL is based on numerical stability and consistency. 
        Deviations from experimental references are documented but not automatically 
        treated as errors.<br/><br/>
        
        <b>Differences from Physics Mode:</b><br/>
        • <b>Focus:</b> Numerical consistency and deterministic execution (not strict 
          physics validation)<br/>
        • <b>Check Interpretation:</b> Algorithmic/numerical success is the focus<br/>
        • <b>Usage:</b> Development, debugging, CI validation<br/><br/>
        
        <b>Difference from Unconventional Suite:</b> This suite contains publication-grade 
        verification modules, while the Unconventional Suite contains search/design tools.
        """.format(module_count=module_count)

    story.append(Paragraph(intro_text, styles["BodyCustom"]))
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph(f"<b>Generated:</b> {now}", styles["BodyCustom"]))
    story.append(Paragraph(f"<b>Output directory:</b> <font name='Courier'>{_relpath_display(out_dir)}</font>", styles["BodyCustom"]))
    story.append(Spacer(1, 0.6 * cm))

    # === OVERVIEW TABLE ===
    story.append(Paragraph("Module Overview", styles["GroupTitle"]))

    overview_rows: list[list[Any]] = [[
        P("Group", style=cell_bold),
        P("Module", style=cell_bold),
        P("Title", style=cell_bold),
        P("Checks", style=cell_bold),
        P("Status", style=cell_bold),
    ]]

    for group_title, _, group_mods in grouped_modules:
        # Use full group name (e.g., "I. Fundamental Invariants" instead of just "I.")
        group_display = group_title
        for m in group_mods:
            passed_strict, total, warned, failed = m.checks_summary()
            status = f"{passed_strict}/{total}" if total else "—"

            # Color-coded status
            if failed > 0:
                status_str = f"FAIL:{failed}"
            elif warned > 0:
                status_str = f"WARN:{warned}"
            else:
                status_str = "OK"

            overview_rows.append([
                P(_replace_unicode_symbols(group_display)),
                P(_replace_unicode_symbols(m.module_id)),
                P(_replace_unicode_symbols(m.title)),
                P(status),
                P(status_str),
            ])

    col_fracs = [0.25, 0.20, 0.35, 0.10, 0.10]  # Increased Group column width for full names
    col_widths = [available_w * f for f in col_fracs]
    tbl = Table(overview_rows, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
    ]))
    story.append(tbl)
    story.append(PageBreak())

    # === GROUPED MODULE SECTIONS ===
    module_counter = 0
    for group_title, group_question, group_mods in grouped_modules:
        # Group header
        story.append(Paragraph(_replace_unicode_symbols(group_title), styles["GroupTitle"]))
        story.append(Paragraph(f"<i>Question: {_replace_unicode_symbols(group_question)}</i>", styles["GroupQuestion"]))

        for m in group_mods:
            module_counter += 1

            # Module header
            story.append(Paragraph(f"{module_counter}. {_replace_unicode_symbols(m.title)}", styles["ModuleTitle"]))
            story.append(Paragraph(
                f"<b>Module ID:</b> <font name='Courier'>{_replace_unicode_symbols(m.module_id)}</font> | "
                f"<b>Output:</b> <font name='Courier'>{_relpath_display(m.out_dir)}</font>",
                styles["BodyCustom"]
            ))
            story.append(Spacer(1, 0.15 * cm))

            # PyR@TE model info for RG modules
            if m.module_id in RG_MODULES:
                pyrate_info = m.get_pyrate_info()
                if pyrate_info:
                    info_lines = ["<b>PyR@TE Model Configuration:</b>"]
                    if pyrate_info.get("model_name"):
                        info_lines.append(f"Model: <font name='Courier'>{pyrate_info['model_name']}</font>")
                    if pyrate_info.get("yaml_source"):
                        info_lines.append(f"YAML: <font name='Courier'>{pyrate_info['yaml_source']}</font>")
                    if pyrate_info.get("yaml_sha256"):
                        info_lines.append(f"SHA256: <font name='Courier'>{pyrate_info['yaml_sha256']}</font>")
                    if pyrate_info.get("pythonoutput_file"):
                        info_lines.append(f"PythonOutput: <font name='Courier'>{pyrate_info['pythonoutput_file']}</font>")
                    elif pyrate_info.get("pythonoutput_dir"):
                        info_lines.append(f"PythonOutput: <font name='Courier'>{pyrate_info['pythonoutput_dir']}</font>")
                    story.append(Paragraph("<br/>".join(info_lines), styles["PyRateInfo"]))
                    story.append(Spacer(1, 0.15 * cm))

            spec = m.spec

            # Try to get docs for unconventional modules
            doc_path = _workspace_root() / "tfpt-suite" / "unconventional" / "docs" / m.module_id / "README.md"
            doc_text = _read_text_if_exists(doc_path) if doc_path.is_file() else None

            def add_text_block(*, title: str, text: str) -> None:
                from xml.sax.saxutils import escape
                # Replace Unicode symbols and format text - do this BEFORE any other processing
                text = _replace_unicode_symbols(str(text))
                # Apply again after sanitization to catch any issues
                text = _sanitize_paths(text)
                text = _replace_unicode_symbols(text)
                safe = escape(_zwsp_break(text)).replace("\n", "<br/>")
                story.append(Paragraph(_replace_unicode_symbols(title), styles["SectionHead"]))
                story.append(Paragraph(safe, styles["BodyCustom"]))
                story.append(Spacer(1, 0.1 * cm))

            def spec_list(key: str) -> list[str]:
                v = spec.get(key)
                if isinstance(v, list):
                    return [str(x) for x in v if str(x).strip()]
                return []

            def spec_str(key: str) -> str | None:
                v = spec.get(key)
                if v is None:
                    return None
                s = str(v).strip()
                return s if s else None

            # Render goal/method from docs or spec
            if doc_text:
                goals_sections = _extract_doc_sections(doc_text, wanted_prefixes=[
                    "purpose", "why this matters", "what it is", "what this suite is",
                ])
                method_sections = _extract_doc_sections(doc_text, wanted_prefixes=[
                    "what the module computes", "what it computes", "what it searches",
                    "holdout split", "toy model used", "scale policy", "inputs", "outputs",
                ])

                if goals_sections:
                    goal_text = "\n\n".join([f"{_replace_unicode_symbols(t)}: {_replace_unicode_symbols(body)}" for t, body in goals_sections if body.strip()])
                    if goal_text:
                        add_text_block(title="Objective", text=goal_text)

                if method_sections:
                    method_text = "\n\n".join([f"{_replace_unicode_symbols(t)}: {_replace_unicode_symbols(body)}" for t, body in method_sections if body.strip()])
                    if method_text:
                        add_text_block(title="Methodology", text=method_text)

            elif spec:
                # Fallback to spec
                goals: list[str] = []
                question = spec_str("question")
                if question:
                    goals.append(_replace_unicode_symbols(question))
                objective = spec_list("objective")
                if objective:
                    goals.extend([_replace_unicode_symbols(str(x)) for x in objective])
                if goals:
                    add_text_block(title="Objective", text="; ".join(goals))

                method_lines: list[str] = []
                if spec.get("inputs"):
                    inputs_text = ", ".join([_replace_unicode_symbols(str(x)) for x in spec.get("inputs", [])])
                    method_lines.append("Inputs: " + inputs_text)
                if spec.get("formulas"):
                    formulas_text = "; ".join([_replace_unicode_symbols(str(x)) for x in spec.get("formulas", [])])
                    method_lines.append("Formulas: " + formulas_text)
                assumptions = spec_list("assumptions")
                if assumptions:
                    assumptions_text = "; ".join([_replace_unicode_symbols(str(x)) for x in assumptions])
                    method_lines.append("Assumptions: " + assumptions_text)
                if method_lines:
                    add_text_block(title="Methodology", text="\n".join(method_lines))

            # Checks table
            story.append(Paragraph("Checks", styles["SectionHead"]))

            check_rows: list[list[Any]] = [[
                P("Check", style=cell_bold),
                P("Severity", style=cell_bold),
                P("Detail", style=cell_bold),
            ]]

            for c in m.checks:
                sev = str(c.get("severity") or "").strip().upper()
                if not sev:
                    sev = "PASS" if _truthy(c.get("passed")) else "FAIL"
                detail = str(c.get("detail", ""))
                # Format detail: replace Unicode, format numbers, truncate if needed
                detail = _replace_unicode_symbols(detail)
                detail = _format_numerical_data(detail)
                if len(detail) > 200:
                    detail = detail[:197] + "..."
                check_id = _replace_unicode_symbols(str(c.get("check_id", "")))
                check_rows.append([
                    P(check_id),
                    P(sev),
                    P(_zwsp_break(detail), style=cell_code),
                ])

            chk_fracs = [0.28, 0.12, 0.60]
            chk_widths = [available_w * f for f in chk_fracs]
            check_tbl = Table(check_rows, colWidths=chk_widths, repeatRows=1)

            # Color rows by severity
            table_style = [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#4a5568")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 7),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]

            # Color code severity
            for row_idx, c in enumerate(m.checks, start=1):
                sev = str(c.get("severity") or "").strip().upper()
                if not sev:
                    sev = "PASS" if _truthy(c.get("passed")) else "FAIL"
                if sev == "FAIL":
                    table_style.append(("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor("#fed7d7")))
                elif sev == "WARN":
                    table_style.append(("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor("#fefcbf")))
                elif sev == "PASS":
                    table_style.append(("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor("#c6f6d5")))

            check_tbl.setStyle(TableStyle(table_style))
            story.append(check_tbl)
            story.append(Spacer(1, 0.15 * cm))

            # Warnings
            if m.warnings:
                warnings_text = "\n".join([f"* {w}" for w in m.warnings])
                warnings_text = _replace_unicode_symbols(warnings_text)
                story.append(Paragraph("Warnings", styles["SectionHead"]))
                story.append(Preformatted(_wrap_lines(_sanitize_paths(warnings_text), 120), styles["CodeSmall"]))
                story.append(Spacer(1, 0.1 * cm))

            # Report.txt (collapsed for brevity)
            if m.report_text.strip():
                report_lines = m.report_text.strip().split("\n")
                # Show first 30 lines, then "..."
                if len(report_lines) > 40:
                    display_text = "\n".join(report_lines[:30]) + "\n\n... (truncated, see report.txt for full output) ..."
                else:
                    display_text = m.report_text
                
                # Format the report text for better readability
                # Apply Unicode replacement FIRST, then format
                display_text = _replace_unicode_symbols(display_text)
                display_text = _format_report_text(display_text)
                # Apply one more time after all formatting
                display_text = _replace_unicode_symbols(display_text)
                
                story.append(Paragraph("Results (report.txt)", styles["SectionHead"]))
                story.append(Preformatted(_wrap_lines(_sanitize_paths(display_text), 130), styles["CodeSmall"]))
                story.append(Spacer(1, 0.1 * cm))

            # Plots
            if m.plots:
                story.append(Paragraph("Plots", styles["SectionHead"]))
                for p in m.plots:
                    try:
                        # Ensure path is absolute and file exists
                        if p.is_absolute():
                            plot_path = p
                        else:
                            # Try relative to module output directory
                            plot_path = m.out_dir / p.name
                        
                        # Resolve to absolute path and check existence
                        plot_path = plot_path.resolve()
                        if not plot_path.exists():
                            story.append(Paragraph(f"(Plot file not found: {p.name} at {plot_path})", styles["BodyCustom"]))
                            continue
                        
                        if not plot_path.is_file():
                            story.append(Paragraph(f"(Plot path is not a file: {p.name})", styles["BodyCustom"]))
                            continue
                        
                        img = Image(str(plot_path))
                        img._restrictSize(16.0 * cm, 20.0 * cm)
                        story.append(Paragraph(f"<font name='Courier' size='7'>{plot_path.name}</font>", styles["BodyCustom"]))
                        story.append(img)
                        story.append(Spacer(1, 0.15 * cm))
                    except Exception as e:
                        # More detailed error message for debugging
                        import traceback
                        error_msg = f"(Plot could not be embedded: {p.name}"
                        if hasattr(e, '__class__'):
                            error_msg += f", error: {e.__class__.__name__}: {str(e)}"
                        error_msg += ")"
                        story.append(Paragraph(error_msg, styles["BodyCustom"]))

            # Add separator between modules
            story.append(Spacer(1, 0.3 * cm))

        # Page break after each group
        story.append(PageBreak())

    # === APPENDIX: JSON DATA ===
    story.append(Paragraph("Appendix: JSON Data", styles["GroupTitle"]))
    story.append(Paragraph(
        "Complete results.json per module for machine processing.",
        styles["GroupQuestion"]
    ))

    for m in mods:
        story.append(Paragraph(f"<b>{m.module_id}</b>", styles["SectionHead"]))
        json_text = json.dumps(m.results_payload, indent=2, sort_keys=True, ensure_ascii=True)  # Use ensure_ascii=True to avoid Unicode issues
        # Format JSON text: replace Unicode and format numbers (though ensure_ascii=True should handle most)
        json_text = _replace_unicode_symbols(json_text)
        json_text = _format_numerical_data(json_text)
        json_text = _replace_unicode_symbols(json_text)  # Apply again after formatting
        story.append(Preformatted(_wrap_lines(_sanitize_paths(json_text), 130), styles["CodeSmall"]))
        story.append(Spacer(1, 0.3 * cm))

    def _add_page_numbers(canvas, _doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#4a5568"))
        page_num = canvas.getPageNumber()
        canvas.drawRightString(A4[0] - 1.8 * cm, 1.2 * cm, f"Page {page_num}")
        canvas.drawString(1.8 * cm, 1.2 * cm, "TFPT Suite Report")
        canvas.restoreState()

    doc.build(story, onFirstPage=_add_page_numbers, onLaterPages=_add_page_numbers)
    return pdf_path
