from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, Paragraph, Preformatted, SimpleDocTemplate, Spacer, Table, TableStyle

from tfpt_suite.report_builder import discover_module_artifacts
from theoryv3_suite.utils import ensure_ascii


MODULE_TEST_MAP: dict[str, dict[str, str]] = {
    "seed_invariants_audit": {
        "idea": "pi seeds the invariants",
        "test": "c3, varphi0_tree, delta_top, varphi0, beta_rad identities (and optional match to core_invariants)",
        "signal": "identities hold and match within numerical tolerance",
    },
    "defect_partition_g5_audit": {
        "idea": "discrete defect partition (g=5) fixes delta2 and alpha_inv(0)",
        "test": "g=5 multiplicity, alpha_inv(0) within 2 sigma, and negative control vs g=4,6",
        "signal": "g equals 5; z <= 2; g=5 minimizes |z|",
    },
    "alpha_backreaction_sensitivity_audit": {
        "idea": "k=2 backreaction exponent is a structural choice, not a fit",
        "test": "sweep k and compute ppm vs CODATA",
        "signal": "k=2 near minimum |ppm|",
    },
    "g5_origin_audit": {
        "idea": "single origin for g from SU(5) holonomy degeneracy",
        "test": "count eigenvalue degeneracies in SU(5) hypercharge spectrum",
        "signal": "degeneracies {3,2} and g=5",
    },
    "dark_energy_exponential_audit": {
        "idea": "exp(-alpha_inv/2) suppression sets dark energy scale",
        "test": "phi_star_base and discrete normalization candidates vs rho_L target",
        "signal": "best candidate within 0.5 dex, n=1/2 preferred",
    },
    "dark_energy_norm_half_origin_audit": {
        "idea": "n=1/2 fixed by double cover (k=2)",
        "test": "n_from_cover=1/k and best candidate label",
        "signal": "n_from_cover=1/2 and best label n=1/2",
    },
    "flavor_pattern_audit": {
        "idea": "Mobius/Z3 flavor anchors from varphi0 and delta_star",
        "test": "lambda (Cabibbo), sin2(theta13), cusp set {1,1/3,2/3}",
        "signal": "lambda and sin2(theta13) within 2 sigma; cusp set matches",
    },
    "pmns_tm1_audit": {
        "idea": "TM1 sum rule for theta12",
        "test": "sin2(theta12) from sin2(theta13)",
        "signal": "sin2(theta12) within conservative band",
    },
    "constant_factory_audit": {
        "idea": "constant factory: hierarchical TFPT constants from simple rules",
        "test": "compute constants, sensitivities, and compare to reference ledger",
        "signal": "derivations + sensitivities documented; references and crosslinks visible",
    },
    "yukawa_exponent_index_audit": {
        "idea": "mass ratios follow rational exponent indices q_ij",
        "test": "q_ij rationalization and ratio reconstruction error",
        "signal": "max relative error <= 2%",
    },
    "yukawa_index_mapping_audit": {
        "idea": "q_ij map to charge-squared index sums",
        "test": "bounded integer sum of charge-squared indices",
        "signal": "mapping errors <= 2%",
    },
    "baryon_consistency_audit": {
        "idea": "beta_rad anchors Omega_b and eta_b",
        "test": "Omega_b identity, eta_b proxy, derived H0",
        "signal": "Omega_b within 2 sigma; eta_b within 0.5 dex; H0 within 2 sigma",
    },
    "axion_dm_audit": {
        "idea": "axion DM target frequency and relic fraction",
        "test": "nu from m_a conversion, Omega_a h^2 vs ref",
        "signal": "frequency matches conversion; relic fraction near ref",
    },
    "g5_crosslink_audit": {
        "idea": "g=5 crosslinks across sectors",
        "test": "delta2 factor g/4, unification patch g/2, gamma0=g/(g+1)",
        "signal": "g=5; delta_b3 includes g/2; gamma0 matches g/(g+1)",
    },
}

# Table layout constants
MAP_COL_FRACS = [0.18, 0.26, 0.36, 0.20]
OVERVIEW_COL_FRACS = [0.42, 0.25, 0.33]
CONST_COL_FRACS = [0.18, 0.36, 0.12, 0.18, 0.08, 0.08]

DEV_Z_THRESHOLDS = (1.0, 2.0, 3.0)
DEV_REL_THRESHOLDS = (0.01, 0.05, 0.1)
DEV_COLOR_PASS = colors.HexColor("#c6f6d5")
DEV_COLOR_WARN = colors.HexColor("#fefcbf")
DEV_COLOR_ALERT = colors.HexColor("#fed7aa")
DEV_COLOR_FAIL = colors.HexColor("#fed7d7")
DEV_COLOR_NA = colors.HexColor("#edf2f7")


EXPECTED_PLOTS: dict[str, list[str]] = {
    "seed_invariants_audit": ["seed_invariants.png"],
    "defect_partition_g5_audit": ["alpha_defect_series.png", "g_negative_control.png"],
    "alpha_backreaction_sensitivity_audit": ["alpha_backreaction_ppm.png"],
    "g5_origin_audit": ["g5_origin.png"],
    "dark_energy_exponential_audit": ["rho_lambda_candidates.png"],
    "dark_energy_norm_half_origin_audit": ["dark_energy_norm_origin.png"],
    "flavor_pattern_audit": ["flavor_anchors.png", "mobius_ratios.png"],
    "yukawa_exponent_index_audit": ["yukawa_q_errors.png"],
    "yukawa_index_mapping_audit": ["yukawa_index_mapping.png"],
    "baryon_consistency_audit": ["baryon_consistency.png"],
    "axion_dm_audit": ["axion_summary.png"],
    "g5_crosslink_audit": ["g5_links.png"],
    "constant_factory_audit": ["constant_factory_summary.png", "constant_factory_sensitivity.png"],
}


def _check_counts(checks: list[dict[str, Any]]) -> tuple[int, int, int]:
    passed = 0
    warned = 0
    failed = 0
    for c in checks:
        sev = str(c.get("severity") or "").strip().upper()
        if not sev:
            sev = "PASS" if bool(c.get("passed")) else "FAIL"
        if sev == "PASS":
            passed += 1
        elif sev == "WARN":
            warned += 1
        elif sev == "FAIL":
            failed += 1
        else:
            passed += 1
    return passed, warned, failed


def _p(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(escape(ensure_ascii(str(text))), style)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_repo_root()))
    except Exception:
        return str(path)


def _sanitize_text(text: str) -> str:
    root = str(_repo_root())
    return ensure_ascii(str(text)).replace(root, ".")


def _fmt_num(value: Any, digits: int = 6) -> str:
    if value is None:
        return "n/a"
    try:
        v = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(v):
        return "n/a"
    return f"{v:.{digits}g}"


def _format_reference(ref: dict[str, Any] | None) -> str:
    if not isinstance(ref, dict):
        return "n/a"
    if ref.get("type") and ref.get("origin"):
        return f"{ref.get('type')}: {ref.get('origin')}"
    value = ref.get("value")
    sigma = ref.get("sigma")
    version = ref.get("version", "n/a")
    units = ref.get("units")
    value_text = _fmt_num(value, digits=6)
    sigma_text = _fmt_num(sigma, digits=4) if sigma not in (None, 0) else None
    parts = [value_text]
    if sigma_text:
        parts.append(f"± {sigma_text}")
    if units:
        parts.append(str(units))
    parts.append(f"({version})")
    return " ".join(parts)


def _deviation_label_color(entry: dict[str, Any]) -> tuple[str, colors.Color]:
    status = str(entry.get("status", "")).lower()
    if status in {"pending", "placeholder", "input"}:
        return "n/a", DEV_COLOR_NA
    comparison = entry.get("comparison") or {}
    if "z" in comparison:
        try:
            z_val = float(comparison["z"])
        except Exception:
            return "n/a", DEV_COLOR_NA
        z_abs = abs(z_val)
        if z_abs <= DEV_Z_THRESHOLDS[0]:
            color = DEV_COLOR_PASS
        elif z_abs <= DEV_Z_THRESHOLDS[1]:
            color = DEV_COLOR_WARN
        elif z_abs <= DEV_Z_THRESHOLDS[2]:
            color = DEV_COLOR_ALERT
        else:
            color = DEV_COLOR_FAIL
        return f"z={_fmt_num(z_val, digits=3)}", color
    if "diff" in comparison and comparison.get("ref_value") not in (None, 0):
        try:
            diff_val = float(comparison["diff"])
            ref_val = float(comparison["ref_value"])
        except Exception:
            return "n/a", DEV_COLOR_NA
        rel = diff_val / ref_val if ref_val else None
        if rel is None or not math.isfinite(rel):
            return "n/a", DEV_COLOR_NA
        rel_abs = abs(rel)
        if rel_abs <= DEV_REL_THRESHOLDS[0]:
            color = DEV_COLOR_PASS
        elif rel_abs <= DEV_REL_THRESHOLDS[1]:
            color = DEV_COLOR_WARN
        elif rel_abs <= DEV_REL_THRESHOLDS[2]:
            color = DEV_COLOR_ALERT
        else:
            color = DEV_COLOR_FAIL
        return f"rel={_fmt_num(rel * 100.0, digits=3)}%", color
    return "n/a", DEV_COLOR_NA


def _render_constant_factory_tables(
    *, story: list[Any], payload: dict[str, Any], styles: dict[str, ParagraphStyle], available_w: float
) -> None:
    results = payload.get("results", {}) if isinstance(payload, dict) else {}
    groups = results.get("groups", [])
    if not isinstance(groups, list) or not groups:
        return
    story.append(Paragraph("Constant tables (grouped)", styles["SubHead"]))
    col_widths = [available_w * f for f in CONST_COL_FRACS]
    for group in groups:
        group_name = str(group.get("group", "n/a"))
        items = group.get("items", [])
        story.append(Paragraph(ensure_ascii(group_name), styles["SubHead"]))
        rows: list[list[Any]] = [[
            _p("Constant", styles["ConstCell"]),
            _p("Formula / Calculation", styles["ConstCell"]),
            _p("Value", styles["ConstCell"]),
            _p("Reference", styles["ConstCell"]),
            _p("Deviation", styles["ConstCell"]),
            _p("Status", styles["ConstCell"]),
        ]]
        table_style_cmds = [
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 6.5),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]

        if isinstance(items, list):
            for idx, item in enumerate(items, start=1):
                label = f"{item.get('label', 'n/a')} ({item.get('key', 'n/a')})"
                formula_math = item.get("formula_math") or item.get("formula") or ""
                formula_text = item.get("formula_text")
                calc_lines = []
                if formula_math:
                    calc_lines.append(formula_math)
                if formula_text:
                    calc_lines.append(formula_text)
                source_module = item.get("source_module_id")
                source_pointer = item.get("source_json_pointer")
                if source_module:
                    detail = f"{source_module}"
                    if source_pointer:
                        detail = f"{detail} {source_pointer}"
                    calc_lines.append(f"source: {detail}")
                calc = "<br/>".join([ensure_ascii(str(x)) for x in calc_lines if str(x).strip()])
                value_text = _fmt_num(item.get("value"), digits=8)
                units = item.get("units")
                if units:
                    value_text = f"{value_text} {units}"
                ref_text = _format_reference(item.get("reference") or item.get("reference_info"))
                dev_label, dev_color = _deviation_label_color(item)
                status = ensure_ascii(str(item.get("status_label") or item.get("status") or "n/a"))

                rows.append([
                    _p(label, styles["ConstCell"]),
                    _p(calc or "n/a", styles["ConstCell"]),
                    _p(value_text, styles["ConstCell"]),
                    _p(ref_text, styles["ConstCell"]),
                    _p(dev_label, styles["ConstCell"]),
                    _p(status, styles["ConstCell"]),
                ])
                table_style_cmds.append(("BACKGROUND", (4, idx), (4, idx), dev_color))

        table = Table(rows, repeatRows=1, colWidths=col_widths)
        table.setStyle(TableStyle(table_style_cmds))
        story.append(table)
        story.append(Spacer(1, 0.2 * cm))

        grammar_flags = []
        for item in items:
            grammar = item.get("grammar") if isinstance(item, dict) else None
            if isinstance(grammar, dict) and not grammar.get("allowed", True):
                issues = ", ".join([str(x) for x in grammar.get("issues", [])])
                grammar_flags.append(f"{item.get('key', 'n/a')}: {issues}")
        if grammar_flags:
            story.append(Paragraph("Grammar warnings: " + "; ".join(grammar_flags), styles["BodySmall"]))
            story.append(Spacer(1, 0.1 * cm))

    crosslinks = results.get("crosslinks", [])
    if isinstance(crosslinks, list) and crosslinks:
        story.append(Paragraph("Crosslink signatures", styles["SubHead"]))
        rows = [[_p("Signature", styles["ConstCell"]), _p("Details", styles["ConstCell"])]]
        for link in crosslinks:
            if not isinstance(link, dict):
                continue
            kind = str(link.get("kind", "signature"))
            detail_parts = []
            for key, value in link.items():
                if key == "kind":
                    continue
                detail_parts.append(f"{key}={value}")
            rows.append([_p(kind, styles["ConstCell"]), _p(", ".join(detail_parts), styles["ConstCell"])])
        table = Table(rows, repeatRows=1, colWidths=[available_w * 0.25, available_w * 0.75])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 6.5),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(table)
        story.append(Spacer(1, 0.2 * cm))

    views = results.get("views", {})
    if isinstance(views, dict) and views.get("anchor") and views.get("tfpt"):
        story.append(Paragraph("Anchor vs TFPT view", styles["SubHead"]))
        anchor = views.get("anchor", {})
        tfpt = views.get("tfpt", {})
        rows = [[_p("Quantity", styles["ConstCell"]), _p("Anchor", styles["ConstCell"]), _p("TFPT", styles["ConstCell"])]]
        for key in ["H0", "Omega_b", "eta_b", "Omega_dm"]:
            rows.append([
                _p(key, styles["ConstCell"]),
                _p(_fmt_num(anchor.get(key), digits=6), styles["ConstCell"]),
                _p(_fmt_num(tfpt.get(key), digits=6), styles["ConstCell"]),
            ])
        note = views.get("note")
        table = Table(rows, repeatRows=1, colWidths=[available_w * 0.2, available_w * 0.4, available_w * 0.4])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 6.5),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(table)
        if note:
            story.append(Paragraph(str(note), styles["BodySmall"]))
        story.append(Spacer(1, 0.2 * cm))

    gaps = []
    for group in groups:
        items = group.get("items", []) if isinstance(group, dict) else []
        for item in items:
            status = str(item.get("status", "")).lower()
            if status in {"pending", "placeholder"}:
                gaps.append(item)
    if gaps:
        story.append(Paragraph("Gap list (pending / placeholder)", styles["SubHead"]))
        rows = [[_p("Constant", styles["ConstCell"]), _p("Reason", styles["ConstCell"])]]
        for item in gaps:
            label = f"{item.get('label', 'n/a')} ({item.get('key', 'n/a')})"
            reason = item.get("formula_text") or item.get("note") or item.get("status") or "pending"
            rows.append([_p(label, styles["ConstCell"]), _p(str(reason), styles["ConstCell"])])
        table = Table(rows, repeatRows=1, colWidths=[available_w * 0.35, available_w * 0.65])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 6.5),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(table)
        story.append(Spacer(1, 0.2 * cm))

def _validate_plots(mod, expected: list[str]) -> dict[str, Any]:
    present: list[str] = []
    missing: list[str] = []
    invalid: list[str] = []
    extra: list[str] = []
    load_check_enabled = True

    by_name = {p.name: p for p in mod.plots}
    for name in expected:
        path = by_name.get(name)
        if not path or not path.exists():
            missing.append(name)
            continue
        try:
            from PIL import Image as PILImage  # type: ignore

            with PILImage.open(path) as img:
                img.verify()
        except Exception:
            invalid.append(name)
        else:
            present.append(name)

    if expected:
        extra = sorted([name for name in by_name.keys() if name not in expected])

    # If Pillow is not available, do not mark plots invalid; treat as present if they exist.
    try:
        import PIL  # type: ignore  # noqa: F401
    except Exception:
        load_check_enabled = False
        invalid = []
        present = [name for name in expected if name not in missing]

    return {
        "expected": expected,
        "present": present,
        "missing": missing,
        "invalid": invalid,
        "extra": extra,
        "load_check_enabled": load_check_enabled,
    }


def build_theoryv3_report(*, out_dir: Path, pdf_path: Path) -> Path:
    mods = discover_module_artifacts(out_dir)
    if not mods:
        raise RuntimeError(f"No module outputs found under: {out_dir}")

    pdf_path = pdf_path.resolve()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.8 * cm,
        title="TFPT theoryv3 analysis",
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="TitleCustom",
        parent=styles["Title"],
        fontSize=20,
        spaceAfter=10,
        textColor=colors.HexColor("#1a365d"),
    ))
    styles.add(ParagraphStyle(
        name="ModuleTitle",
        parent=styles["Heading2"],
        fontSize=12,
        spaceBefore=12,
        spaceAfter=6,
        textColor=colors.HexColor("#2d3748"),
    ))
    styles.add(ParagraphStyle(
        name="SectionTitle",
        parent=styles["Heading2"],
        fontSize=12,
        spaceBefore=10,
        spaceAfter=6,
        textColor=colors.HexColor("#2c5282"),
    ))
    styles.add(ParagraphStyle(
        name="SubHead",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
        textColor=colors.HexColor("#2d3748"),
        spaceBefore=4,
        spaceAfter=2,
    ))
    styles.add(ParagraphStyle(
        name="BodySmall",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
    ))
    styles.add(ParagraphStyle(
        name="CodeSmall",
        parent=styles["BodyText"],
        fontName="Courier",
        fontSize=7.5,
        leading=9,
    ))
    styles.add(ParagraphStyle(
        name="MapCell",
        parent=styles["BodyText"],
        fontSize=7,
        leading=8.5,
        wordWrap="CJK",
    ))
    styles.add(ParagraphStyle(
        name="ConstCell",
        parent=styles["BodyText"],
        fontSize=6.5,
        leading=8,
        wordWrap="CJK",
    ))

    plot_checks: dict[str, dict[str, Any]] = {}
    plot_expected_total = 0
    plot_present_total = 0
    plot_missing_total = 0
    plot_invalid_total = 0
    plot_extra_total = 0
    plot_load_check_enabled = True

    pass_total = 0
    warn_total = 0
    fail_total = 0
    for m in mods:
        passed, warned, failed = _check_counts(m.checks)
        pass_total += passed
        warn_total += warned
        fail_total += failed
        expected = EXPECTED_PLOTS.get(m.module_id, [])
        status = _validate_plots(m, expected)
        plot_checks[m.module_id] = status
        plot_expected_total += len(status["expected"])
        plot_present_total += len(status["present"])
        plot_missing_total += len(status["missing"])
        plot_invalid_total += len(status["invalid"])
        plot_extra_total += len(status["extra"])
        plot_load_check_enabled = plot_load_check_enabled and bool(status["load_check_enabled"])

    story: list[Any] = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    story.append(Paragraph("TFPT theoryv3 analysis report", styles["TitleCustom"]))
    story.append(Paragraph(f"Generated: {now}", styles["BodySmall"]))
    story.append(Paragraph(f"Output directory: {_rel_path(out_dir)}", styles["BodySmall"]))
    story.append(Spacer(1, 0.4 * cm))

    story.append(Paragraph("Executive Summary", styles["SectionTitle"]))
    intro_text = (
        "<b>Purpose</b>: This report tests the TFPT theoryv3 hypothesis that complexity collapses into "
        "discrete building blocks: pi seeds the invariants, a discrete defect partition (g=5) fixes delta2 "
        "and alpha_inv(0), exponential suppression exp(-alpha_inv/2) sets the dark energy scale, and Mobius/Z3 "
        "structures fix flavor and mass hierarchies.<br/><br/>"
        "<b>What we test</b>: Each module targets one piece of that thesis and turns it into an explicit, "
        "deterministic check with declared tolerances (no continuous fitting).<br/><br/>"
        "<b>What we are searching for</b>: A small set of discrete identities that reproduces multiple observables "
        "across sectors and cross-links them consistently (especially the g=5 signature).<br/><br/>"
        "<b>How we proceed</b>: Modules consume existing out_physics results when available; otherwise they compute "
        "from the TFPT invariants. Each test is a finite enumeration or closed-form calculation. Plot generation is "
        "validated by file existence and image load checks (when available).<br/><br/>"
        "<b>Result</b>: {modules} modules, {checks} total checks, "
        "PASS={passc}, WARN={warnc}, FAIL={failc}. Plot verification: expected={pexp}, present={ppresent}, "
        "missing={pmiss}, invalid={pinvalid}, extra={pextra}, load_check={pload}."
    ).format(
        modules=len(mods),
        checks=pass_total + warn_total + fail_total,
        passc=pass_total,
        warnc=warn_total,
        failc=fail_total,
        pexp=plot_expected_total,
        ppresent=plot_present_total,
        pmiss=plot_missing_total,
        pinvalid=plot_invalid_total,
        pextra=plot_extra_total,
        pload="on" if plot_load_check_enabled else "off",
    )
    story.append(Paragraph(ensure_ascii(intro_text), styles["BodySmall"]))
    story.append(Spacer(1, 0.3 * cm))

    # Test map table
    story.append(Paragraph("Test Map (idea -> test -> signal)", styles["SectionTitle"]))
    map_rows = [[
        _p("Module", styles["MapCell"]),
        _p("Idea segment", styles["MapCell"]),
        _p("What is checked", styles["MapCell"]),
        _p("Signal we seek", styles["MapCell"]),
    ]]
    for m in mods:
        meta = MODULE_TEST_MAP.get(m.module_id, {})
        map_rows.append([
            _p(m.module_id, styles["MapCell"]),
            _p(meta.get("idea", "n/a"), styles["MapCell"]),
            _p(meta.get("test", "n/a"), styles["MapCell"]),
            _p(meta.get("signal", "n/a"), styles["MapCell"]),
        ])
    available_w = float(A4[0] - doc.leftMargin - doc.rightMargin)
    map_col_widths = [available_w * f for f in MAP_COL_FRACS]
    map_table = Table(map_rows, repeatRows=1, colWidths=map_col_widths)
    map_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
    ]))
    story.append(map_table)
    story.append(Spacer(1, 0.4 * cm))

    # Overview table
    rows = [["Module", "Checks", "Status"]]
    for m in mods:
        passed, warned, failed = _check_counts(m.checks)
        status = "OK"
        if failed:
            status = f"FAIL:{failed}"
        elif warned:
            status = f"WARN:{warned}"
        rows.append([ensure_ascii(m.module_id), f"{passed}/{len(m.checks)}", status])

    overview_col_widths = [available_w * f for f in OVERVIEW_COL_FRACS]
    table = Table(rows, repeatRows=1, colWidths=overview_col_widths)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
    ]))
    story.append(table)
    story.append(Spacer(1, 0.4 * cm))

    story.append(Paragraph("Validated Results Tables", styles["SectionTitle"]))
    story.append(Spacer(1, 0.2 * cm))

    appendix_entries: list[dict[str, str]] = []

    # Module sections
    for m in mods:
        story.append(Paragraph(f"{ensure_ascii(m.title)}", styles["ModuleTitle"]))
        story.append(Paragraph(f"Module ID: {ensure_ascii(m.module_id)}", styles["BodySmall"]))
        meta = MODULE_TEST_MAP.get(m.module_id, {})
        if meta:
            story.append(Paragraph(f"Idea segment: {ensure_ascii(meta.get('idea', 'n/a'))}", styles["BodySmall"]))
            story.append(Paragraph(f"Test focus: {ensure_ascii(meta.get('test', 'n/a'))}", styles["BodySmall"]))
            story.append(Paragraph(f"Signal: {ensure_ascii(meta.get('signal', 'n/a'))}", styles["BodySmall"]))
        plot_status = plot_checks.get(m.module_id, {})
        if plot_status.get("expected"):
            story.append(Paragraph(
                "Plot check: expected={exp}, present={present}, missing={missing}, invalid={invalid}, extra={extra}".format(
                    exp=len(plot_status.get("expected", [])),
                    present=len(plot_status.get("present", [])),
                    missing=",".join(plot_status.get("missing", [])) or "none",
                    invalid=",".join(plot_status.get("invalid", [])) or "none",
                    extra=",".join(plot_status.get("extra", [])) or "none",
                ),
                styles["BodySmall"],
            ))
        story.append(Spacer(1, 0.1 * cm))

        spec = m.spec or {}
        question = spec.get("question") or ""
        objective = spec.get("objective") or []
        inputs = spec.get("inputs") or []
        validation = spec.get("validation") or []
        formulas = spec.get("formulas") or []
        assumptions = spec.get("assumptions") or []
        determinism = spec.get("determinism") or ""

        def _block(title: str, lines: list[str]) -> None:
            story.append(Paragraph(title, styles["SubHead"]))
            if not lines:
                lines = ["n/a"]
            text = "<br/>".join([ensure_ascii(x) for x in lines if str(x).strip()])
            story.append(Paragraph(text, styles["BodySmall"]))

        _block("A) What this test is intended to find", [question] + [f"Objective: {x}" for x in objective])
        declared_inputs = [str(x) for x in inputs]
        results_inputs = []
        results_payload = m.results_payload or {}
        results_inputs_raw = results_payload.get("results", {}).get("inputs", {})
        if isinstance(results_inputs_raw, dict):
            results_inputs = [f"{k}={v}" for k, v in results_inputs_raw.items()]
        refs_block = []
        refs_raw = results_payload.get("results", {}).get("references", {})
        if isinstance(refs_raw, dict):
            for key, ref in refs_raw.items():
                if isinstance(ref, dict):
                    version = ref.get("version", "n/a")
                    value = ref.get("value", "n/a")
                    sigma = ref.get("sigma", "n/a")
                    refs_block.append(f"reference[{key}]: {version} value={value} sigma={sigma}")
        _block("B) Inputs", declared_inputs + results_inputs + refs_block)
        _block("C) Expected result (signal)", [str(x) for x in validation])
        method_lines = [str(x) for x in formulas]
        if assumptions:
            method_lines.append("Assumptions: " + "; ".join([str(x) for x in assumptions]))
        if determinism:
            method_lines.append(f"Determinism: {determinism}")
        _block("D) Method (how the test proceeds)", method_lines)

        # E) Results and JSON outputs
        passed, warned, failed = _check_counts(m.checks)
        check_lines = [f"Summary: PASS={passed}, WARN={warned}, FAIL={failed}"]
        for c in m.checks:
            sev = str(c.get("severity") or "").strip().upper()
            if not sev:
                sev = "PASS" if bool(c.get("passed")) else "FAIL"
            check_lines.append(f"{sev}: {c.get('check_id')} - {c.get('detail')}")
        warnings_text = [f"Warning: {w}" for w in m.warnings] if m.warnings else []
        plot_line = []
        if plot_status.get("expected"):
            plot_line = [f"Plot check: missing={plot_status.get('missing', [])}, invalid={plot_status.get('invalid', [])}, extra={plot_status.get('extra', [])}"]
        _block("E) Results and JSON outputs", check_lines + warnings_text + plot_line)

        if m.module_id == "constant_factory_audit":
            _render_constant_factory_tables(
                story=story,
                payload=m.results_payload or {},
                styles=styles,
                available_w=available_w,
            )

        if m.report_text.strip():
            appendix_entries.append({
                "module_id": ensure_ascii(m.module_id),
                "report_text": _sanitize_text(m.report_text),
            })

        json_text = json.dumps(m.results_payload, indent=2, ensure_ascii=True)
        appendix_entries.append({
            "module_id": ensure_ascii(m.module_id),
            "results_path": _rel_path(m.out_dir / "results.json"),
            "json_text": _sanitize_text(json_text),
        })

        if m.plots:
            story.append(Paragraph("Plots", styles["BodySmall"]))
            for p in m.plots:
                try:
                    img = Image(str(p))
                    img._restrictSize(16.0 * cm, 20.0 * cm)
                    story.append(Paragraph(ensure_ascii(p.name), styles["BodySmall"]))
                    story.append(img)
                    story.append(Spacer(1, 0.1 * cm))
                except Exception:
                    story.append(Paragraph(f"(plot could not be embedded: {ensure_ascii(p.name)})", styles["BodySmall"]))

        story.append(Spacer(1, 0.3 * cm))

    story.append(Paragraph("Appendix: Raw JSON", styles["SectionTitle"]))
    story.append(Spacer(1, 0.2 * cm))
    for entry in appendix_entries:
        module_id = entry.get("module_id", "module")
        if entry.get("report_text"):
            story.append(Paragraph(f"{module_id}: report.txt", styles["SubHead"]))
            story.append(Preformatted(entry["report_text"], styles["CodeSmall"]))
            story.append(Spacer(1, 0.1 * cm))
        if entry.get("json_text"):
            path_text = entry.get("results_path", "results.json")
            story.append(Paragraph(f"{module_id}: results.json ({path_text})", styles["SubHead"]))
            story.append(Preformatted(entry["json_text"], styles["CodeSmall"]))
            story.append(Spacer(1, 0.2 * cm))

    doc.build(story)
    return pdf_path
