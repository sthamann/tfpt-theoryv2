from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import MPL_REDUCED_GEV, V_star_from_As_r
from tfpt_suite.defect_partition import derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import coerce_float, ensure_ascii, extract_nested, load_tfpt_results, read_json_if_exists, safe_log10


@dataclass(frozen=True)
class ConstantItem:
    key: str
    label: str
    formula_math: str
    compute: Callable[[dict[str, float]], float | None]
    formula_text: str | None = None
    reference_id: str | None = None
    reference_info: dict[str, str] | None = None
    note: str | None = None
    inputs: list[str] | None = None
    depends_on: list[str] | None = None
    status: str = "derived"
    units: str | None = None
    source: str | None = None
    source_module_id: str | None = None
    source_artifact: str | None = None
    source_json_pointer: str | None = None
    lineage: str = "axiom"
    fixed_point: bool | None = None
    latex_ref: str | None = None
    primitives: list[str] | None = None
    sensitivity_mode: str | None = None
    sensitivity_override: tuple[float | None, float | None] | None = None
    approximation: str | None = None
    reference_role: str | None = None
    unit_system: str = "natural_units"


def _load_out_result(config, module_id: str) -> dict[str, Any]:
    path = Path(getattr(config, "output_dir", "")) / module_id / "results.json"
    return read_json_if_exists(path) or load_tfpt_results(module_id, prefer_physics=True) or {}


def _plot_constant_hierarchy(
    *, out_dir: Path, entries: list[dict[str, Any]]
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"constant_factory_summary.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels: list[str] = []
        values: list[float] = []
        for row in entries:
            value = row.get("value")
            if value is None:
                continue
            try:
                v = float(value)
            except Exception:
                continue
            if v <= 0:
                continue
            labels.append(row.get("key", ""))
            values.append(safe_log10(v))

        if not labels:
            return plot, warnings

        fig, ax = plt.subplots(figsize=(9.5, 4.0))
        ax.bar(labels, values, color="#2c5282")
        ax.set_ylabel("log10(value)")
        ax.set_title("TFPT constant factory (computed entries)")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.tick_params(axis="x", rotation=20)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "constant_factory_summary.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["constant_factory_summary.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


def _plot_sensitivity_heatmap(
    *, out_dir: Path, entries: list[dict[str, Any]]
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"constant_factory_sensitivity.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        rows: list[tuple[str, float, float]] = []
        for row in entries:
            sens = row.get("sensitivity") or {}
            sc3 = _num(sens.get("Sc3"))
            sphi0 = _num(sens.get("Sphi0"))
            if sc3 is None and sphi0 is None:
                continue
            status = str(row.get("status", ""))
            if status in {"pending", "placeholder", "input"}:
                continue
            rows.append((str(row.get("key", "")), abs(sc3 or 0.0), abs(sphi0 or 0.0)))

        if not rows:
            return plot, warnings

        labels = [r[0] for r in rows]
        data = [[r[1], r[2]] for r in rows]
        fig_h = max(3.0, 0.25 * len(labels))
        fig, ax = plt.subplots(figsize=(6.0, fig_h))
        im = ax.imshow(data, aspect="auto", cmap="viridis")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["|Sc3|", "|Sphi0|"], fontsize=7)
        ax.set_title("Sensitivity heatmap (abs values)", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "constant_factory_sensitivity.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["constant_factory_sensitivity.png"] = str(path)
    except Exception as exc:
        warnings.append(f"sensitivity_plot_failed: {exc}")
    return plot, warnings


def _prediction_value(payload: dict[str, Any], key: str) -> float | None:
    results = payload.get("results", {}) if isinstance(payload, dict) else {}
    predictions = results.get("predictions", [])
    if isinstance(predictions, list):
        for row in predictions:
            if isinstance(row, dict) and row.get("key") == key:
                try:
                    return float(row.get("value"))
                except Exception:
                    return None
    return None


def _const(value: float | None) -> Callable[[dict[str, float]], float | None]:
    return lambda _ctx, value=value: value


def _ctx_value(key: str) -> Callable[[dict[str, float]], float | None]:
    return lambda ctx: ctx.get(key)


def _num(value: Any) -> float | None:
    v = coerce_float(value, default=float("nan"))
    if not math.isfinite(v):
        return None
    return float(v)


def _fmt(value: float | None, *, digits: int = 8) -> str:
    if value is None:
        return "n/a"
    try:
        v = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(v):
        return "n/a"
    return f"{v:.{digits}g}"


def _calc_text(item: ConstantItem, entry: dict[str, Any]) -> str:
    lines: list[str] = [f"formula_math: {item.formula_math}"]
    if item.formula_text:
        lines.append(f"formula_text: {item.formula_text}")
    inputs = entry.get("inputs", {})
    if isinstance(inputs, dict) and inputs:
        parts = [f"{k}={_fmt(_num(v), digits=6)}" for k, v in inputs.items()]
        lines.append("inputs: " + ", ".join(parts))
    depends_on = entry.get("depends_on", [])
    if depends_on:
        lines.append("depends_on: " + ", ".join([str(x) for x in depends_on]))
    if entry.get("value") is not None:
        lines.append(f"value: {_fmt(_num(entry.get('value')), digits=10)}")
    if item.note:
        lines.append(f"note: {item.note}")
    if item.source:
        lines.append(f"source: {item.source}")
    return "\n".join(lines)


def _dataset_value(dataset_id: str) -> float | None:
    try:
        ref = get_dataset(dataset_id)
    except Exception:
        return None
    return _num(ref.get("value"))


def _alpha_inv_from_c3_varphi0(
    *, c3: float, varphi0: float, delta_top: float, delta2: float, mp_dps: int
) -> float | None:
    old = mp.dps
    try:
        mp.dps = int(mp_dps)
        c3_m = mp.mpf(c3)
        varphi0_m = mp.mpf(varphi0)
        delta_top_m = mp.mpf(delta_top)
        delta2_m = mp.mpf(delta2)
        b1 = mp.mpf(41) / 10
        varphi_tree = varphi0_m - delta_top_m

        if not mp.isfinite(varphi_tree) or varphi_tree <= 0:
            return None

        def cfe(alpha: mp.mpf, varphi: mp.mpf) -> mp.mpf:
            return alpha**3 - mp.mpf(2) * (c3_m**3) * alpha**2 - mp.mpf(8) * b1 * (c3_m**6) * mp.log(mp.mpf(1) / varphi)

        def solve_cfe_for(varphi: mp.mpf) -> mp.mpf:
            f = lambda a: cfe(a, varphi)
            return mp.findroot(f, (mp.mpf("0.006"), mp.mpf("0.010")))

        a = solve_cfe_for(varphi0_m)
        for _ in range(80):
            v = varphi_tree + delta_top_m * mp.e ** (-mp.mpf(2) * a) + delta2_m * mp.e ** (-mp.mpf(4) * a)
            nxt = solve_cfe_for(v)
            if abs(nxt - a) < mp.mpf("1e-30"):
                return float(mp.mpf(1) / nxt)
            a = nxt
        return float(mp.mpf(1) / a)
    except Exception:
        return None
    finally:
        mp.dps = old


def _build_context(
    *,
    c3: float,
    varphi0: float,
    pi: float,
    g_value: float,
    k_value: float,
    delta_top: float,
    delta2: float,
    mp_dps: int,
    external: dict[str, float | None],
) -> dict[str, float]:
    varphi0_tree = varphi0 - delta_top
    beta_rad = varphi0 / (4 * pi)
    beta_deg = (180.0 / pi) * beta_rad
    delta_star = (3.0 / 5.0) + varphi0 / 6.0
    gamma0 = g_value / (g_value + 1.0)
    alpha_inv0 = _alpha_inv_from_c3_varphi0(c3=c3, varphi0=varphi0, delta_top=delta_top, delta2=delta2, mp_dps=mp_dps)
    alpha0 = (1.0 / alpha_inv0) if alpha_inv0 and alpha_inv0 > 0 else None
    phi_star_base = math.exp(-alpha_inv0 / 2.0) if alpha_inv0 else None
    v_over_mpl = None
    v_ew = None
    if alpha_inv0 is not None and math.isfinite(alpha_inv0):
        v_over_mpl = (4 * math.pi) * g_value * (c3**2) * math.exp(-(alpha_inv0 / 4.0))
        v_ew = float(v_over_mpl * float(MPL_REDUCED_GEV))
    y_tau = math.pi * (varphi0**2)
    y_mu = math.pi * (varphi0**3)
    y_e = 2 * math.pi * (varphi0**5)
    m_tau = (v_ew * y_tau / math.sqrt(2)) if v_ew is not None else None
    m_mu = (v_ew * y_mu / math.sqrt(2)) if v_ew is not None else None
    m_e = (v_ew * y_e / math.sqrt(2)) if v_ew is not None else None
    sin2_theta13 = varphi0 * math.exp(-(5.0 / 6.0))
    sin2_theta12 = (1.0 / 3.0) * (1.0 - 2.0 * sin2_theta13)
    lambda_pred = math.sqrt(varphi0) * (1.0 - varphi0 / 2.0)
    delta_cp_deg_formula = 180.0 * (1.0 - delta_star)
    Mpl_unreduced = float(math.sqrt(8 * math.pi) * float(MPL_REDUCED_GEV))
    M_over_Mpl = float(math.sqrt(8 * math.pi) * (float(c3) ** 4))

    ctx: dict[str, float] = {
        "pi": float(pi),
        "c3": float(c3),
        "varphi0": float(varphi0),
        "varphi0_tree": float(varphi0_tree),
        "delta_top": float(delta_top),
        "delta2": float(delta2),
        "beta_rad": float(beta_rad),
        "beta_deg": float(beta_deg),
        "delta_star": float(delta_star),
        "gamma0": float(gamma0),
        "g": float(g_value),
        "k": float(k_value),
        "alpha_inv0": float(alpha_inv0) if alpha_inv0 is not None else float("nan"),
        "alpha0": float(alpha0) if alpha0 is not None else float("nan"),
        "phi_star_base": float(phi_star_base) if phi_star_base is not None else float("nan"),
        "v_over_mpl": float(v_over_mpl) if v_over_mpl is not None else float("nan"),
        "v_ew": float(v_ew) if v_ew is not None else float("nan"),
        "y_tau": float(y_tau),
        "y_mu": float(y_mu),
        "y_e": float(y_e),
        "m_tau": float(m_tau) if m_tau is not None else float("nan"),
        "m_mu": float(m_mu) if m_mu is not None else float("nan"),
        "m_e": float(m_e) if m_e is not None else float("nan"),
        "sin2_theta13": float(sin2_theta13),
        "sin2_theta12": float(sin2_theta12),
        "lambda_pred": float(lambda_pred),
        "delta_cp_deg": float(delta_cp_deg_formula),
        "Mpl_reduced": float(MPL_REDUCED_GEV),
        "Mpl_unreduced": Mpl_unreduced,
        "M_over_Mpl": M_over_Mpl,
    }
    for key, value in external.items():
        if value is None:
            continue
        ctx[key] = float(value)
    return ctx


def _log_value(value: float) -> float | None:
    if not math.isfinite(value) or value == 0:
        return None
    return math.log(abs(value))


def _numeric_sensitivity(
    *,
    compute: Callable[[dict[str, float]], float | None],
    ctx_base: dict[str, float],
    ctx_c3_up: dict[str, float],
    ctx_c3_down: dict[str, float],
    ctx_phi_up: dict[str, float],
    ctx_phi_down: dict[str, float],
    eps: float,
) -> tuple[float | None, float | None, bool]:
    base = compute(ctx_base)
    if base is None:
        return None, None, False
    base_ln = _log_value(float(base))
    if base_ln is None:
        return None, None, False
    use_abs = float(base) < 0

    def ln_or_none(ctx: dict[str, float]) -> float | None:
        val = compute(ctx)
        if val is None:
            return None
        return _log_value(float(val))

    ln_c3_up = ln_or_none(ctx_c3_up)
    ln_c3_down = ln_or_none(ctx_c3_down)
    ln_phi_up = ln_or_none(ctx_phi_up)
    ln_phi_down = ln_or_none(ctx_phi_down)

    sc3 = None
    sphi0 = None
    if ln_c3_up is not None and ln_c3_down is not None:
        sc3 = (ln_c3_up - ln_c3_down) / (2.0 * eps)
    if ln_phi_up is not None and ln_phi_down is not None:
        sphi0 = (ln_phi_up - ln_phi_down) / (2.0 * eps)
    return sc3, sphi0, use_abs


def _derive_hierarchy_levels(items: Iterable[dict[str, Any]]) -> dict[str, int]:
    levels: dict[str, int] = {}

    def level_for(key: str, deps: list[str]) -> int:
        if key in levels:
            return levels[key]
        if not deps:
            levels[key] = 0
            return 0
        max_dep = 0
        for dep in deps:
            dep_level = levels.get(dep)
            if dep_level is None:
                levels[dep] = 0
                dep_level = 0
            max_dep = max(max_dep, dep_level)
        levels[key] = max_dep + 1
        return levels[key]

    for item in items:
        key = str(item.get("key", ""))
        deps = item.get("depends_on", []) if isinstance(item.get("depends_on"), list) else []
        level_for(key, deps)
    return levels


def _check_formula_grammar(formula: str, allowed_symbols: set[str]) -> dict[str, Any]:
    issues: list[str] = []
    cleaned = str(formula)
    fraction_matches = re.findall(r"\b\d+\s*/\s*\d+\b", cleaned)
    cleaned = re.sub(r"\b\d+\s*/\s*\d+\b", "", cleaned)
    for token in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", cleaned):
        if token not in allowed_symbols:
            issues.append(f"unknown_symbol:{token}")
    for num in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", cleaned):
        if "." in num or "e" in num.lower():
            issues.append(f"non_rational_literal:{num}")
    return {
        "allowed": not issues,
        "issues": issues,
        "fractions": fraction_matches,
    }


def _load_symbol_registry(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        registry: dict[str, dict[str, str]] = {}
        current = None
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            if not line.startswith(" ") and line.endswith(":"):
                current = line[:-1].strip()
                registry[current] = {}
                continue
            if current and ":" in line:
                key, value = line.split(":", 1)
                registry[current][key.strip()] = value.strip()
        return registry


def _best_k_from_sensitivity(payload: dict[str, Any]) -> float | None:
    abs_ppm = extract_nested(payload, ["results", "abs_ppm"])
    if isinstance(abs_ppm, dict) and abs_ppm:
        best_key = min(abs_ppm, key=lambda k: float(abs_ppm[k]))
        return _num(best_key)
    return None


def _mobius_ratio(y: float, delta: float) -> float:
    return (y + delta) / (y - delta)


def _wolfenstein_from_matrix(matrix: dict[str, Any]) -> dict[str, float | None]:
    if not isinstance(matrix, dict):
        return {"lambda": None, "A": None, "rho": None, "eta": None}
    Vus = _num(matrix.get("Vus"))
    Vcb = _num(matrix.get("Vcb"))
    Vub = _num(matrix.get("Vub"))
    Vtd = _num(matrix.get("Vtd"))
    if Vus is None or Vcb is None:
        return {"lambda": Vus, "A": None, "rho": None, "eta": None}
    if Vus <= 0:
        return {"lambda": Vus, "A": None, "rho": None, "eta": None}
    A = Vcb / (Vus**2) if Vcb is not None else None
    if A is None or Vub is None or Vtd is None:
        return {"lambda": Vus, "A": A, "rho": None, "eta": None}
    denom = A * (Vus**3)
    if denom <= 0:
        return {"lambda": Vus, "A": A, "rho": None, "eta": None}
    r = Vub / denom
    t = Vtd / denom
    rho = (1.0 - (t * t - r * r)) / 2.0
    eta_sq = r * r - rho * rho
    eta = math.sqrt(eta_sq) if eta_sq > 0 else None
    return {"lambda": Vus, "A": A, "rho": rho, "eta": eta}


def _planck_x_max() -> float:
    # Solve d/dx [x^3 / (exp(x)-1)] = 0 (dimensionless)
    old = mp.dps
    try:
        mp.dps = 50
        f = lambda x: (3 * x**2 * (mp.e**x - 1) - x**3 * mp.e**x) / ((mp.e**x - 1) ** 2)
        root = mp.findroot(f, (mp.mpf("2.0"), mp.mpf("3.0")))
        return float(root)
    finally:
        mp.dps = old


def _latex_snippet(entry: dict[str, Any]) -> str:
    key = entry.get("key", "constant")
    label = entry.get("label", "")
    formula = entry.get("formula_math", "") or entry.get("formula", "")
    formula_text = entry.get("formula_text")
    value = entry.get("value", "n/a")
    units = entry.get("units", "")
    sens = entry.get("sensitivity", {}) if isinstance(entry.get("sensitivity"), dict) else {}
    sc3 = sens.get("Sc3", "n/a")
    sphi0 = sens.get("Sphi0", "n/a")
    lines = [
        f"\\paragraph{{{key} ({label})}}",
        f"Definition: ${{{formula}}}$.\\\\",
    ]
    if formula_text:
        lines.append(f"Meaning: {formula_text}.\\\\")
    lines.extend(
        [
            f"Value: ${{{value}}}\\,\\mathrm{{{units}}}$.\\\\",
            f"Sensitivities: $S_{{c_3}}={sc3},\\;S_{{\\varphi_0}}={sphi0}$.",
        ]
    )
    return "\n".join(lines) + "\n"


class ConstantFactoryAuditModule(TfptModule):
    module_id = "constant_factory_audit"
    title = "Constant factory audit (hierarchical TFPT constants)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT constants (pi, c3, varphi0, delta_top, beta_rad)",
                "defect partition delta2 (g=5)",
                "theoryv3 outputs (baryon, axion, dark energy, flavor)",
                "reference ledger (references.json)",
                "source list: theoryv3/constantfactory.md",
            ],
            outputs=[
                "hierarchical list of constants with formulas",
                "grouped table data with calculation traces",
                "computed values and references (z-score where applicable)",
                "summary plot of log10 magnitudes",
                "sensitivity ledger (Sc3, Sphi0) with propagated sigmas",
                "sensitivity heatmap plot",
                "crosslink summary blocks (g, k)",
            ],
            formulas=[
                "c3 = 1/(8*pi)",
                "varphi0 = 1/(6*pi) + 3/(256*pi^4)",
                "beta_rad = varphi0/(4*pi)",
                "delta2 = (g/4) * delta_top^2 (g=5)",
                "alpha_inv0 from CFE + backreaction",
                "alpha_inv0_simple = 4*pi*exp(1/varphi0)",
                "v/Mpl candidate = (4*pi) * g * c3^2 * exp(-alpha_inv/4)",
                "y_tau = pi*varphi0^2, y_mu = pi*varphi0^3, y_e = 2*pi*varphi0^5",
                "m_l = v * y_l / sqrt(2)",
                "sin2_theta13 = varphi0 * exp(-5/6)",
                "sin2_theta12 = (1/3) * (1 - 2*sin2_theta13)",
            ],
            validation=[
                "reports computed entries with reference comparisons where available",
                "flags missing references as INFO (not a failure)",
                "records grammar violations for non-primitive literals",
            ],
            determinism="Deterministic (closed forms + deterministic module outputs).",
            question="Can TFPT generate a hierarchical constant list from a small discrete grammar without fitting?",
            objective=[
                "Compute a structured set of constants using the simplest TFPT rules.",
                "Document derivations and reference comparisons in one report.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        pi_val = float(c.pi)
        c3_base = float(c.c3)
        varphi0_base = float(c.varphi0)
        delta_top = float(c.delta_top)

        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(delta_top))
        g_value = float(d2.justification.get("effective_multiplicity", 5))
        delta2 = float(d2.delta2)

        mp_dps = int(getattr(config, "mp_dps", 80))
        eps = float(getattr(config, "constant_factory_sensitivity_eps", 1e-6))
        sigma_c3_rel = float(getattr(config, "constant_factory_sigma_c3_rel", 0.0))
        sigma_phi0_rel = float(getattr(config, "constant_factory_sigma_phi0_rel", 0.0))

        alpha_sens = _load_out_result(config, "alpha_backreaction_sensitivity_audit")
        k_value = _best_k_from_sensitivity(alpha_sens) or 2.0

        baryon = _load_out_result(config, "baryon_consistency_audit")
        axion = _load_out_result(config, "axion_dm_audit")
        axion_pipeline = _load_out_result(config, "axion_dm_pipeline")
        dark_energy = _load_out_result(config, "dark_energy_exponential_audit")
        dark_energy_norm = _load_out_result(config, "dark_energy_norm_half_origin_audit")
        flavor = _load_out_result(config, "flavor_pattern_audit")
        tm1 = _load_out_result(config, "pmns_tm1_audit")
        alpha_precision = _load_out_result(config, "alpha_precision_audit")
        rg = _load_out_result(config, "two_loop_rg_fingerprints")
        unification = _load_out_result(config, "unification_gate")
        seesaw = _load_out_result(config, "seesaw_block")
        mass_min = _load_out_result(config, "mass_spectrum_minimal")
        torsion = _load_out_result(config, "torsion_bounds_mapping")
        stability = _load_out_result(config, "stability_unitarity_audit")
        bounce = _load_out_result(config, "bounce_perturbations")
        ckm = _load_out_result(config, "ckm_full_pipeline")
        g5_crosslink = _load_out_result(config, "g5_crosslink_audit")
        defect_partition = _load_out_result(config, "defect_partition_g5_audit")

        omega_b_pred = _num(extract_nested(baryon, ["results", "omega_b_pred"]))
        eta_b_pred = _num(extract_nested(baryon, ["results", "eta_b_pred"]))
        H0_pred = _num(extract_nested(baryon, ["results", "H0_pred"]))

        omega_a_h2 = _num(extract_nested(axion, ["results", "relic", "Omega_a_h2"]))
        nu_ghz = _num(extract_nested(axion, ["results", "axion_claim", "nu_GHz"])) or _num(extract_nested(axion_pipeline, ["results", "axion_claim", "nu_GHz"]))
        f_a_GeV = _num(extract_nested(axion, ["results", "axion_claim", "f_a_GeV"])) or _num(extract_nested(axion_pipeline, ["results", "axion_claim", "f_a_GeV"]))
        g_phys = _num(extract_nested(axion_pipeline, ["results", "axion_claim", "g_phys_GeV_inv_candidate"]))

        rho_L_target = _num(extract_nested(dark_energy, ["results", "targets", "rho_L_target"]))

        alpha_bar5_inv_MZ = _num(extract_nested(alpha_precision, ["results", "secondary_alpha_bar5_MZ", "alpha_bar5_inv_MZ_pred"]))
        alpha3_at_mu_star = _num(extract_nested(rg, ["results", "alpha3", "crossing_varphi0", "alpha3"]))
        mu_star_varphi0 = _num(extract_nested(rg, ["results", "alpha3", "crossing_varphi0", "mu_GeV"]))

        mu_gut = _num(extract_nested(unification, ["results", "best_unification_point", "mu_GeV"]))
        mr_GeV = _num(extract_nested(seesaw, ["results", "paper_v1_06_anchor", "MR_GeV"]))
        mnu3_eV = _num(extract_nested(seesaw, ["results", "derived", "mnu3_from_v_sm_eV"]))

        proton_mass = _dataset_value("m_p_pdg") or _num(extract_nested(mass_min, ["results", "ledger", "placeholders", "proton_mass_GeV", "value"]))
        mnu1_placeholder = _num(extract_nested(mass_min, ["results", "ledger", "placeholders", "neutrino_mass_scale_eV", "value"]))

        torsion_max = None
        bounds = extract_nested(torsion, ["results", "inferred_bounds"], default=[])
        if isinstance(bounds, list):
            for row in bounds:
                v = _num(row.get("inferred_abs_max_S_mu_GeV"))
                if v is not None:
                    torsion_max = v if torsion_max is None else max(torsion_max, v)

        lambda_cross = _num(extract_nested(stability, ["results", "lambda_crossing", "mu_cross_GeV"]))
        k_bounce_t = _num(extract_nested(bounce, ["results", "diagnostics", "k_bounce_t_est_raw"]))

        N_star = float(getattr(config, "starobinsky_N", 56.0))
        M_over_Mpl = float(math.sqrt(8 * math.pi) * (c3_base**4))
        A_s_pred = (N_star**2) / (24 * (math.pi**2)) * (M_over_Mpl**2)
        n_s_pred = 1.0 - 2.0 / N_star if N_star > 0 else None
        r_pred = 12.0 / (N_star**2) if N_star > 0 else None
        V_star_GeV4 = None
        if A_s_pred is not None and r_pred is not None:
            V_star_Mpl4 = V_star_from_As_r(A_s_pred, r_pred)
            if V_star_Mpl4 and math.isfinite(V_star_Mpl4) and V_star_Mpl4 > 0:
                V_star_GeV4 = float(V_star_Mpl4) * (float(MPL_REDUCED_GEV) ** 4)

        omega_b_h2_ref = get_dataset("omega_b_h2_planck2018")
        H0_planck_ref = get_dataset("H0_planck2018")
        omega_dm_h2_ref = get_dataset("Omega_dm_h2_planck2018")
        h_ref = float(H0_planck_ref["value"]) / 100.0
        omega_b_ref_value = float(omega_b_h2_ref["value"]) / (h_ref * h_ref)
        eta10_ref = 273.9 * float(omega_b_h2_ref["value"])
        eta_b_ref_value = eta10_ref * 1e-10
        omega_dm_ref_value = float(omega_dm_h2_ref["value"]) / (h_ref * h_ref)
        omega_dm_pred = omega_a_h2 / (h_ref * h_ref) if omega_a_h2 is not None and h_ref > 0 else None
        h_pred = H0_pred / 100.0 if H0_pred is not None and H0_pred > 0 else None
        omega_dm_tfpt = omega_a_h2 / (h_pred * h_pred) if omega_a_h2 is not None and h_pred is not None and h_pred > 0 else None
        omega_b_tfpt = (4 * math.pi - 1) * float(c.beta_rad)
        views = {
            "anchor": {
                "H0": float(H0_planck_ref["value"]),
                "Omega_b": omega_b_ref_value,
                "eta_b": eta_b_ref_value,
                "Omega_dm": omega_dm_ref_value,
            },
            "tfpt": {
                "H0": H0_pred,
                "Omega_b": omega_b_tfpt,
                "eta_b": eta_b_pred,
                "Omega_dm": omega_dm_tfpt,
            },
            "note": "Anchor uses Planck H0/h; TFPT uses predicted H0 where available.",
        }

        m_e_ref = _dataset_value("m_e_pdg")
        m_mu_ref = _dataset_value("m_mu_pdg")
        m_tau_ref = _dataset_value("m_tau_pdg")
        m_t_ref = _dataset_value("m_t_pdg")
        m_b_ref = _dataset_value("m_b_pdg")
        m_c_ref = _dataset_value("m_c_pdg")
        sin2_thetaW_ref = _dataset_value("sin2_thetaW_mz_pdg")
        sin2_theta23_ref = _dataset_value("pmns_sin2_theta23_nufit53_no")

        me_over_mp_ref = (m_e_ref / proton_mass) if m_e_ref is not None and proton_mass is not None and proton_mass > 0 else None
        mu_over_e_ref = _dataset_value("mass_ratio_mu_over_e_pdg")
        tau_over_mu_ref = _dataset_value("mass_ratio_tau_over_mu_pdg")

        ckm_matrix = extract_nested(ckm, ["results", "alpha_s_sensitivity", "matrix_abs_mu_uv", "central"], default={})
        wolf = _wolfenstein_from_matrix(ckm_matrix if isinstance(ckm_matrix, dict) else {})

        Lambda_target = rho_L_target / (float(MPL_REDUCED_GEV) ** 4) if rho_L_target is not None else None

        external = {
            "N_star": N_star,
            "omega_b_pred": omega_b_pred,
            "eta_b_pred": eta_b_pred,
            "H0_pred": H0_pred,
            "omega_a_h2": omega_a_h2,
            "nu_ghz": nu_ghz,
            "f_a_GeV": f_a_GeV,
            "g_phys": g_phys,
            "alpha_bar5_inv_MZ": alpha_bar5_inv_MZ,
            "alpha3_at_mu_star": alpha3_at_mu_star,
            "mu_star_varphi0": mu_star_varphi0,
            "mu_gut": mu_gut,
            "M_R": mr_GeV,
            "mnu3_eV": mnu3_eV,
            "mnu1_placeholder": mnu1_placeholder,
            "m_p_ref": proton_mass,
            "m_e_ref": m_e_ref,
            "m_mu_ref": m_mu_ref,
            "m_tau_ref": m_tau_ref,
            "torsion_max": torsion_max,
            "lambda_cross": lambda_cross,
            "k_bounce_t": k_bounce_t,
            "A_s": A_s_pred,
            "n_s": n_s_pred,
            "r": r_pred,
            "V_star_GeV4": V_star_GeV4,
            "lambda_ext": _num(extract_nested(flavor, ["results", "lambda_pred"])),
            "sin2_theta13_ext": _num(extract_nested(flavor, ["results", "sin2_theta13"])),
            "sin2_theta12_ext": _num(extract_nested(tm1, ["results", "sin2_theta12"])),
            "delta_cp_deg_ext": _num(extract_nested(flavor, ["results", "delta_cp_deg"])),
            "ckm_lambda_ext": wolf.get("lambda"),
            "ckm_A_ext": wolf.get("A"),
            "ckm_rho_ext": wolf.get("rho"),
            "ckm_eta_ext": wolf.get("eta"),
        }

        ctx_base = _build_context(
            c3=c3_base,
            varphi0=varphi0_base,
            pi=pi_val,
            g_value=g_value,
            k_value=k_value,
            delta_top=delta_top,
            delta2=delta2,
            mp_dps=mp_dps,
            external=external,
        )
        ctx_c3_up = _build_context(
            c3=c3_base * (1 + eps),
            varphi0=varphi0_base,
            pi=pi_val,
            g_value=g_value,
            k_value=k_value,
            delta_top=delta_top,
            delta2=delta2,
            mp_dps=mp_dps,
            external=external,
        )
        ctx_c3_down = _build_context(
            c3=c3_base * (1 - eps),
            varphi0=varphi0_base,
            pi=pi_val,
            g_value=g_value,
            k_value=k_value,
            delta_top=delta_top,
            delta2=delta2,
            mp_dps=mp_dps,
            external=external,
        )
        ctx_phi_up = _build_context(
            c3=c3_base,
            varphi0=varphi0_base * (1 + eps),
            pi=pi_val,
            g_value=g_value,
            k_value=k_value,
            delta_top=delta_top,
            delta2=delta2,
            mp_dps=mp_dps,
            external=external,
        )
        ctx_phi_down = _build_context(
            c3=c3_base,
            varphi0=varphi0_base * (1 - eps),
            pi=pi_val,
            g_value=g_value,
            k_value=k_value,
            delta_top=delta_top,
            delta2=delta2,
            mp_dps=mp_dps,
            external=external,
        )

        symbol_registry_path = Path(__file__).resolve().parent.parent / "symbols.yaml"
        symbol_registry = _load_symbol_registry(symbol_registry_path)
        allowed_symbols = set(symbol_registry.keys()) | {"exp", "sqrt", "log", "ln"}

        def finite(value: float | None) -> bool:
            return value is not None and math.isfinite(value)

        def pick(ctx: dict[str, float], key: str, fallback: str | None = None) -> float | None:
            val = ctx.get(key)
            if not finite(val) and fallback:
                val = ctx.get(fallback)
            return val if finite(val) else None

        def me_value(ctx: dict[str, float]) -> float | None:
            return pick(ctx, "m_e", "m_e_ref")

        def item(
            key: str,
            label: str,
            formula_math: str,
            compute: Callable[[dict[str, float]], float | None],
            *,
            formula_text: str | None = None,
            **kwargs: Any,
        ) -> ConstantItem:
            primitives = kwargs.pop("primitives", ["c3", "varphi0", "pi"])
            return ConstantItem(
                key=key,
                label=label,
                formula_math=formula_math,
                formula_text=formula_text,
                compute=compute,
                primitives=primitives,
                **kwargs,
            )

        groups: dict[str, list[ConstantItem]] = {
            "Group 0: Discrete anchors": [
                item(
                    "pi",
                    "pi",
                    "pi",
                    _ctx_value("pi"),
                    formula_text="circle constant",
                    status="input",
                    lineage="anchor",
                    fixed_point=True,
                    primitives=["pi"],
                    sensitivity_mode="external",
                    units="dimensionless",
                    reference_info={"type": "analytic", "origin": "geometry constant"},
                ),
                item(
                    "g",
                    "g",
                    "g",
                    _ctx_value("g"),
                    formula_text="SU(5) holonomy multiplicity (two-defect)",
                    status="derived",
                    lineage="discrete_scan",
                    source="defect_partition_g5_audit",
                    source_module_id="defect_partition_g5_audit",
                    source_json_pointer="/results/delta2/g_value",
                    sensitivity_mode="external",
                    units="dimensionless",
                    reference_info={"type": "internal_derivation", "origin": "defect_partition_g5_audit"},
                ),
                item(
                    "k",
                    "k",
                    "k",
                    _ctx_value("k"),
                    formula_text="backreaction cover degree (min |ppm|)",
                    status="derived",
                    lineage="discrete_scan",
                    source="alpha_backreaction_sensitivity_audit",
                    source_module_id="alpha_backreaction_sensitivity_audit",
                    source_json_pointer="/results/abs_ppm",
                    sensitivity_mode="external",
                    units="dimensionless",
                    reference_info={"type": "internal_derivation", "origin": "alpha_backreaction_sensitivity_audit"},
                ),
                item(
                    "gamma0",
                    "gamma0",
                    "g/(g+1)",
                    lambda ctx: ctx["g"] / (ctx["g"] + 1.0),
                    formula_text="crosslink signature",
                    depends_on=["g"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                    reference_info={"type": "internal_derivation", "origin": "g5_crosslink_audit"},
                ),
                item(
                    "delta2",
                    "delta2",
                    "(g/4)*delta_top^2",
                    lambda ctx: (ctx["g"] / 4.0) * (ctx["delta_top"] ** 2),
                    depends_on=["g", "delta_top"],
                    status="derived",
                    lineage="discrete_scan",
                    units="dimensionless",
                    reference_info={"type": "internal_derivation", "origin": "defect_partition_g5_audit"},
                ),
            ],
            "Group 1: Kernel invariants": [
                item(
                    "c3",
                    "c3",
                    "1/(8*pi)",
                    _ctx_value("c3"),
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    primitives=["pi"],
                    units="dimensionless",
                ),
                item(
                    "varphi0_tree",
                    "varphi0_tree",
                    "1/(6*pi)",
                    _ctx_value("varphi0_tree"),
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    primitives=["pi"],
                    units="dimensionless",
                ),
                item(
                    "delta_top",
                    "delta_top",
                    "3/(256*pi^4)",
                    _ctx_value("delta_top"),
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    primitives=["pi"],
                    units="dimensionless",
                ),
                item(
                    "varphi0",
                    "varphi0",
                    "varphi0_tree+delta_top",
                    _ctx_value("varphi0"),
                    depends_on=["varphi0_tree", "delta_top"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "delta_star",
                    "delta_star",
                    "3/5+varphi0/6",
                    _ctx_value("delta_star"),
                    depends_on=["varphi0"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "beta_rad",
                    "beta_rad",
                    "varphi0/(4*pi)",
                    _ctx_value("beta_rad"),
                    depends_on=["varphi0", "pi"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="radian",
                ),
            ],
            "Group 2: Gauge sectors (couplings)": [
                item(
                    "alpha_inv0",
                    "alpha_inv(0)",
                    "alpha_inv0",
                    _ctx_value("alpha_inv0"),
                    formula_text=(
                        "Solve CFE(alpha,varphi)=0 with varphi=varphi_tree+delta_top*exp(-2*alpha)"
                        "+delta2*exp(-4*alpha), then alpha_inv0=1/alpha (k=2)."
                    ),
                    reference_id="alpha_inv_codata_2022",
                    inputs=["c3", "varphi0", "delta2"],
                    depends_on=["c3", "varphi0", "delta2"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "alpha0",
                    "alpha(0)",
                    "1/alpha_inv0",
                    lambda ctx: (1.0 / ctx["alpha_inv0"]) if finite(ctx.get("alpha_inv0")) else None,
                    depends_on=["alpha_inv0"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "alpha_inv0_simple",
                    "alpha_inv(0) (simple)",
                    "4*pi*exp(1/varphi0)",
                    lambda ctx: 4 * math.pi * math.exp(1.0 / ctx["varphi0"]),
                    depends_on=["pi", "varphi0"],
                    status="candidate",
                    lineage="candidate",
                    note="Simplified closure candidate.",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "alpha_bar5_inv_MZ",
                    "alpha_bar5_inv(MZ)",
                    "alpha_bar5_inv_MZ",
                    _ctx_value("alpha_bar5_inv_MZ"),
                    formula_text="alpha(0) -> MZ running (alpha_precision_audit)",
                    reference_id="alpha_bar5_inv_MZ_pdg2024",
                    status="derived_external",
                    lineage="derived_external",
                    source="alpha_precision_audit",
                    source_module_id="alpha_precision_audit",
                    source_json_pointer="/results/secondary_alpha_bar5_MZ/alpha_bar5_inv_MZ_pred",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "alpha_s_MZ",
                    "alpha_s(MZ)",
                    "alpha_s_MZ",
                    _const(_dataset_value("alpha_s_mz_pdg")),
                    formula_text="SM input (MSbar)",
                    reference_id="alpha_s_mz_pdg",
                    status="input",
                    lineage="anchor",
                    source="sm_inputs_mz.json",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "sin2_thetaW_MZ",
                    "sin^2(thetaW) (MZ)",
                    "sin2_thetaW_MZ",
                    _const(sin2_thetaW_ref),
                    formula_text="SM input (MSbar)",
                    reference_id="sin2_thetaW_mz_pdg",
                    status="input",
                    lineage="anchor",
                    source="sm_inputs_mz.json",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "sin2_thetaW_0",
                    "sin^2(thetaW)(0)",
                    "sin2_thetaW_0",
                    _const(None),
                    formula_text="G(c3,varphi0) (pending)",
                    status="pending",
                    lineage="pending",
                    units="dimensionless",
                ),
                item(
                    "alpha3_at_mu_star",
                    "alpha3(mu_star)",
                    "alpha3_at_mu_star",
                    _ctx_value("alpha3_at_mu_star"),
                    formula_text="alpha3(mu)=varphi0 crossing",
                    status="derived_external",
                    lineage="derived_external",
                    source="two_loop_rg_fingerprints",
                    source_module_id="two_loop_rg_fingerprints",
                    source_json_pointer="/results/alpha3/crossing_varphi0/alpha3",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "mu_star_varphi0",
                    "mu_star (alpha3=varphi0)",
                    "mu_star_varphi0",
                    _ctx_value("mu_star_varphi0"),
                    formula_text="alpha3(mu)=varphi0 crossing",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="two_loop_rg_fingerprints",
                    source_module_id="two_loop_rg_fingerprints",
                    source_json_pointer="/results/alpha3/crossing_varphi0/mu_GeV",
                    sensitivity_mode="external",
                ),
                item(
                    "g_a_gamma_gamma",
                    "g_a_gamma_gamma",
                    "-4*c3",
                    lambda ctx: -4.0 * ctx["c3"],
                    depends_on=["c3"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
            ],
            "Group 3: Architecture (block scales)": [
                item(
                    "Mpl_reduced",
                    "Mpl_reduced",
                    "Mpl_reduced",
                    _ctx_value("Mpl_reduced"),
                    formula_text="input (reduced Planck mass)",
                    reference_id="Mpl_reduced",
                    status="input",
                    lineage="anchor",
                    units="GeV",
                    sensitivity_mode="external",
                ),
                item(
                    "Mpl_unreduced",
                    "Mpl (unreduced)",
                    "sqrt(8*pi)*Mpl_reduced",
                    _ctx_value("Mpl_unreduced"),
                    depends_on=["pi", "Mpl_reduced"],
                    status="derived",
                    lineage="axiom",
                    units="GeV",
                ),
                item(
                    "M_over_Mpl",
                    "M/Mpl (Starobinsky)",
                    "sqrt(8*pi)*c3^4",
                    _ctx_value("M_over_Mpl"),
                    depends_on=["c3", "pi"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "M_GUT",
                    "M_GUT",
                    "M_GUT",
                    _ctx_value("mu_gut"),
                    formula_text="mu where alpha1=alpha2=alpha3",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="unification_gate",
                    source_module_id="unification_gate",
                    source_json_pointer="/results/best_unification_point/mu_GeV",
                    sensitivity_mode="external",
                ),
                item(
                    "M_R",
                    "M_R",
                    "M_R",
                    _ctx_value("M_R"),
                    formula_text="seesaw anchor",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="seesaw_block",
                    source_module_id="seesaw_block",
                    source_json_pointer="/results/paper_v1_06_anchor/MR_GeV",
                    sensitivity_mode="external",
                ),
                item(
                    "f_a",
                    "f_a",
                    "f_a",
                    _ctx_value("f_a_GeV"),
                    formula_text="axion_dm_pipeline",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="axion_dm_pipeline",
                    source_module_id="axion_dm_pipeline",
                    source_json_pointer="/results/axion_claim/f_a_GeV",
                    sensitivity_mode="external",
                ),
                item(
                    "v_ew",
                    "v (EW scale)",
                    "4*pi*g*c3^2*exp(-alpha_inv0/4)*Mpl_reduced",
                    _ctx_value("v_ew"),
                    formula_text="candidate EW-scale closure",
                    reference_id="v_ew_246",
                    depends_on=["g", "c3", "alpha_inv0", "Mpl_reduced", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    note="Candidate EW-scale formula pending discrete coefficient closure.",
                ),
                item(
                    "Lambda_QCD",
                    "Lambda_QCD",
                    "Mpl_reduced*exp(-1/(2*varphi0))",
                    lambda ctx: float(MPL_REDUCED_GEV) * math.exp(-1.0 / (2.0 * ctx["varphi0"])),
                    formula_text="placeholder scaling (RG-based derivation pending)",
                    depends_on=["Mpl_reduced", "varphi0"],
                    status="pending",
                    lineage="pending",
                    units="GeV",
                    reference_info={"type": "literature", "origin": "PDG Lambda_QCD ~0.2 GeV (numeric bound pending)"},
                    reference_role="pending_reference",
                ),
                item(
                    "m_p",
                    "m_p",
                    "m_p",
                    _ctx_value("m_p_ref"),
                    formula_text="PDG placeholder",
                    reference_id="m_p_pdg",
                    status="placeholder",
                    lineage="anchor",
                    units="GeV",
                    source="mass_spectrum_minimal",
                    sensitivity_mode="external",
                ),
                item(
                    "m_pi",
                    "m_pi",
                    "m_pi",
                    _const(None),
                    formula_text="pending TFPT derivation",
                    status="pending",
                    lineage="pending",
                    units="GeV",
                ),
                item(
                    "G_F_Mpl2",
                    "G_F * Mpl^2",
                    "Mpl_reduced^2/(sqrt(2)*v_ew^2)",
                    lambda ctx: (ctx["Mpl_reduced"] ** 2) / (math.sqrt(2) * (ctx["v_ew"] ** 2)) if finite(ctx.get("v_ew")) else None,
                    depends_on=["Mpl_reduced", "v_ew"],
                    status="candidate",
                    lineage="candidate",
                    units="dimensionless",
                    formula_text="dimensionless Fermi coupling proxy",
                    approximation="scaling_proxy",
                    reference_role="order_of_magnitude_only",
                ),
            ],
            "Group 4: Matter (fermion masses & ratios)": [
                item(
                    "y_tau",
                    "y_tau",
                    "pi*varphi0^2",
                    _ctx_value("y_tau"),
                    depends_on=["pi", "varphi0"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "y_mu",
                    "y_mu",
                    "pi*varphi0^3",
                    _ctx_value("y_mu"),
                    depends_on=["pi", "varphi0"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "y_e",
                    "y_e",
                    "2*pi*varphi0^5",
                    _ctx_value("y_e"),
                    depends_on=["pi", "varphi0"],
                    status="derived",
                    lineage="axiom",
                    fixed_point=True,
                    units="dimensionless",
                ),
                item(
                    "m_e",
                    "m_e",
                    "v_ew*(2*pi*varphi0^5)/sqrt(2)",
                    _ctx_value("m_e"),
                    formula_text="candidate Yukawa ladder",
                    reference_id="m_e_pdg",
                    depends_on=["v_ew", "varphi0", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    note="Depends on EW-scale candidate.",
                ),
                item(
                    "m_mu",
                    "m_mu",
                    "v_ew*(pi*varphi0^3)/sqrt(2)",
                    _ctx_value("m_mu"),
                    formula_text="candidate Yukawa ladder",
                    reference_id="m_mu_pdg",
                    depends_on=["v_ew", "varphi0", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    note="Depends on EW-scale candidate.",
                ),
                item(
                    "m_tau",
                    "m_tau",
                    "v_ew*(pi*varphi0^2)/sqrt(2)",
                    _ctx_value("m_tau"),
                    formula_text="candidate Yukawa ladder",
                    reference_id="m_tau_pdg",
                    depends_on=["v_ew", "varphi0", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    note="Depends on EW-scale candidate.",
                ),
                item(
                    "m_e_over_m_p",
                    "m_e/m_p",
                    "alpha0^2*exp(-1/varphi0)",
                    lambda ctx: (ctx["alpha0"] ** 2) * math.exp(-1.0 / ctx["varphi0"]) if finite(ctx.get("alpha0")) else None,
                    depends_on=["alpha0", "varphi0"],
                    status="candidate",
                    lineage="candidate",
                    note="Candidate ratio from constantfactory.md.",
                    units="dimensionless",
                ),
                item(
                    "m_mu_over_m_e",
                    "m_mu/m_e",
                    "((1+delta_star)/(1-delta_star))^2*((1/3+delta_star)/(1/3-delta_star))^2",
                    lambda ctx: (_mobius_ratio(1.0, ctx["delta_star"]) * abs(_mobius_ratio(1.0 / 3.0, ctx["delta_star"]))) ** 2,
                    reference_id="mass_ratio_mu_over_e_pdg",
                    depends_on=["delta_star"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "m_tau_over_m_mu",
                    "m_tau/m_mu",
                    "((1+delta_star)/(1-delta_star))^2",
                    lambda ctx: (_mobius_ratio(1.0, ctx["delta_star"]) ** 2),
                    reference_id="mass_ratio_tau_over_mu_pdg",
                    depends_on=["delta_star"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "m_t",
                    "m_t (top)",
                    "m_t",
                    _const(m_t_ref),
                    formula_text="input (scheme dependent)",
                    reference_id="m_t_pdg",
                    status="input",
                    lineage="anchor",
                    units="GeV",
                    reference_role="needs_rg",
                    sensitivity_mode="external",
                ),
                item(
                    "m_b",
                    "m_b (bottom)",
                    "m_b",
                    _const(m_b_ref),
                    formula_text="input (scheme dependent)",
                    reference_id="m_b_pdg",
                    status="input",
                    lineage="anchor",
                    units="GeV",
                    reference_role="needs_rg",
                    sensitivity_mode="external",
                ),
                item(
                    "m_c",
                    "m_c (charm)",
                    "m_c",
                    _const(m_c_ref),
                    formula_text="input (scheme dependent)",
                    reference_id="m_c_pdg",
                    status="input",
                    lineage="anchor",
                    units="GeV",
                    reference_role="needs_rg",
                    sensitivity_mode="external",
                ),
                item(
                    "m_s",
                    "m_s (strange)",
                    "m_s",
                    _const(None),
                    formula_text="pending (scheme/scale)",
                    status="pending",
                    lineage="pending",
                    units="GeV",
                ),
                item(
                    "m_u",
                    "m_u (up)",
                    "m_u",
                    _const(None),
                    formula_text="pending (scheme/scale)",
                    status="pending",
                    lineage="pending",
                    units="GeV",
                ),
                item(
                    "m_d",
                    "m_d (down)",
                    "m_d",
                    _const(None),
                    formula_text="pending (scheme/scale)",
                    status="pending",
                    lineage="pending",
                    units="GeV",
                ),
                item(
                    "m_nu1",
                    "m_nu1",
                    "m_nu1",
                    _const(mnu1_placeholder),
                    formula_text="placeholder (order-of-mag)",
                    status="placeholder",
                    lineage="anchor",
                    units="eV",
                    sensitivity_mode="external",
                ),
                item(
                    "m_nu3",
                    "m_nu3",
                    "v_ew^2/M_R",
                    _const(mnu3_eV),
                    depends_on=["v_ew", "M_R"],
                    status="derived_external",
                    lineage="derived_external",
                    units="eV",
                    source="seesaw_block",
                    source_module_id="seesaw_block",
                    source_json_pointer="/results/derived/mnu3_from_v_sm_eV",
                    sensitivity_mode="external",
                ),
                item(
                    "sum_mnu",
                    "sum m_nu",
                    "sum_mnu",
                    _const(None),
                    formula_text="pending (needs consistent neutrino spectrum)",
                    status="pending",
                    lineage="pending",
                    units="eV",
                ),
                item(
                    "delta_m_np",
                    "Delta m_np",
                    "delta_m_np",
                    _const(None),
                    formula_text="pending (EM + isospin)",
                    status="pending",
                    lineage="pending",
                    units="GeV",
                ),
                item(
                    "lambda_C_ratio",
                    "lambda_C^p/lambda_C^e",
                    "m_e/m_p",
                    _const(me_over_mp_ref),
                    depends_on=["m_e", "m_p"],
                    status="derived_external",
                    lineage="derived_external",
                    note="Ratio from PDG masses.",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
            ],
            "Group 5: Flavor code (mixing angles)": [
                item(
                    "lambda",
                    "Cabibbo lambda",
                    "sqrt(varphi0)*(1-varphi0/2)",
                    _ctx_value("lambda_pred"),
                    reference_id="cabibbo_lambda_pdg2024",
                    depends_on=["varphi0"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "CKM_A",
                    "CKM A",
                    "Vcb/lambda^2",
                    _const(wolf.get("A")),
                    formula_text="Wolfenstein from CKM pipeline",
                    status="derived_external",
                    lineage="derived_external",
                    source="ckm_full_pipeline",
                    source_module_id="ckm_full_pipeline",
                    source_json_pointer="/results/alpha_s_sensitivity/matrix_abs_mu_uv/central",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "CKM_rho",
                    "CKM rho",
                    "rho_CKM",
                    _const(wolf.get("rho")),
                    formula_text="Wolfenstein from |Vub|,|Vtd|",
                    status="derived_external",
                    lineage="derived_external",
                    source="ckm_full_pipeline",
                    source_module_id="ckm_full_pipeline",
                    source_json_pointer="/results/alpha_s_sensitivity/matrix_abs_mu_uv/central",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "CKM_eta",
                    "CKM eta",
                    "eta_CKM",
                    _const(wolf.get("eta")),
                    formula_text="Wolfenstein from |Vub|,|Vtd|",
                    status="derived_external",
                    lineage="derived_external",
                    source="ckm_full_pipeline",
                    source_module_id="ckm_full_pipeline",
                    source_json_pointer="/results/alpha_s_sensitivity/matrix_abs_mu_uv/central",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "sin2_theta13",
                    "sin^2(theta13)",
                    "varphi0*exp(-5/6)",
                    _ctx_value("sin2_theta13"),
                    reference_id="pmns_sin2_theta13_nufit53_no",
                    depends_on=["varphi0"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "sin2_theta12",
                    "sin^2(theta12)",
                    "(1/3)*(1-2*sin2_theta13)",
                    _ctx_value("sin2_theta12"),
                    reference_id="pmns_sin2_theta12_nufit53_no",
                    depends_on=["sin2_theta13"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "sin2_theta23",
                    "sin^2(theta23)",
                    "sin2_theta23",
                    _const(sin2_theta23_ref),
                    formula_text="input (NuFIT)",
                    reference_id="pmns_sin2_theta23_nufit53_no",
                    status="input",
                    lineage="anchor",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "delta_cp",
                    "delta_CP (deg)",
                    "180*(1-delta_star)",
                    _ctx_value("delta_cp_deg"),
                    formula_text="placeholder until Z3-breaking is wired",
                    reference_id="pmns_delta_cp_deg_nufit53_no",
                    status="placeholder",
                    lineage="pending",
                    units="degrees",
                    note="Placeholder until Z3-breaking is wired.",
                ),
            ],
            "Group 6: Cosmos (space-time & energy)": [
                item(
                    "phi_star_base",
                    "phi_star_base",
                    "exp(-alpha_inv0/2)",
                    _ctx_value("phi_star_base"),
                    depends_on=["alpha_inv0"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "rho_L",
                    "rho_L",
                    "(1/2*phi_star_base*Mpl_reduced)^4",
                    lambda ctx: (ctx["Mpl_reduced"] * (0.5 * ctx["phi_star_base"])) ** 4 if finite(ctx.get("phi_star_base")) else None,
                    depends_on=["phi_star_base", "Mpl_reduced"],
                    status="derived",
                    lineage="scan",
                    units="GeV^4",
                ),
                item(
                    "Lambda_Mpl2",
                    "Lambda Mpl^-2",
                    "rho_L/Mpl_reduced^4",
                    lambda ctx: ((ctx["Mpl_reduced"] * (0.5 * ctx["phi_star_base"])) ** 4) / (ctx["Mpl_reduced"] ** 4) if finite(ctx.get("phi_star_base")) else None,
                    depends_on=["rho_L", "Mpl_reduced"],
                    status="derived",
                    lineage="scan",
                    units="dimensionless",
                ),
                item(
                    "Omega_b",
                    "Omega_b",
                    "(4*pi-1)*beta_rad",
                    lambda ctx: (4 * math.pi - 1) * ctx["beta_rad"],
                    depends_on=["beta_rad", "pi"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "eta_b",
                    "eta_b",
                    "eta_b_pred",
                    _const(eta_b_pred),
                    formula_text="baryon_consistency_audit proxy",
                    status="derived_external",
                    lineage="derived_external",
                    source="baryon_consistency_audit",
                    source_module_id="baryon_consistency_audit",
                    source_json_pointer="/results/eta_b_pred",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "H0",
                    "H0",
                    "H0_pred",
                    _const(H0_pred),
                    formula_text="baryon_consistency_audit derived",
                    reference_id="H0_planck2018",
                    status="derived_external",
                    lineage="derived_external",
                    units="km/s/Mpc",
                    unit_system="cosmology_units",
                    source="baryon_consistency_audit",
                    source_module_id="baryon_consistency_audit",
                    source_json_pointer="/results/H0_pred",
                    sensitivity_mode="external",
                ),
                item(
                    "Omega_dm",
                    "Omega_dm",
                    "Omega_dm",
                    _const(omega_dm_pred),
                    formula_text="(Omega_a h^2)/h^2",
                    status="derived_external",
                    lineage="derived_external",
                    source="axion_dm_audit",
                    source_module_id="axion_dm_audit",
                    source_json_pointer="/results/relic/Omega_a_h2",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "beta_deg",
                    "beta (deg)",
                    "180/pi*beta_rad",
                    _ctx_value("beta_deg"),
                    reference_id="beta_deg_minami_komatsu_2020",
                    depends_on=["beta_rad", "pi"],
                    status="derived",
                    lineage="axiom",
                    units="degrees",
                    unit_system="cosmology_units",
                ),
                item(
                    "beta_deg_PR",
                    "beta (deg, PR)",
                    "beta_deg",
                    _ctx_value("beta_deg"),
                    formula_text="frequency independent",
                    depends_on=["beta_deg"],
                    status="derived",
                    lineage="axiom",
                    units="degrees",
                    unit_system="cosmology_units",
                ),
                item(
                    "nu_ghz",
                    "axion frequency (GHz)",
                    "nu_ghz",
                    _const(nu_ghz),
                    formula_text="axion_dm_audit pipeline",
                    status="derived_external",
                    lineage="derived_external",
                    source="axion_dm_audit",
                    source_module_id="axion_dm_audit",
                    source_json_pointer="/results/axion_claim/nu_GHz",
                    sensitivity_mode="external",
                    units="GHz",
                    unit_system="cosmology_units",
                ),
                item(
                    "g_phys_axion",
                    "g_phys (axion)",
                    "g_phys",
                    _const(g_phys),
                    formula_text="g_coeff/f_a",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV^-1",
                    source="axion_dm_pipeline",
                    source_module_id="axion_dm_pipeline",
                    source_json_pointer="/results/axion_claim/g_phys_GeV_inv_candidate",
                    sensitivity_mode="external",
                ),
            ],
            "Group 7: Inflation & origin": [
                item(
                    "N_star",
                    "N_star",
                    "N_star",
                    _const(N_star),
                    formula_text="Starobinsky e-folds (policy)",
                    status="input",
                    lineage="anchor",
                    units="dimensionless",
                    reference_info={"type": "internal_policy", "origin": "theoryv3 config starobinsky_N"},
                ),
                item(
                    "A_s",
                    "A_s",
                    "N_star^2/(24*pi^2)*M_over_Mpl^2",
                    _const(A_s_pred),
                    formula_text="Starobinsky amplitude",
                    reference_id="A_s_planck2018",
                    depends_on=["N_star", "M_over_Mpl", "pi"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "n_s",
                    "n_s",
                    "1-2/N_star",
                    _const(n_s_pred),
                    formula_text="Starobinsky tilt",
                    reference_id="n_s_planck2018",
                    depends_on=["N_star"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                ),
                item(
                    "r",
                    "r",
                    "12/N_star^2",
                    _const(r_pred),
                    formula_text="Starobinsky tensor ratio",
                    depends_on=["N_star"],
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                    reference_info={"type": "upper_limit", "origin": "Planck/BICEP (numeric bound pending)"},
                    reference_role="bound",
                ),
                item(
                    "V_star",
                    "V_star",
                    "(3/2)*pi^2*A_s*r*Mpl_reduced^4",
                    _const(V_star_GeV4),
                    depends_on=["A_s", "r", "Mpl_reduced"],
                    status="derived",
                    lineage="axiom",
                    units="GeV^4",
                ),
                item(
                    "bounce_scale",
                    "bounce scale",
                    "k_bounce_t",
                    _const(k_bounce_t),
                    formula_text="bounce_perturbations k_bounce_t (raw)",
                    status="derived_external",
                    lineage="derived_external",
                    source="bounce_perturbations",
                    source_module_id="bounce_perturbations",
                    source_json_pointer="/results/diagnostics/k_bounce_t_est_raw",
                    sensitivity_mode="external",
                    units="dimensionless",
                ),
                item(
                    "x_max",
                    "x_max (Planck spectrum)",
                    "x_max",
                    _const(_planck_x_max()),
                    formula_text="Planck spectrum extremum",
                    status="derived",
                    lineage="axiom",
                    units="dimensionless",
                    reference_info={"type": "analytic", "origin": "Planck spectrum extremum"},
                ),
            ],
            "Group 8: Exotic (bounds & tests)": [
                item(
                    "torsion_max",
                    "max torsion coupling",
                    "torsion_max",
                    _const(torsion_max),
                    formula_text="max inferred |S_mu| bound",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="torsion_bounds_mapping",
                    source_module_id="torsion_bounds_mapping",
                    source_json_pointer="/results/inferred_bounds",
                    sensitivity_mode="external",
                ),
                item(
                    "lambda_cross",
                    "vacuum stability scale",
                    "lambda_cross",
                    _const(lambda_cross),
                    formula_text="lambda crossing scale",
                    status="derived_external",
                    lineage="derived_external",
                    units="GeV",
                    source="stability_unitarity_audit",
                    source_module_id="stability_unitarity_audit",
                    source_json_pointer="/results/lambda_crossing/mu_cross_GeV",
                    sensitivity_mode="external",
                ),
            ],
            "Group 9: QED derived constants": [
                item(
                    "a_e",
                    "a_e",
                    "alpha0/(2*pi)",
                    lambda ctx: (ctx["alpha0"] / (2 * math.pi)) if finite(ctx.get("alpha0")) else None,
                    depends_on=["alpha0", "pi"],
                    status="candidate",
                    lineage="candidate",
                    reference_info={"type": "literature", "origin": "CODATA 2022 QED lowest order"},
                    units="dimensionless",
                ),
                item(
                    "Rydberg_E",
                    "R_infty (energy)",
                    "1/2*alpha0^2*m_e",
                    lambda ctx: 0.5 * (ctx["alpha0"] ** 2) * (me_value(ctx) or 0.0) if finite(ctx.get("alpha0")) and finite(me_value(ctx)) else None,
                    depends_on=["alpha0", "m_e"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    reference_info={"type": "literature", "origin": "CODATA 2022"},
                ),
                item(
                    "a0",
                    "Bohr radius (a0)",
                    "1/(alpha0*m_e)",
                    lambda ctx: (1.0 / (ctx["alpha0"] * (me_value(ctx) or 0.0))) if finite(ctx.get("alpha0")) and finite(me_value(ctx)) else None,
                    depends_on=["alpha0", "m_e"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV^-1",
                    reference_info={"type": "literature", "origin": "CODATA 2022"},
                ),
                item(
                    "sigma_T",
                    "Thomson cross section",
                    "8*pi/3*alpha0^2/m_e^2",
                    lambda ctx: (8 * math.pi / 3.0) * (ctx["alpha0"] ** 2) / ((me_value(ctx) or 0.0) ** 2) if finite(ctx.get("alpha0")) and finite(me_value(ctx)) else None,
                    depends_on=["alpha0", "m_e", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV^-2",
                    reference_info={"type": "literature", "origin": "classical electrodynamics"},
                ),
                item(
                    "fs_2p",
                    "fine-structure 2p (H)",
                    "alpha0^4*m_e",
                    lambda ctx: (ctx["alpha0"] ** 4) * (me_value(ctx) or 0.0) if finite(ctx.get("alpha0")) and finite(me_value(ctx)) else None,
                    depends_on=["alpha0", "m_e"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    approximation="scaling_proxy",
                    reference_role="order_of_magnitude_only",
                ),
                item(
                    "lamb_shift",
                    "Lamb shift (proxy)",
                    "alpha0^5*m_e*log(1/alpha0)",
                    lambda ctx: (ctx["alpha0"] ** 5) * (me_value(ctx) or 0.0) * math.log(1.0 / ctx["alpha0"]) if finite(ctx.get("alpha0")) and finite(me_value(ctx)) else None,
                    depends_on=["alpha0", "m_e"],
                    status="candidate",
                    lineage="candidate",
                    units="GeV",
                    approximation="scaling_proxy",
                    reference_role="order_of_magnitude_only",
                ),
                item(
                    "g_p",
                    "proton g-factor",
                    "2*(1+c3/pi)",
                    lambda ctx: 2.0 * (1.0 + ctx["c3"] / ctx["pi"]),
                    depends_on=["c3", "pi"],
                    status="candidate",
                    lineage="candidate",
                    units="dimensionless",
                    reference_info={"type": "literature", "origin": "classical electrodynamics"},
                ),
            ],
        }

        hierarchy_levels = _derive_hierarchy_levels(
            [
                {"key": item.key, "depends_on": item.depends_on or []}
                for items in groups.values()
                for item in items
            ]
        )

        flat_entries: list[dict[str, Any]] = []
        ledger_entries: list[dict[str, Any]] = []
        references: dict[str, Any] = {}
        checks: list[Check] = []
        status_counts: dict[str, int] = {}
        grammar_violations = 0

        for group_name, items in groups.items():
            for item in items:
                value = item.compute(ctx_base)
                status = item.status
                if value is None or (isinstance(value, float) and not math.isfinite(value)):
                    if status not in {"pending", "placeholder"}:
                        status = "pending"
                status_counts[status] = status_counts.get(status, 0) + 1
                lineage = item.lineage

                ref = None
                comparison: dict[str, Any] | None = None
                deviation: dict[str, Any] | None = None
                delta_ppm = None
                if item.reference_id:
                    try:
                        ref = get_dataset(item.reference_id)
                        references[item.reference_id] = ref
                    except Exception:
                        ref = None

                if item.key == "Omega_b":
                    ref = {
                        "dataset_id": "omega_b_planck_derived",
                        "version": f"{omega_b_h2_ref.get('version')} + {H0_planck_ref.get('version')}",
                        "value": omega_b_ref_value,
                        "sigma": None,
                        "units": "dimensionless",
                        "notes": "Derived Omega_b from Omega_b h^2 and H0 (Planck).",
                    }
                    references["omega_b_planck_derived"] = ref
                if item.key == "eta_b":
                    ref = {
                        "dataset_id": "eta_b_planck_derived",
                        "version": omega_b_h2_ref.get("version"),
                        "value": eta_b_ref_value,
                        "sigma": None,
                        "units": "dimensionless",
                        "notes": "Derived eta_b from Omega_b h^2 (Planck) using eta10=273.9*Omega_b h^2.",
                    }
                    references["eta_b_planck_derived"] = ref
                if item.key == "Omega_dm":
                    ref = {
                        "dataset_id": "omega_dm_planck_derived",
                        "version": f"{omega_dm_h2_ref.get('version')} + {H0_planck_ref.get('version')}",
                        "value": omega_dm_ref_value,
                        "sigma": None,
                        "units": "dimensionless",
                        "notes": "Derived Omega_dm from Omega_dm h^2 and H0 (Planck).",
                    }
                    references["omega_dm_planck_derived"] = ref
                if item.key == "rho_L" and rho_L_target is not None:
                    ref = {
                        "dataset_id": "rho_L_planck_derived",
                        "version": "Planck 2018 (k_calibration)",
                        "value": rho_L_target,
                        "sigma": None,
                        "units": "GeV^4",
                        "notes": "Derived target rho_L from k_calibration cosmology snapshot.",
                    }
                    references["rho_L_planck_derived"] = ref
                if item.key == "Lambda_Mpl2" and Lambda_target is not None:
                    ref = {
                        "dataset_id": "Lambda_planck_derived",
                        "version": "Planck 2018 (k_calibration)",
                        "value": Lambda_target,
                        "sigma": None,
                        "units": "dimensionless",
                        "notes": "Derived Lambda/Mpl^2 from rho_L_target.",
                    }
                    references["Lambda_planck_derived"] = ref
                if item.key == "m_e_over_m_p" and me_over_mp_ref is not None:
                    ref = {
                        "dataset_id": "mass_ratio_me_over_mp",
                        "version": "PDG pole masses",
                        "value": me_over_mp_ref,
                        "sigma": None,
                        "units": "dimensionless",
                        "notes": "Derived from m_e and m_p.",
                    }
                    references["mass_ratio_me_over_mp"] = ref

                if ref and value is not None and ref.get("value") is not None and status not in {"input", "pending", "placeholder"}:
                    ref_value = float(ref["value"])
                    ref_sigma = float(ref["sigma"]) if ref.get("sigma") not in (None, 0) else None
                    if status == "candidate":
                        diff = float(value) - ref_value
                        comparison = {"ref_value": ref_value, "ref_sigma": None, "diff": diff}
                        rel = diff / ref_value if ref_value != 0 else None
                        deviation = {"metric": "rel", "value": rel}
                    elif ref_sigma and ref_sigma > 0:
                        z = (float(value) - ref_value) / ref_sigma
                        comparison = {"ref_value": ref_value, "ref_sigma": ref_sigma, "z": z}
                        deviation = {"metric": "z", "value": z, "abs": abs(z)}
                    else:
                        diff = float(value) - ref_value
                        comparison = {"ref_value": ref_value, "ref_sigma": ref_sigma, "diff": diff}
                        rel = diff / ref_value if ref_value != 0 else None
                        deviation = {"metric": "rel", "value": rel}
                    if ref_value != 0:
                        delta_ppm = 1e6 * (float(value) - ref_value) / ref_value

                mode = item.sensitivity_mode
                if mode is None:
                    if status in {"derived", "candidate"}:
                        mode = "numeric"
                    elif status in {"derived_external"}:
                        mode = "external"
                    else:
                        mode = "external"

                sc3 = None
                sphi0 = None
                use_abs = False
                if item.sensitivity_override:
                    sc3, sphi0 = item.sensitivity_override
                    mode = "analytic"
                elif mode == "numeric":
                    sc3, sphi0, use_abs = _numeric_sensitivity(
                        compute=item.compute,
                        ctx_base=ctx_base,
                        ctx_c3_up=ctx_c3_up,
                        ctx_c3_down=ctx_c3_down,
                        ctx_phi_up=ctx_phi_up,
                        ctx_phi_down=ctx_phi_down,
                        eps=eps,
                    )

                sigma_rel = None
                sigma_abs = None
                if sc3 is not None or sphi0 is not None:
                    sc3_val = float(sc3 or 0.0)
                    sphi0_val = float(sphi0 or 0.0)
                    sigma_rel = math.sqrt((sc3_val * sigma_c3_rel) ** 2 + (sphi0_val * sigma_phi0_rel) ** 2)
                    if value is not None and math.isfinite(float(value)):
                        sigma_abs = abs(float(value)) * sigma_rel

                grammar = {"allowed": True, "issues": [], "fractions": []}
                if item.lineage != "derived_external" and status not in {"input", "pending", "placeholder"}:
                    grammar = _check_formula_grammar(item.formula_math, allowed_symbols)
                    label_text = f"{item.label} {item.formula_text or ''}".lower()
                    if "proxy" in label_text and not item.approximation:
                        grammar["allowed"] = False
                        grammar["issues"].append("proxy_without_flag")
                if not grammar.get("allowed", True) and status not in {"input", "pending", "placeholder"}:
                    grammar_violations += 1

                units = item.units or "dimensionless"
                source_module_id = item.source_module_id or item.source
                source_artifact = item.source_artifact or ("results.json" if source_module_id else None)
                source_json_pointer = item.source_json_pointer

                reference_info = item.reference_info
                if not reference_info and item.lineage == "derived_external" and source_module_id:
                    reference_info = {"type": "pipeline_output", "origin": str(source_module_id)}

                reference_role = item.reference_role
                if reference_role is None:
                    if status == "input":
                        reference_role = "anchor"
                    elif status == "derived_external":
                        reference_role = "pipeline"
                    elif ref is not None:
                        reference_role = "compare"
                    elif reference_info:
                        reference_role = "internal"
                    else:
                        reference_role = "prediction"
                if item.key in {"Omega_b", "eta_b", "Omega_dm"} and ref is not None:
                    reference_role = "derived_from_anchor"

                entry = {
                    "group": group_name,
                    "key": item.key,
                    "label": item.label,
                    "formula_math": item.formula_math,
                    "formula_text": item.formula_text,
                    "value": value,
                    "units": units,
                    "unit_system": item.unit_system,
                    "reference": ref,
                    "reference_info": reference_info,
                    "comparison": comparison,
                    "deviation": deviation,
                    "delta_ppm": delta_ppm,
                    "status": status,
                    "status_label": f"{status}/{lineage}",
                    "lineage": lineage,
                    "fixed_point": item.fixed_point,
                    "note": item.note,
                    "inputs": {k: ctx_base.get(k) for k in (item.inputs or []) if k in ctx_base},
                    "depends_on": item.depends_on or [],
                    "hierarchy_level": hierarchy_levels.get(item.key),
                    "calculation": "",
                    "source": item.source,
                    "source_module_id": source_module_id,
                    "source_artifact": source_artifact,
                    "source_json_pointer": source_json_pointer,
                    "sensitivity": {"mode": mode, "Sc3": sc3, "Sphi0": sphi0, "log_abs": use_abs},
                    "sigma": {"sigma_rel": sigma_rel, "sigma_abs": sigma_abs, "sigma_c3_rel": sigma_c3_rel, "sigma_phi0_rel": sigma_phi0_rel},
                    "grammar": grammar,
                    "latex_ref": item.latex_ref,
                    "approximation": item.approximation,
                    "reference_role": reference_role,
                }
                entry["calculation"] = _calc_text(item, entry)
                flat_entries.append(entry)

                ledger_entry = {
                    "name": item.key,
                    "label": item.label,
                    "definition": item.formula_math,
                    "formula_text": item.formula_text,
                    "primitives": item.primitives or [],
                    "fixed_point": bool(item.fixed_point) if item.fixed_point is not None else False,
                    "sensitivity": {"Sc3": sc3, "Sphi0": sphi0, "mode": mode},
                    "value_pred": value,
                    "value_ref": ref.get("value") if isinstance(ref, dict) else None,
                    "delta_ppm": delta_ppm,
                    "units": units,
                    "unit_system": item.unit_system,
                    "status": status,
                    "lineage": lineage,
                    "proof_ref": item.latex_ref,
                    "source": item.source,
                    "latex_snippet": _latex_snippet(entry),
                    "reference_info": reference_info,
                    "approximation": item.approximation,
                    "reference_role": reference_role,
                    "source_module_id": source_module_id,
                    "source_artifact": source_artifact,
                    "source_json_pointer": source_json_pointer,
                }
                ledger_entries.append(ledger_entry)

        checks.append(mk_check_pass("constant_factory_ran", f"entries={len(flat_entries)}"))
        checks.append(mk_check_info(
            "constant_factory_statuses",
            " ".join([f"{k}={v}" for k, v in sorted(status_counts.items())]),
        ))
        if symbol_registry:
            checks.append(mk_check_info("constant_factory_symbol_registry", f"symbols={len(symbol_registry)}"))
        else:
            checks.append(mk_check_warn("constant_factory_symbol_registry", "registry_empty"))
        if grammar_violations:
            checks.append(mk_check_warn("constant_factory_grammar", f"violations={grammar_violations}"))

        negative_controls = []
        g_neg = extract_nested(defect_partition, ["results", "negative_control"], default=None)
        if g_neg:
            negative_controls.append({"kind": "g", "payload": g_neg})
        k_abs = extract_nested(alpha_sens, ["results", "abs_ppm"], default=None)
        if k_abs:
            negative_controls.append({"kind": "k", "payload": k_abs})
        de_candidates = extract_nested(dark_energy, ["results", "candidates"], default=None)
        if de_candidates:
            negative_controls.append({"kind": "n", "payload": de_candidates})

        crosslinks = []
        g_links = extract_nested(g5_crosslink, ["results"], default=None)
        if isinstance(g_links, dict):
            crosslinks.append({
                "kind": "g_signature",
                "g_value": g_links.get("g_value"),
                "g_over_4": g_links.get("g_over_4"),
                "g_over_2": g_links.get("g_over_2"),
                "g_over_g_plus_1": g_links.get("g_over_g_plus_1"),
            })
        n_from_cover = _num(extract_nested(dark_energy_norm, ["results", "n_from_cover"]))
        best_label = extract_nested(dark_energy_norm, ["results", "best_candidate_label"])
        if n_from_cover is not None:
            crosslinks.append({
                "kind": "k_signature",
                "k": k_value,
                "n_from_cover": n_from_cover,
                "best_label": best_label,
            })

        warnings: list[str] = []
        plot: dict[str, str | None] = {"constant_factory_summary.png": None, "constant_factory_sensitivity.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot_summary, plot_warnings = _plot_constant_hierarchy(out_dir=out_dir, entries=flat_entries)
            plot.update(plot_summary)
            warnings.extend(plot_warnings)
            plot_sens, plot_warnings = _plot_sensitivity_heatmap(out_dir=out_dir, entries=flat_entries)
            plot.update(plot_sens)
            warnings.extend(plot_warnings)

        out_dir = self.output_dir(config)
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "constant_factory_ledger.json").write_text(
                ensure_ascii(json.dumps(ledger_entries, indent=2)), encoding="utf-8"
            )
            (out_dir / "constant_factory_ledger.tex").write_text(
                ensure_ascii("\n".join([entry["latex_snippet"] for entry in ledger_entries])), encoding="utf-8"
            )
        except Exception as exc:
            warnings.append(f"ledger_write_failed: {exc}")

        lines: list[str] = []
        lines.append("Constant factory audit (hierarchical)")
        lines.append("Source: theoryv3/constantfactory.md")
        lines.append(f"Status counts: {status_counts}")
        lines.append(f"Grammar violations: {grammar_violations}")
        lines.append("")
        for group_name, items in groups.items():
            lines.append(group_name)
            lines.append("-" * len(group_name))
            for item in items:
                entry = next(e for e in flat_entries if e["key"] == item.key and e["group"] == group_name)
                lines.append(f"{item.key}: {item.label}")
                lines.append(f"  status: {entry['status_label']}")
                if entry.get("calculation"):
                    lines.append("  calculation:")
                    for line in str(entry["calculation"]).splitlines():
                        lines.append(f"    {line}")
                if entry["reference"]:
                    ref = entry["reference"]
                    lines.append(f"  reference: {ref.get('value')} ± {ref.get('sigma')} ({ref.get('version')})")
                if entry["comparison"]:
                    cmp = entry["comparison"]
                    if "z" in cmp:
                        lines.append(f"  z: {cmp['z']}")
                    else:
                        lines.append(f"  diff: {cmp.get('diff')}")
            lines.append("")
        lines.append("Checks:")
        lines.extend(
            [
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ]
        )

        return ModuleResult(
            results={
                "source_document": "tfpt-suite/theoryv3/constantfactory.md",
                "status_counts": status_counts,
                "policy": {
                    "sensitivity_eps": eps,
                    "sigma_c3_rel": sigma_c3_rel,
                    "sigma_phi0_rel": sigma_phi0_rel,
                    "symbol_registry": {
                        "path": str(symbol_registry_path),
                        "count": len(symbol_registry),
                    },
                },
                "groups": [
                    {
                        "group": group_name,
                        "items": [e for e in flat_entries if e["group"] == group_name],
                    }
                    for group_name in groups.keys()
                ],
                "ledger": ledger_entries,
                "hierarchy": {
                    "levels": hierarchy_levels,
                    "max_level": max(hierarchy_levels.values()) if hierarchy_levels else 0,
                },
                "references": references,
                "views": views,
                "crosslinks": crosslinks,
                "negative_controls": negative_controls,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
