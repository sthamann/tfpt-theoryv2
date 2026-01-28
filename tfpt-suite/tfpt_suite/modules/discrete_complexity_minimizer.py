from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from mpmath import mp
from sympy import Symbol, sqrt as sym_sqrt  # type: ignore
from sympy import simplify as sym_simplify  # type: ignore
from sympy import sympify  # type: ignore
from sympy import pi as sym_pi  # type: ignore

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


@dataclass(frozen=True)
class Expr:
    expr: str
    value: mp.mpf
    cost: int


def _safe_key(x: mp.mpf, digits: int = 30) -> str:
    # key for deduplication (precision-limited string)
    return mp.nstr(x, digits)


def _generate_expression_pool(*, atoms: list[Expr], max_cost: int, magnitude_cap: mp.mpf) -> list[Expr]:
    """
    Exhaustive (but deduplicated) expression generator on a small grammar.

    Grammar:
      constants: user-provided `atoms` (keep this small per target)
      unary: sqrt, pow2, pow3, pow4
      binary: +, -, *, /

    Deduplication by numeric key mp.nstr(value, digits).
    """
    # expressions grouped by cost
    by_cost: dict[int, list[Expr]] = {1: []}
    best: dict[str, Expr] = {}

    def add(e: Expr) -> None:
        if e.cost > max_cost:
            return
        if not mp.isfinite(e.value):
            return
        if abs(e.value) > magnitude_cap:
            return
        k = _safe_key(e.value)
        prev = best.get(k)
        if prev is None or e.cost < prev.cost or (e.cost == prev.cost and e.expr < prev.expr):
            best[k] = e

    for e in atoms:
        add(e)

    by_cost[1] = [e for e in best.values() if e.cost == 1]

    # Dynamic programming by increasing cost
    for new_cost in range(2, max_cost + 1):
        # unary ops from lower costs
        base_cost = new_cost - 1
        for e in by_cost.get(base_cost, []):
            # sqrt
            if e.value >= 0:
                add(Expr(f"sqrt({e.expr})", mp.sqrt(e.value), base_cost + 1))
            # powers
            add(Expr(f"({e.expr})^2", e.value**2, base_cost + 1))
            add(Expr(f"({e.expr})^3", e.value**3, base_cost + 1))
            add(Expr(f"({e.expr})^4", e.value**4, base_cost + 1))

        # binary ops by partitions (ca+cb+1 = new_cost)
        for ca in range(1, new_cost):
            cb = new_cost - ca - 1
            if cb < 1:
                continue
            A = by_cost.get(ca, [])
            B = by_cost.get(cb, [])
            for ea in A:
                for eb in B:
                    # + and * are commutative: enforce ordering to reduce duplicates
                    if ea.expr <= eb.expr:
                        add(Expr(f"({ea.expr})+({eb.expr})", ea.value + eb.value, new_cost))
                        add(Expr(f"({ea.expr})*({eb.expr})", ea.value * eb.value, new_cost))
                    # - and / are non-commutative
                    add(Expr(f"({ea.expr})-({eb.expr})", ea.value - eb.value, new_cost))
                    if eb.value != 0:
                        add(Expr(f"({ea.expr})/({eb.expr})", ea.value / eb.value, new_cost))

        by_cost[new_cost] = [e for e in best.values() if e.cost == new_cost]

    pool = list(best.values())
    pool.sort(key=lambda e: (e.cost, e.expr))
    return pool


def _find_exact_matches(pool: Iterable[Expr], *, target: mp.mpf, rel_tol: mp.mpf) -> list[Expr]:
    matches: list[Expr] = []
    t = target
    scale = abs(t) if t != 0 else mp.mpf(1)
    for e in pool:
        if abs(e.value - t) / scale <= rel_tol:
            matches.append(e)
    matches.sort(key=lambda e: (e.cost, e.expr))
    return matches


def _sympy_env() -> dict[str, object]:
    """
    Symbolic environment matching TFPT constants (paper v2.4).
    """
    pi = sym_pi
    c3 = 1 / (8 * pi)
    delta_top = 3 / (256 * pi**4)
    varphi0 = 1 / (6 * pi) + delta_top
    b1 = 41 / 10
    return {
        "pi": pi,
        "c3": c3,
        "delta_top": delta_top,
        "varphi0": varphi0,
        "b1": b1,
        "sqrt": sym_sqrt,
    }


def _verify_symbolic(expr_str: str, target_expr_str: str) -> bool:
    """
    Verify expr == target by symbolic simplification in a sympy environment.

    Note: This is conservative; expressions with fractional powers are still handled,
    but simplification may be limited. For TFPT constants (positive), principal sqrt
    is the intended branch.
    """
    env = _sympy_env()
    # Our generated strings use '^' for exponent; sympy uses '**'.
    expr_s = expr_str.replace("^", "**")
    target_s = target_expr_str.replace("^", "**")
    try:
        lhs = sympify(expr_s, locals=env)
        rhs = sympify(target_s, locals=env)
        return sym_simplify(lhs - rhs) == 0
    except Exception:
        return False


def _search_min_cost_matches(
    *,
    atoms: list[Expr],
    target: mp.mpf,
    rel_tol: mp.mpf,
    magnitude_cap: mp.mpf,
    max_cost: int,
) -> tuple[int | None, list[Expr], int]:
    """
    Incremental search: stop once minimal-cost matches are found.

    Returns: (min_cost, matches_at_min_cost, explored_unique_count)
    """
    by_cost: dict[int, list[Expr]] = {}
    best: dict[str, Expr] = {}

    def add(e: Expr) -> None:
        if e.cost > max_cost:
            return
        if not mp.isfinite(e.value):
            return
        if abs(e.value) > magnitude_cap:
            return
        k = _safe_key(e.value)
        prev = best.get(k)
        if prev is None or e.cost < prev.cost or (e.cost == prev.cost and e.expr < prev.expr):
            best[k] = e

    for e in atoms:
        add(e)
    by_cost[1] = [e for e in best.values() if e.cost == 1]

    def is_match(val: mp.mpf) -> bool:
        scale = abs(target) if target != 0 else mp.mpf(1)
        return abs(val - target) / scale <= rel_tol

    # Check cost=1
    matches = [e for e in by_cost[1] if is_match(e.value)]
    if matches:
        matches.sort(key=lambda e: e.expr)
        return 1, matches, len(best)

    for new_cost in range(2, max_cost + 1):
        # unary ops: from cost-1 only (keeps branching under control)
        base_cost = new_cost - 1
        for e in by_cost.get(base_cost, []):
            if e.value >= 0:
                add(Expr(f"sqrt({e.expr})", mp.sqrt(e.value), new_cost))
            add(Expr(f"({e.expr})^2", e.value**2, new_cost))
            add(Expr(f"({e.expr})^3", e.value**3, new_cost))
            add(Expr(f"({e.expr})^4", e.value**4, new_cost))

        # binary ops
        for ca in range(1, new_cost):
            cb = new_cost - ca - 1
            if cb < 1:
                continue
            A = by_cost.get(ca, [])
            B = by_cost.get(cb, [])
            for ea in A:
                for eb in B:
                    if ea.expr <= eb.expr:
                        add(Expr(f"({ea.expr})+({eb.expr})", ea.value + eb.value, new_cost))
                        add(Expr(f"({ea.expr})*({eb.expr})", ea.value * eb.value, new_cost))
                    add(Expr(f"({ea.expr})-({eb.expr})", ea.value - eb.value, new_cost))
                    if eb.value != 0:
                        add(Expr(f"({ea.expr})/({eb.expr})", ea.value / eb.value, new_cost))

        by_cost[new_cost] = [e for e in best.values() if e.cost == new_cost]

        matches = [e for e in by_cost[new_cost] if is_match(e.value)]
        if matches:
            # return all matches at this minimal cost (limited)
            matches.sort(key=lambda e: e.expr)
            return new_cost, matches[:50], len(best)

    return None, [], len(best)


class DiscreteComplexityMinimizerModule(TfptModule):
    module_id = "discrete_complexity_minimizer"
    title = "Discrete complexity minimizer (bounded grammar search; rediscover identities)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["TFPT constants + bounded grammar (max_cost=8)"],
            outputs=["minimal-cost expressions matching selected TFPT identities"],
            formulas=[
                "grammar constants: {pi, c3, varphi0, b1, small integers}",
                "grammar ops: +,-,*,/,sqrt,powers",
                "match criterion: relative error <= rel_tol",
                "report minimal-cost exact/near-exact identities (proxy for MDL compression)",
            ],
            validation=[
                "recover beta_rad = varphi0/(4*pi)",
                "recover g_aγγ = -4*c3 (as an equivalent expression if reachable)",
                "recover M/Mpl = sqrt(8*pi)*c3^4",
                "recover Cabibbo λ = sqrt(varphi0)*(1-varphi0/2)",
            ],
            determinism="Deterministic exhaustive search with numeric deduplication.",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        max_cost = 8
        rel_tol = mp.mpf("1e-30")
        magnitude_cap = mp.mpf("1e6")

        # Targets (paper v2.4)
        targets: dict[str, tuple[mp.mpf, list[Expr]]] = {
            "beta_rad": (
                c.beta_rad,
                [
                    Expr("varphi0", c.varphi0, 1),
                    Expr("pi", mp.pi, 1),
                    Expr("4", mp.mpf(4), 1),
                    Expr("1", mp.mpf(1), 1),
                    Expr("2", mp.mpf(2), 1),
                ],
            ),
            "M_over_Mpl": (
                c.M_over_Mpl,
                [
                    Expr("c3", c.c3, 1),
                    Expr("pi", mp.pi, 1),
                    Expr("8", mp.mpf(8), 1),
                    Expr("4", mp.mpf(4), 1),
                    Expr("2", mp.mpf(2), 1),
                    Expr("1", mp.mpf(1), 1),
                ],
            ),
            "cabibbo_lambda": (
                mp.sqrt(c.varphi0) * (mp.mpf(1) - mp.mpf(1) / 2 * c.varphi0),
                [
                    Expr("varphi0", c.varphi0, 1),
                    Expr("1", mp.mpf(1), 1),
                    Expr("2", mp.mpf(2), 1),
                ],
            ),
            "delta_top_over_c3_4": (
                c.delta_top / (c.c3**4),
                [
                    Expr("delta_top", c.delta_top, 1),
                    Expr("c3", c.c3, 1),
                    Expr("1", mp.mpf(1), 1),
                    Expr("2", mp.mpf(2), 1),
                    Expr("3", mp.mpf(3), 1),
                    Expr("4", mp.mpf(4), 1),
                    Expr("6", mp.mpf(6), 1),
                    Expr("8", mp.mpf(8), 1),
                ],
            ),
        }

        matches_summary: dict[str, dict[str, object]] = {}
        checks: list[Check] = []

        explored_total = 0
        symbolic_targets: dict[str, str] = {
            "beta_rad": "(varphi0)/(4*pi)",
            "M_over_Mpl": "sqrt(8*pi)*(c3)^4",
            "cabibbo_lambda": "sqrt(varphi0)*(1 - varphi0/2)",
            "delta_top_over_c3_4": "(delta_top)/(c3)^4",
        }

        for name, (target, atoms) in targets.items():
            min_cost, matches, explored = _search_min_cost_matches(
                atoms=atoms,
                target=target,
                rel_tol=rel_tol,
                magnitude_cap=magnitude_cap,
                max_cost=max_cost,
            )
            explored_total += explored

            target_sym = symbolic_targets.get(name)
            match_rows: list[dict[str, object]] = []
            for m in matches[:10]:
                sym_ok = _verify_symbolic(m.expr, target_sym) if target_sym else False
                match_rows.append({"expr": m.expr, "value": m.value, "cost": m.cost, "symbolic_ok": sym_ok})

            matches_summary[name] = {
                "target": target,
                "min_cost": min_cost,
                "matches": match_rows,
                "explored_unique": explored,
                "symbolic_target": target_sym,
            }

        # sanity checks: we should find something for all targets under this grammar
        checks.append(
            Check(
                check_id="all_targets_matched",
                passed=bool(all(matches_summary[n]["min_cost"] is not None for n in targets.keys())),  # type: ignore[index]
                detail="all selected targets have at least one match within max_cost",
            )
        )
        checks.append(
            Check(
                check_id="symbolic_verification_for_key_identities",
                passed=bool(
                    all(
                        any(m.get("symbolic_ok") for m in matches_summary[name]["matches"])  # type: ignore[index]
                        for name in ("beta_rad", "M_over_Mpl")
                    )
                ),
                detail="symbolic verification succeeded for beta_rad and M_over_Mpl (at least one minimal match each)",
            )
        )

        report_lines = [
            "Discrete complexity minimizer (bounded grammar search)",
            "",
            f"mp.dps = {mp.dps}",
            f"max_cost = {max_cost}",
            f"rel_tol = {rel_tol}",
            f"magnitude_cap = {magnitude_cap}",
            f"explored_unique_total = {explored_total} (sum of unique numeric values explored per target; early-stop at first match cost)",
            "",
            "Targets and minimal-cost matching expressions:",
        ]
        for name, (target, _atoms) in targets.items():
            block = matches_summary[name]
            report_lines.append("")
            report_lines.append(f"- {name}: target={block['target']}")
            report_lines.append(f"  min_cost={block['min_cost']}")
            for m in block["matches"]:
                report_lines.append(f"  * cost={m['cost']}  expr={m['expr']}  value={m['value']}  symbolic_ok={m.get('symbolic_ok')}")

        report_lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This implements a deterministic, bounded search as a proxy for an MDL/description-length minimizer.",
            "- Numeric deduplication is used for performance; matches are additionally verified symbolically where feasible.",
            "- For larger grammars and many observables, the next step would be beam search + holdout + randomization tests (tasks.md).",
            "",
        ]

        return ModuleResult(
            results={
                "settings": {"max_cost": max_cost, "rel_tol": rel_tol, "magnitude_cap": magnitude_cap, "explored_unique_total": explored_total},
                "targets": {k: t[0] for k, t in targets.items()},
                "matches": matches_summary,
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=[],
        )

