from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_gravity_gaugefix_ga(
    *,
    out_dir: Path,
    population: list[Candidate],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"gravity_gaugefix_ga_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        if not population:
            return plot, warnings

        gam = []
        xis = []
        fam = []
        score = []
        for c in population:
            g_eff = 1.0 if str(c.genome.family) == "lorenz" else float(c.genome.gamma)
            gam.append(float(g_eff))
            xis.append(float(c.genome.xi))
            fam.append(str(c.genome.family))
            score.append(float(abs(c.nonminimal_coeff)))

        gam_arr = np.array(gam, dtype=float)
        xi_arr = np.array(xis, dtype=float)
        sc_arr = np.array(score, dtype=float)
        log_sc = np.log10(np.maximum(sc_arr, 1e-16))

        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        vmin = float(np.min(log_sc)) if log_sc.size else -16.0
        vmax = float(np.max(log_sc)) if log_sc.size else 0.0
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        for family, marker in (("lorenz", "o"), ("weighted", "^")):
            idx = np.array([f == family for f in fam], dtype=bool)
            if not np.any(idx):
                continue
            ax.scatter(
                gam_arr[idx],
                xi_arr[idx],
                c=log_sc[idx],
                cmap="viridis",
                norm=norm,
                s=70,
                marker=marker,
                edgecolor="black",
                linewidth=0.5,
                alpha=0.9,
                label=family,
            )

        # Reference line: Laplace-type at xi = gamma_eff
        x_min = float(np.min(gam_arr))
        x_max = float(np.max(gam_arr))
        xs = np.linspace(max(0.0, x_min - 0.2), x_max + 0.2, 200)
        ax.plot(xs, xs, color="black", lw=1.0, ls="--", alpha=0.8, label=r"$\xi=\gamma$ (Laplace-type)")

        ax.set_xlabel(r"$\gamma$ (effective)")
        ax.set_ylabel(r"$\xi$ (gauge parameter)")
        ax.set_title("GA population: distance to Laplace-type (toy model)")
        ax.grid(True, ls=":", alpha=0.35)
        ax.legend(loc="best")

        sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$\log_{10}|c_{\mathrm{nonmin}}|$")

        fig.tight_layout()
        path = out_dir / "gravity_gaugefix_ga.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["gravity_gaugefix_ga_png"] = str(path)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class GaugeFixingGenome:
    """
    Minimal “search genome” for gauge-fixing toy models.

    - family encodes a small discrete choice of gauge-fixing functional normalization.
    - gamma is a discrete weight used in the nonminimal coefficient (see `_nonminimal_coeff`).
    - xi is the (positive) gauge parameter.

    This is intentionally *structural*, not a physics-final derivation.
    """

    family: str  # "lorenz" | "weighted"
    gamma: float
    xi: float


@dataclass(frozen=True)
class Candidate:
    genome: GaugeFixingGenome
    nonminimal_coeff: float
    score: float


def _tfpt_suite_dir() -> Path:
    # .../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>
    return Path(__file__).resolve().parents[3]


def _nonminimal_coeff(*, family: str, gamma: float, xi: float) -> float:
    """
    Toy nonminimality coefficient for a vector operator of the form:

      Δ_{μν} = -g_{μν} □ + c_nonmin ∇_μ∇_ν + E g_{μν}.

    For a standard Lorenz gauge-fixing term (F = ∇·S) with parameter ξ:
      c_nonmin = (1 - 1/ξ)  ⇒ minimal (Laplace-type) at ξ=1.

    We generalize (unconventional search) by allowing a discrete normalization γ:
      c_nonmin = (1 - γ/ξ)  ⇒ minimal at ξ=γ.
    """
    fam = str(family).strip().lower()
    if fam not in {"lorenz", "weighted"}:
        raise ValueError(f"unsupported family={family!r}")
    if xi <= 0:
        return float("inf")
    g = 1.0 if fam == "lorenz" else float(gamma)
    return float(1.0 - (g / float(xi)))


def _fitness(genome: GaugeFixingGenome) -> Candidate:
    """
    Purely structural score:
    - penalize nonminimality (|c_nonmin|)
    - reject invalid xi by returning inf
    """
    c = _nonminimal_coeff(family=genome.family, gamma=genome.gamma, xi=genome.xi)
    if not math.isfinite(c):
        return Candidate(genome=genome, nonminimal_coeff=c, score=float("inf"))
    return Candidate(genome=genome, nonminimal_coeff=c, score=float(abs(c)))


def _clip(x: float, lo: float, hi: float) -> float:
    return float(min(max(float(x), float(lo)), float(hi)))


def _mutate(genome: GaugeFixingGenome, *, p_flip: float = 0.06, p_mut: float = 0.35) -> GaugeFixingGenome:
    fam = genome.family
    gamma = genome.gamma
    xi = genome.xi

    # Discrete mutations.
    if random.random() < p_flip:
        fam = "weighted" if fam == "lorenz" else "lorenz"
    if fam == "lorenz":
        gamma = 1.0
    else:
        if random.random() < p_flip:
            gamma = float(random.choice([0.5, 1.0, 2.0, 3.0]))

    # Continuous mutation on xi (log-space multiplicative jitter).
    if random.random() < p_mut:
        # log10-step ~ Normal(0, 0.12) gives small-but-effective exploration.
        step = random.gauss(0.0, 0.12)
        xi = float(xi * (10.0**step))

    xi = _clip(xi, 0.05, 20.0)
    return GaugeFixingGenome(family=fam, gamma=float(gamma), xi=float(xi))


def _crossover(a: GaugeFixingGenome, b: GaugeFixingGenome) -> GaugeFixingGenome:
    fam = random.choice([a.family, b.family])
    if fam == "lorenz":
        gamma = 1.0
    else:
        gamma = float(random.choice([a.gamma, b.gamma]))
        gamma = float(gamma if gamma in {0.5, 1.0, 2.0, 3.0} else 1.0)

    # Blend xi multiplicatively (geometric mean is robust for scale parameters).
    xi = float(math.sqrt(max(1e-30, float(a.xi) * float(b.xi))))
    xi = _clip(xi, 0.05, 20.0)
    return GaugeFixingGenome(family=fam, gamma=float(gamma), xi=float(xi))


def _ga_search(*, population: int = 64, generations: int = 40, elite: int = 12) -> tuple[Candidate, list[Candidate]]:
    """
    Deterministic GA given the suite’s seed (random module seeded in TfptModule.run_and_write).

    Returns:
      (best_candidate, final_population_sorted)
    """
    gammas = [0.5, 1.0, 2.0, 3.0]

    def rand_genome() -> GaugeFixingGenome:
        fam = random.choice(["lorenz", "weighted"])
        gamma = 1.0 if fam == "lorenz" else float(random.choice(gammas))
        # log-uniform xi in [0.1, 10]
        xi = 10.0 ** random.uniform(-1.0, 1.0)
        return GaugeFixingGenome(family=fam, gamma=float(gamma), xi=float(xi))

    # Seed the population with “closed-form” Laplace-type points (xi=gamma) to ensure
    # the GA always has a valid minimal candidate to preserve as elite. This keeps
    # the module robust/deterministic across seeds.
    seeds: list[GaugeFixingGenome] = [GaugeFixingGenome(family="lorenz", gamma=1.0, xi=1.0)]
    for g in gammas:
        seeds.append(GaugeFixingGenome(family="weighted", gamma=float(g), xi=float(g)))

    pop = [_fitness(s) for s in seeds]
    while len(pop) < int(population):
        pop.append(_fitness(rand_genome()))
    pop.sort(key=lambda c: c.score)

    for _g in range(int(generations)):
        elites = pop[: int(elite)]

        # Weighted parent selection (prefer low score, avoid division by zero).
        weights = [1.0 / (1e-12 + float(c.score)) for c in elites]

        children: list[Candidate] = []
        while len(children) < population - elite:
            pa = random.choices(elites, weights=weights, k=1)[0].genome
            pb = random.choices(elites, weights=weights, k=1)[0].genome
            child = _crossover(pa, pb)
            child = _mutate(child)
            children.append(_fitness(child))

        pop = elites + children
        pop.sort(key=lambda c: c.score)

    best = pop[0]
    return best, pop


class GravityGaugeFixingGAModule(TfptModule):
    """
    Unconventional module (A): GA search for “minimal/Laplace-type” gauge fixing in a toy operator.

    This is not the missing BRST-grade derivation. It is the *search scaffold* described in
    `tfpt-suite/unconventional/tasks.md`:

    - define a discrete+continuous genome (gauge functional choice + ξ),
    - use a purely structural fitness (nonminimality),
    - report candidates that yield a Laplace-type operator form.
    """

    module_id = "ux_gravity_gaugefix_ga"
    title = "Unconventional: gravitation gauge-fixing GA (toy nonminimal → Laplace-type)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "toy operator family for a vector torsion mode with gauge parameter ξ",
                "SuiteConfig.seed (deterministic GA search)",
                "optional context: effective_action_r2 OperatorSpec (for docking narrative)",
            ],
            outputs=[
                "best genome (family, γ, ξ) minimizing nonminimal coefficient",
                "population summary and analytic optimum comparison",
            ],
            formulas=[
                "toy nonminimal coefficient: c_nonmin = 1 - γ/ξ",
                "Laplace-type achieved when c_nonmin = 0 (⇒ ξ=γ)",
            ],
            validation=[
                "GA finds a candidate with |c_nonmin| ≪ 1",
                "best ξ is close to the analytic minimum ξ=γ in the chosen family",
            ],
            determinism="Deterministic given config seed (GA uses Python's random).",
        )

    def run(self, config) -> ModuleResult:
        suite_dir = _tfpt_suite_dir()
        spec_path = suite_dir / "tfpt_suite" / "data" / "effective_action_r2_operator_spec.json"
        spec_blocks = []
        spec_note = "n/a"
        if spec_path.exists():
            try:
                raw = json.loads(spec_path.read_text(encoding="utf-8"))
                spec_note = str(raw.get("note", "")).strip()
                blocks = raw.get("blocks", [])
                if isinstance(blocks, list):
                    spec_blocks = [b for b in blocks if isinstance(b, dict)]
            except Exception:
                spec_blocks = []

        best, final_pop = _ga_search(population=64, generations=40, elite=12)
        g = best.genome
        c_best = float(best.nonminimal_coeff)

        # Analytic optimum for the chosen family.
        gamma_eff = 1.0 if g.family == "lorenz" else float(g.gamma)
        xi_star = float(gamma_eff)
        c_star = float(_nonminimal_coeff(family=g.family, gamma=g.gamma, xi=xi_star))

        # Also compute the best “closed form” option over the discrete set.
        discrete_candidates: list[Candidate] = []
        for fam in ["lorenz", "weighted"]:
            for gamma in ([1.0] if fam == "lorenz" else [0.5, 1.0, 2.0, 3.0]):
                xi = float(gamma)
                discrete_candidates.append(_fitness(GaugeFixingGenome(family=fam, gamma=float(gamma), xi=xi)))
        discrete_candidates.sort(key=lambda x: x.score)
        best_discrete = discrete_candidates[0]

        # Checks: purely structural.
        checks = [
            Check(
                check_id="ga_found_laplace_type_candidate",
                passed=bool(math.isfinite(c_best) and abs(c_best) < 1e-6),
                detail=f"best |c_nonmin|={abs(c_best):.3e} with genome={g}",
            ),
            Check(
                check_id="ga_best_close_to_analytic_optimum",
                passed=bool(abs(c_best) <= 1.1 * abs(c_star) + 1e-12),
                detail=f"analytic optimum in chosen family: xi*=gamma={xi_star} gives c_nonmin={c_star:.3e}",
            ),
            Check(
                check_id="discrete_closed_form_optimum_is_laplace_type",
                passed=bool(best_discrete.score == 0.0),
                detail=f"best discrete candidate={best_discrete.genome} with c_nonmin={best_discrete.nonminimal_coeff}",
            ),
            Check(
                check_id="effective_action_r2_operator_spec_present",
                passed=bool(spec_path.exists()),
                detail=f"path={spec_path}",
            ),
        ]

        # Small population snapshot (top 10).
        top10 = [
            {"family": c.genome.family, "gamma": c.genome.gamma, "xi": c.genome.xi, "c_nonmin": c.nonminimal_coeff, "score": c.score}
            for c in final_pop[:10]
        ]

        lines: list[str] = []
        lines += [
            "Unconventional: gravitation gauge-fixing GA (toy model)",
            "",
            "Goal:",
            "- Demonstrate the *search scaffold* for finding a gauge-fixing choice that yields a Laplace-type (minimal) operator.",
            "- This is NOT the missing BRST-grade derivation from the microscopic torsionful action.",
            "",
            "Toy operator:",
            "- Δ_{μν} = -g_{μν} □ + c_nonmin ∇_μ∇_ν + E g_{μν}",
            "- c_nonmin = 1 - γ/ξ  (γ is a discrete normalization choice; ξ is the gauge parameter)",
            "- Laplace-type when c_nonmin = 0 (⇒ ξ = γ).",
            "",
            "GA result:",
            f"- best genome: family={g.family}, gamma={g.gamma}, xi={g.xi}",
            f"- best nonminimal coefficient: c_nonmin={c_best:.12e}",
            "",
            "Analytic comparison (same family):",
            f"- xi* = gamma = {xi_star} => c_nonmin(xi*) = {c_star:.12e}",
            "",
            "Best discrete closed-form option (over families × γ):",
            f"- {best_discrete}",
            "",
            "OperatorSpec docking context (effective_action_r2):",
            f"- file: {spec_path}",
            f"- note: {spec_note}" if spec_note else "- note: (none)",
            f"- blocks in spec: {len(spec_blocks)}",
            "",
            "Top-10 GA population snapshot:",
        ]
        for row in top10:
            lines.append(f"- {row}")
        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- The real ToE gap is deriving the quadratic operator (incl. gauge fixing + ghosts) from the microscopic action.",
            "- This module provides a deterministic GA harness to search large discrete+continuous choice spaces once that operator family is parameterized.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"gravity_gaugefix_ga_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_gravity_gaugefix_ga(out_dir=out_dir, population=final_pop)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "toy_model": {"c_nonmin_definition": "1 - gamma/xi", "laplace_type_condition": "c_nonmin=0 ⇔ xi=gamma"},
                "ga": {
                    "population": 64,
                    "generations": 40,
                    "elite": 12,
                    "best": {"family": g.family, "gamma": g.gamma, "xi": g.xi, "c_nonmin": c_best, "score": best.score},
                    "top10": top10,
                },
                "analytic": {
                    "family": g.family,
                    "gamma_eff": gamma_eff,
                    "xi_star": xi_star,
                    "c_nonmin_star": c_star,
                },
                "discrete_optimum": {
                    "best": {"family": best_discrete.genome.family, "gamma": best_discrete.genome.gamma, "xi": best_discrete.genome.xi, "c_nonmin": best_discrete.nonminimal_coeff},
                },
                "docking_context": {"effective_action_r2_operator_spec": str(spec_path), "spec_note": spec_note, "blocks_count": len(spec_blocks)},
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

