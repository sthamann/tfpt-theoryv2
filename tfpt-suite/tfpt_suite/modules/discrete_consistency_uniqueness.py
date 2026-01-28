from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_cfe_uniqueness(
    *,
    out_dir: Path,
    c3: float,
    b1: float,
    varphi0: float,
    alpha_root: float,
    alpha_star: float,
    f_alpha_star: float,
    iter_log: list[dict[str, str]],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"cfe_uniqueness_png": None}
    warnings: list[str] = []

    try:
        import sys

        import matplotlib  # type: ignore
        # Headless-safe backend for CI / non-interactive runs.
        if "matplotlib.pyplot" not in sys.modules:
            matplotlib.use("Agg")  # type: ignore[attr-defined]
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        # Baseline cubic: f(a)=a^3 - 2 c3^3 a^2 - 8 b1 c3^6 ln(1/varphi0)
        C = 8.0 * float(b1) * (float(c3) ** 6) * float(np.log(1.0 / float(varphi0)))
        a_grid = np.linspace(0.0, 0.02, 500)
        f_grid = (a_grid**3) - (2.0 * (float(c3) ** 3) * (a_grid**2)) - float(C)

        # Convergence trace from iter_log
        it = []
        d = []
        for row in iter_log:
            try:
                it.append(int(row.get("iter", "0")))
                d.append(float(row.get("abs_delta", "nan")))
            except Exception:
                continue

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 4.2))

        ax1.plot(a_grid, f_grid, lw=2.0)
        ax1.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        root_pt = ax1.scatter([alpha_root], [0.0], s=50, label="root (baseline)")
        star_pt = ax1.scatter([alpha_star], [f_alpha_star], s=50, label=r"$\alpha_{\ast}$ (local min)")
        ax1.set_xlabel(r"$\alpha$")
        ax1.set_ylabel(r"$f(\alpha)$")
        ax1.set_title("CFE baseline cubic + uniqueness certificate")
        ax1.grid(True, ls=":", alpha=0.4)
        ax1.legend(loc="best")

        if it and d:
            ax2.semilogy(it, d, marker="o", lw=2.0)
            ax2.set_xlabel("fixed-point iteration")
            ax2.set_ylabel(r"$|\Delta\alpha|$")
            ax2.set_title("Self-consistent fixed-point convergence")
            ax2.grid(True, ls=":", alpha=0.4)
        else:
            ax2.text(0.5, 0.5, "No iteration log available", ha="center", va="center")
            ax2.set_axis_off()

        path = out_dir / "cfe_uniqueness.png"

        fig.tight_layout()
        try:
            fig.savefig(path, dpi=160)
        except Exception as e:
            # Safety net: if mathtext parsing fails for any label, fallback to plaintext so
            # tests and artifacts stay reproducible even under stricter mathtext parsers.
            star_pt.set_label("alpha_* (local min)")
            root_pt.set_label("root (baseline)")
            ax1.set_xlabel("alpha")
            ax1.set_ylabel("f(alpha)")
            ax1.set_title("CFE baseline cubic + uniqueness certificate")
            if it and d:
                ax2.set_ylabel("|Delta alpha|")
            ax1.legend(loc="best")
            fig.tight_layout()
            fig.savefig(path, dpi=160)
            warnings.append(f"plot_mathtext_failed_fallback_plaintext: {e}")
        finally:
            plt.close(fig)
        plot["cfe_uniqueness_png"] = str(path)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class CfeRoot:
    alpha: mp.mpf
    alpha_inv: mp.mpf
    varphi0_used: mp.mpf


def _cfe_polynomial(alpha: mp.mpf, *, c3: mp.mpf, varphi0: mp.mpf, b1: mp.mpf) -> mp.mpf:
    # Paper v2.4: α^3 - 2 c3^3 α^2 - 8 b1 c3^6 ln(1/varphi0) = 0
    return alpha**3 - mp.mpf(2) * (c3**3) * alpha**2 - mp.mpf(8) * b1 * (c3**6) * mp.log(mp.mpf(1) / varphi0)


def _solve_cfe_root(*, c3: mp.mpf, varphi0: mp.mpf, b1: mp.mpf) -> CfeRoot:
    """
    Solve the TFPT cubic fixed point equation (CFE) for the unique positive root.
    Uses a robust bracketed secant method via two initial guesses.
    """
    f = lambda a: _cfe_polynomial(a, c3=c3, varphi0=varphi0, b1=b1)

    # The physical root is near alpha ~ 1/137 ~ 0.0073
    x0 = mp.mpf("0.006")
    x1 = mp.mpf("0.010")

    alpha = mp.findroot(f, (x0, x1))
    return CfeRoot(alpha=alpha, alpha_inv=mp.mpf(1) / alpha, varphi0_used=varphi0)


def _unique_positive_root_certificate(*, c3: mp.mpf, varphi0: mp.mpf, b1: mp.mpf) -> dict[str, mp.mpf]:
    """
    A small analytic certificate for uniqueness of the positive cubic root.

    f(α)=α^3 - 2 c3^3 α^2 - C, with C>0.
    f(0)=-C<0 and f→+∞ as α→+∞.

    f'(α)=α(3α-4c3^3) has a local minimum at α* = 4c3^3/3.
    Since f(α*) < 0, f crosses zero exactly once for α>0.
    """
    C = mp.mpf(8) * b1 * (c3**6) * mp.log(mp.mpf(1) / varphi0)
    alpha_star = mp.mpf(4) * (c3**3) / mp.mpf(3)
    f_star = _cfe_polynomial(alpha_star, c3=c3, varphi0=varphi0, b1=b1)
    return {"C": C, "alpha_star": alpha_star, "f(alpha_star)": f_star}


class DiscreteConsistencyUniquenessModule(TfptModule):
    module_id = "discrete_consistency_uniqueness"
    title = "Consistency uniqueness: CFE root (baseline + self-consistent backreaction)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["TFPT invariants (c3, varphi0_tree, delta_top, b1)"],
            outputs=[
                "baseline CFE root (fixed varphi0)",
                "self-consistent CFE root (varphi0(alpha)=varphi_tree+delta_top*exp(-2alpha))",
                "uniqueness certificate (local minimum negative)",
            ],
            formulas=[
                "CFE: α^3 - 2 c3^3 α^2 - 8 b1 c3^6 ln(1/varphi0) = 0",
                "backreaction closure: varphi0(alpha)=varphi_tree + delta_top e^{-2α}",
                "uniqueness (baseline): f(0)<0, f(∞)>0, and f has only one local minimum on α>0 with f(α*)<0",
            ],
            validation=[
                "baseline α^{-1} matches paper baseline (~137.0365)",
                "self-consistent α^{-1} matches paper (~137.03599410...)",
                "iteration converges (difference decreases and reaches tolerance)",
            ],
            determinism="Deterministic given mp.dps (no fitted parameters).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        # Baseline root: fixed varphi0
        baseline = _solve_cfe_root(c3=c.c3, varphi0=c.varphi0, b1=c.b1)
        cert = _unique_positive_root_certificate(c3=c.c3, varphi0=c.varphi0, b1=c.b1)

        # Self-consistent backreaction closure (paper v2.4):
        # varphi0(alpha)=varphi_tree + delta_top e^{-2α}
        def varphi0_of(alpha: mp.mpf) -> mp.mpf:
            return c.varphi0_tree + c.delta_top * mp.e ** (-mp.mpf(2) * alpha)

        # Fixed-point iteration: α_{n+1} = CFE_root(varphi0(alpha_n))
        max_iter = 50
        tol = mp.mpf("1e-30")
        alpha = baseline.alpha
        iter_log: list[dict[str, str]] = []
        converged = False
        for i in range(max_iter):
            v = varphi0_of(alpha)
            nxt = _solve_cfe_root(c3=c.c3, varphi0=v, b1=c.b1).alpha
            diff = abs(nxt - alpha)
            iter_log.append(
                {
                    "iter": str(i),
                    "alpha": str(alpha),
                    "varphi0(alpha)": str(v),
                    "alpha_next": str(nxt),
                    "abs_delta": str(diff),
                }
            )
            if diff < tol:
                alpha = nxt
                converged = True
                break
            alpha = nxt

        self_consistent = CfeRoot(alpha=alpha, alpha_inv=mp.mpf(1) / alpha, varphi0_used=varphi0_of(alpha))

        # Numerical monotonicity evidence for the self-consistent equation:
        # F(α)=α^3 - 2 c3^3 α^2 - 8 b1 c3^6 ln(1/varphi0(α))
        def F(a: mp.mpf) -> mp.mpf:
            return _cfe_polynomial(a, c3=c.c3, varphi0=varphi0_of(a), b1=c.b1)

        a_lo = mp.mpf("0.004")
        a_hi = mp.mpf("0.012")
        grid_n = 200
        # Check strict increase on a coarse grid
        vals = [F(a_lo + (a_hi - a_lo) * mp.mpf(j) / mp.mpf(grid_n)) for j in range(grid_n + 1)]
        mono = all(vals[j + 1] > vals[j] for j in range(grid_n))

        # Paper reference values (v2.4)
        paper_alpha_inv_baseline = mp.mpf("137.03650146")
        paper_alpha_inv_self = mp.mpf("137.03599410")

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="baseline_unique_root_certificate",
                passed=bool(cert["f(alpha_star)"] < 0),
                detail=f"f(alpha*)={cert['f(alpha_star)']} at alpha*={cert['alpha_star']} < 0",
            )
        )
        checks.append(
            Check(
                check_id="baseline_alpha_inv_close",
                passed=bool(abs(baseline.alpha_inv - paper_alpha_inv_baseline) < mp.mpf("5e-6")),
                detail=f"alpha_inv={baseline.alpha_inv} vs paper {paper_alpha_inv_baseline}",
            )
        )
        checks.append(
            Check(
                check_id="self_consistent_converged",
                passed=bool(converged),
                detail=f"converged={converged} in <= {max_iter} iterations, tol={tol}",
            )
        )
        checks.append(
            Check(
                check_id="self_consistent_alpha_inv_close",
                passed=bool(abs(self_consistent.alpha_inv - paper_alpha_inv_self) < mp.mpf("5e-6")),
                detail=f"alpha_inv={self_consistent.alpha_inv} vs paper {paper_alpha_inv_self}",
            )
        )
        checks.append(
            Check(
                check_id="self_consistent_monotone_evidence",
                passed=bool(mono),
                detail=f"F(α) strictly increasing on [{a_lo}, {a_hi}] (grid={grid_n})",
            )
        )

        report_lines = [
            "TFPT consistency uniqueness: CFE root",
            "",
            f"mp.dps = {mp.dps}",
            "",
            "Baseline (fixed varphi0):",
            f"- varphi0 = {baseline.varphi0_used}",
            f"- alpha = {baseline.alpha}",
            f"- alpha^{-1} = {baseline.alpha_inv}",
            "",
            "Uniqueness certificate for baseline cubic:",
            f"- C = 8 b1 c3^6 ln(1/varphi0) = {cert['C']}",
            f"- alpha* = 4 c3^3 / 3 = {cert['alpha_star']}",
            f"- f(alpha*) = {cert['f(alpha_star)']}  (must be < 0 for exactly one positive root)",
            "",
            "Self-consistent backreaction closure (fixed-point iteration):",
            f"- varphi_tree = {c.varphi0_tree}",
            f"- delta_top   = {c.delta_top}",
            f"- varphi0(alpha) = varphi_tree + delta_top * exp(-2 alpha)",
            f"- converged = {converged}",
            f"- alpha = {self_consistent.alpha}",
            f"- alpha^{-1} = {self_consistent.alpha_inv}",
            f"- varphi0_used = {self_consistent.varphi0_used}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"cfe_uniqueness_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_cfe_uniqueness(
                out_dir=out_dir,
                c3=float(c.c3),
                b1=float(c.b1),
                varphi0=float(c.varphi0),
                alpha_root=float(baseline.alpha),
                alpha_star=float(cert["alpha_star"]),
                f_alpha_star=float(cert["f(alpha_star)"]),
                iter_log=iter_log,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "baseline": {"alpha": baseline.alpha, "alpha_inv": baseline.alpha_inv, "varphi0": baseline.varphi0_used},
                "baseline_uniqueness_certificate": cert,
                "self_consistent": {
                    "alpha": self_consistent.alpha,
                    "alpha_inv": self_consistent.alpha_inv,
                    "varphi0": self_consistent.varphi0_used,
                    "converged": converged,
                    "iterations": iter_log,
                },
                "self_consistent_monotonicity": {"a_lo": a_lo, "a_hi": a_hi, "grid_n": grid_n, "strictly_increasing": mono},
                "plot": plot,
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=warnings,
        )

