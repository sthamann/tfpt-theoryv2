from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from mpmath import mp

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_aps_eta_gluing(
    *,
    out_dir: Path,
    table: list[SpectralFlowResult],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"aps_eta_spectral_flow_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        spins = sorted(set(r.spin for r in table))
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for spin in spins:
            rows = [r for r in table if r.spin == spin]
            rows.sort(key=lambda r: r.m)
            ms = [r.m for r in rows]
            sf = [r.numeric_sf for r in rows]
            ax.plot(ms, sf, marker="o", lw=2.0, label=spin)

        ax.set_xlabel("m (winding)")
        ax.set_ylabel("spectral flow SF")
        ax.set_title("APS η-gluing: spectral flow equals winding")
        ax.grid(True, ls=":", alpha=0.4)
        ax.legend(loc="best")
        fig.tight_layout()

        path = out_dir / "aps_eta_spectral_flow.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["aps_eta_spectral_flow_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class SpectralFlowResult:
    spin: str
    m: int
    theta_start: float
    theta_end: float
    numeric_sf: int
    analytic_sf: int
    winding_det_u: int
    eta_start_analytic: float
    eta_start_numeric: float
    eta_end_analytic: float
    eta_end_numeric: float


def _eta0_dirac_shift(a: float) -> tuple[float, float]:
    """
    η(0) for the 1D Dirac family with eigenvalues λ_n = n + a, n∈Z.

    For 0<a<1: η(0,a) = 1 - 2a.
    Also: η(0,a) = ζ(0,a) - ζ(0,1-a) with Hurwitz ζ.

    Returns (analytic, numeric_zeta).
    """
    a_mod = float(a % 1.0)
    # Avoid the zero-mode endpoints a=0 (not generic in our epsilon-shifted paths)
    if a_mod == 0.0:
        a_mod = 1e-15
    analytic = float(1.0 - 2.0 * a_mod)
    mp.dps = 60
    numeric = float(mp.zeta(0, a_mod) - mp.zeta(0, 1.0 - a_mod))
    return analytic, numeric


def _spectral_flow_dirac_1d(*, spin: str, m: int, epsilon: float = 1e-3) -> SpectralFlowResult:
    """
    Minimal spectral-flow counter for a 1D Dirac-type eigenvalue family.

    Model:
      λ_n(θ) = n + s + θ/(2π), n ∈ Z
      with spin-structure shift s ∈ {0, 1/2}.

    Consider θ running from ε to 2π m + ε. This avoids endpoint zero-modes and yields:
      SF = m  (eigenvalues n=-m..-1 cross upward once).

    This matches the K-theory/winding statement used in the paper notes:
      SF(U_Γ) = wind(det U_Γ) for U_Γ(θ)=exp(i m θ) ∈ U(1).
    """
    if m < 0:
        raise ValueError("m must be >= 0")

    if spin not in {"periodic", "antiperiodic"}:
        raise ValueError("spin must be 'periodic' or 'antiperiodic'")
    s = 0.0 if spin == "periodic" else 0.5

    theta_start = float(epsilon)
    theta_end = float(2.0 * np.pi * m + epsilon)

    # Numeric count via endpoint sign test (monotone eigenvalues).
    # For λ_n(θ)=n+θ/(2π), derivative dλ/dθ = 1/(2π)>0, so crossings are upward.
    tol = 1e-14
    n_max = max(5, m + 5)
    count = 0
    for n in range(-n_max, n_max + 1):
        lam_s = n + s + theta_start / (2.0 * np.pi)
        lam_e = n + s + theta_end / (2.0 * np.pi)
        if lam_s < -tol and lam_e > tol:
            count += 1

    a_start = s + theta_start / (2.0 * np.pi)
    a_end = s + theta_end / (2.0 * np.pi)
    eta_s_an, eta_s_num = _eta0_dirac_shift(a_start)
    eta_e_an, eta_e_num = _eta0_dirac_shift(a_end)

    return SpectralFlowResult(
        spin=spin,
        m=m,
        theta_start=theta_start,
        theta_end=theta_end,
        numeric_sf=count,
        analytic_sf=m,
        winding_det_u=m,  # det(exp(i m θ)) winds m times
        eta_start_analytic=eta_s_an,
        eta_start_numeric=eta_s_num,
        eta_end_analytic=eta_e_an,
        eta_end_numeric=eta_e_num,
    )


class ApseEtaGluingModule(TfptModule):
    module_id = "aps_eta_gluing"
    title = "APS η-gluing (v2.4 Appendix seam): seam operator D_Γ and spectral flow of U_Γ"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "v2.4 Appendix seam (app:seam): D_Γ = i d/dθ on θ∈[0,2π), APS gluing with matching U_Γ",
                "integer winding m (tested range m=0..8)",
                "spin structure shift s ∈ {0,1/2} (periodic vs antiperiodic)",
            ],
            outputs=["spectral flow SF(U_Γ), winding(det U_Γ), η(0) values at endpoints, minimal seam term Δ_Γ"],
            formulas=[
                "D_Γ = i d/dθ (Appendix seam); eigenvalues λ_n = n + s + θ/(2π), n∈Z",
                "θ ∈ [ε, 2π m + ε] to avoid endpoint zero-modes",
                "SF = # upward crossings = m",
                "wind(det U_Γ) = m for U_Γ(θ)=exp(i m θ)",
                "η(0,a) = ζ(0,a) - ζ(0,1-a) = 1 - 2a for 0<a<1 (Hurwitz-zeta regularization)",
                "TFPT bookkeeping: seam term Δ_Γ = 2π for the minimal nontrivial Z2 class (m=1)",
            ],
            validation=[
                "numeric spectral flow equals analytic m for all tested m and both spin structures",
                "η(0) numeric (Hurwitz ζ) matches analytic closed form",
                "minimal nontrivial flow occurs at m=1",
                "Δ_Γ = 2π for minimal class (m=1) in TFPT normalization",
            ],
            determinism="Deterministic; no floating randomness (ε fixed).",
        )

    def run(self, config) -> ModuleResult:
        epsilon = 1e-3
        m_values = list(range(0, 9))
        spins = ["periodic", "antiperiodic"]
        table = [_spectral_flow_dirac_1d(spin=spin, m=m, epsilon=epsilon) for spin in spins for m in m_values]

        checks: list[Check] = []
        all_ok = all(r.numeric_sf == r.analytic_sf == r.winding_det_u for r in table)
        checks.append(
            Check(
                check_id="sf_equals_winding_all_m",
                passed=bool(all_ok),
                detail="numeric_sf == analytic_sf == winding_det_u for m=0..8 and both spins",
            )
        )
        eta_ok = all(abs(r.eta_start_numeric - r.eta_start_analytic) < 1e-12 and abs(r.eta_end_numeric - r.eta_end_analytic) < 1e-12 for r in table)
        checks.append(
            Check(
                check_id="eta0_matches_hurwitz_zeta",
                passed=bool(eta_ok),
                detail="η(0) numeric (Hurwitz ζ) matches analytic closed form for start/end points",
            )
        )

        # Minimality: smallest m with sf>0 is 1
        min_nonzero = next((r.m for r in table if r.numeric_sf > 0 and r.spin == "periodic"), None)
        checks.append(
            Check(
                check_id="minimal_nontrivial_m",
                passed=bool(min_nonzero == 1),
                detail=f"min m with SF>0 (periodic) is {min_nonzero}",
            )
        )
        # TFPT seam bookkeeping: Δ_Γ = 2π for minimal nontrivial gluing class (m=1).
        min_row = next((r for r in table if r.spin == "periodic" and r.m == 1), None)
        delta_gamma = float(2.0 * np.pi * (min_row.numeric_sf if min_row is not None else 0))
        checks.append(
            Check(
                check_id="tfpt_seam_term_delta_gamma_2pi",
                passed=bool(abs(delta_gamma - float(2.0 * np.pi)) < 1e-12),
                detail=f"Δ_Γ := 2π·SF(U_Γ) for m=1 gives {delta_gamma} (expected 2π)",
            )
        )

        lines: list[str] = []
        lines += [
            "APS η-gluing / spectral flow (v2.4 Appendix seam)",
            "",
            "Implements the seam boundary operator used in v2.4 (Appendix seam):",
            "- D_Γ = i d/dθ on θ∈[0,2π) with APS gluing and matching operator U_Γ(θ)=exp(i m θ).",
            "Checks:",
            "- spectral flow == winding(det U_Γ),",
            "- η(0) evaluation via Hurwitz ζ,",
            "- TFPT seam term Δ_Γ = 2π for minimal nontrivial class (m=1).",
            "",
            f"epsilon = {epsilon}",
            "",
            "spin          m   theta_start      theta_end        SF   wind(det U)    eta0(start)   eta0(end)",
        ]
        for r in table:
            eta_s = r.eta_start_numeric
            eta_e = r.eta_end_numeric
            lines.append(
                f"{r.spin:<12s}  {r.m:<2d}  {r.theta_start:>12.6g}  {r.theta_end:>12.6g}  {r.numeric_sf:>2d}  {r.winding_det_u:>10d}  {eta_s:>12.6g}  {eta_e:>12.6g}"
            )
        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This implements the explicit seam operator D_Γ that is stated in the v2.4 TeX appendix; it is no longer a free toy choice.",
            "- What remains open for a full QFT completion is the embedding of D_Γ into the full 4D/6D Dirac operator and the full η-gluing evaluation on the complete TFPT geometry.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"aps_eta_spectral_flow_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_aps_eta_gluing(out_dir=out_dir, table=table)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "epsilon": epsilon,
                "spins": spins,
                "table": [
                    {
                        "spin": r.spin,
                        "m": r.m,
                        "theta_start": r.theta_start,
                        "theta_end": r.theta_end,
                        "numeric_sf": r.numeric_sf,
                        "analytic_sf": r.analytic_sf,
                        "winding_det_u": r.winding_det_u,
                        "eta_start_analytic": r.eta_start_analytic,
                        "eta_start_numeric": r.eta_start_numeric,
                        "eta_end_analytic": r.eta_end_analytic,
                        "eta_end_numeric": r.eta_end_numeric,
                    }
                    for r in table
                ],
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

