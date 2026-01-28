from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.alpha_running import AlphaRunningInputs, alpha_bar5_MZ_from_alpha0, alpha_running_inputs_from_pdg
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_alpha_k_sensitivity(
    *,
    out_dir: Path,
    k_sensitivity: list[dict[str, Any]],
    k_match: mp.mpf | None,
    alpha_inv_codata_ref: mp.mpf,
    alpha_inv_two_defect_k2: mp.mpf,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"alpha_k_sensitivity_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        ks: list[float] = []
        ppm: list[float] = []
        alpha_inv: list[float] = []
        for row in k_sensitivity:
            try:
                ks.append(float(row["k"]))
                ppm.append(float(row["ppm_vs_codata"]))
                alpha_inv.append(float(row["alpha_inv"]))
            except Exception:
                continue

        if not ks:
            return plot, warnings

        # Sort by k
        order = sorted(range(len(ks)), key=lambda i: ks[i])
        ks = [ks[i] for i in order]
        ppm = [ppm[i] for i in order]
        alpha_inv = [alpha_inv[i] for i in order]

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(9, 6), sharex=True)

        # ppm residuals
        ax1.plot(ks, ppm, marker="o", lw=2.0, label="single-defect (δ2=0)")
        ax1.scatter([2.0], [float((alpha_inv_two_defect_k2 - alpha_inv_codata_ref) / alpha_inv_codata_ref * mp.mpf(1_000_000))], marker="*", s=140, label="two-defect @ k=2 (δ2=δ_top^2)")
        ax1.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        ax1.axvline(2.0, color="black", lw=1.0, ls="--", alpha=0.7, label="k=2 (double cover)")
        if k_match is not None:
            ax1.axvline(float(k_match), color="gray", lw=1.0, ls=":", alpha=0.8, label=f"k_match≈{float(k_match):.4f}")
        ax1.set_ylabel("ppm vs CODATA (α⁻¹(0))")
        ax1.set_title("α backreaction exponent sensitivity")
        ax1.grid(True, ls=":", alpha=0.4)
        ax1.legend(loc="best")

        # alpha_inv values
        ax2.plot(ks, alpha_inv, marker="o", lw=2.0, color="#1f77b4")
        ax2.axhline(float(alpha_inv_codata_ref), color="black", lw=1.0, alpha=0.7, label="CODATA")
        ax2.set_xlabel("k (backreaction exponent)")
        ax2.set_ylabel("α⁻¹(0)")
        ax2.grid(True, ls=":", alpha=0.4)
        ax2.legend(loc="best")

        fig.tight_layout()
        path = out_dir / "alpha_k_sensitivity.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["alpha_k_sensitivity_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


def _plot_alpha_defect_zscore_overview(
    *,
    out_dir: Path,
    alpha_inv_codata_ref: mp.mpf,
    alpha_inv_sigma: mp.mpf,
    alpha_inv_baseline: mp.mpf,
    alpha_inv_k2: mp.mpf,
    alpha_inv_two_defect_half: mp.mpf,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"alpha_defect_zscore_overview_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        if not (alpha_inv_sigma and alpha_inv_sigma > 0):
            return plot, warnings

        variants = [
            ("baseline (fixed $\\varphi_0$)", alpha_inv_baseline),
            ("k=2 (single-defect)", alpha_inv_k2),
            ("k=2 (two-defect, $\\delta_2=\\frac{1}{2}\\delta_{\\mathrm{top}}^2$)", alpha_inv_two_defect_half),
        ]
        labels = [v[0] for v in variants]
        z_scores = [float((v[1] - alpha_inv_codata_ref) / alpha_inv_sigma) for v in variants]

        fig, ax = plt.subplots(figsize=(9, 4.8))
        x = list(range(len(labels)))
        ax.bar(x, z_scores, color="#1f77b4", alpha=0.85)
        ax.axhline(2.0, color="black", lw=1.2, ls="--", alpha=0.7, label="|z|=2")
        ax.axhline(-2.0, color="black", lw=1.2, ls="--", alpha=0.7)
        ax.set_ylabel("z-score vs CODATA (α⁻¹(0))")
        ax.set_title("Defect-expansion variants: α⁻¹(0) z-score overview")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=10, ha="right")
        ax.legend(loc="best")
        fig.tight_layout()

        path = out_dir / "alpha_defect_zscore_overview.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["alpha_defect_zscore_overview_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class AlphaSolveResult:
    alpha: mp.mpf
    alpha_inv: mp.mpf
    iterations: int
    converged: bool


def _cfe(alpha: mp.mpf, *, c3: mp.mpf, b1: mp.mpf, varphi0: mp.mpf) -> mp.mpf:
    """
    CFE polynomial as implemented in global_consistency_test.
    """
    return alpha**3 - mp.mpf(2) * (c3**3) * alpha**2 - mp.mpf(8) * b1 * (c3**6) * mp.log(mp.mpf(1) / varphi0)


def _solve_cfe_for(*, c3: mp.mpf, b1: mp.mpf, varphi0: mp.mpf) -> mp.mpf:
    f = lambda a: _cfe(a, c3=c3, b1=b1, varphi0=varphi0)
    return mp.findroot(f, (mp.mpf("0.006"), mp.mpf("0.010")))


def _fixed_point_alpha(
    *,
    c3: mp.mpf,
    b1: mp.mpf,
    varphi_tree: mp.mpf,
    delta_top: mp.mpf,
    k: mp.mpf,
    delta2: mp.mpf = mp.mpf(0),
    max_iter: int = 60,
    tol: mp.mpf = mp.mpf("1e-30"),
) -> AlphaSolveResult:
    """
    Solve self-consistent alpha by iterating:
      alpha_{n+1} = solve_CFE(varphi(alpha_n))

    with the explicit response model:
      varphi(alpha) = varphi_tree + delta_top * exp(-k alpha) + delta2 * exp(-2k alpha)

    (delta2=0 corresponds to the v2.4 minimal model; delta2 term is a concrete "next correction"
     placeholder whose required magnitude can be inferred from CODATA.)
    """
    alpha = _solve_cfe_for(c3=c3, b1=b1, varphi0=(varphi_tree + delta_top))
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        varphi = varphi_tree + delta_top * mp.e ** (-k * alpha) + delta2 * mp.e ** (-mp.mpf(2) * k * alpha)
        nxt = _solve_cfe_for(c3=c3, b1=b1, varphi0=varphi)
        if abs(nxt - alpha) < tol:
            alpha = nxt
            converged = True
            break
        alpha = nxt
    return AlphaSolveResult(alpha=alpha, alpha_inv=mp.mpf(1) / alpha, iterations=it, converged=converged)


def _ppm(pred: mp.mpf, ref: mp.mpf) -> mp.mpf:
    if ref == 0:
        return mp.mpf(0)
    return (pred - ref) / ref * mp.mpf(1_000_000)


class AlphaPrecisionAuditModule(TfptModule):
    module_id = "alpha_precision_audit"
    title = "α precision audit: canonical self-consistency + defect-expansion diagnostics (assumption-explicit)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3, varphi0_tree, delta_top, b1)",
                "CODATA 2022 alpha_inv reference (from global_reference_minimal.json)",
                "canonical backreaction model (suite/paper v2.5 policy): varphi(alpha)=varphi_tree+delta_top*exp(-k alpha), with k=2 (double cover) as the baseline exponent",
            ],
            outputs=[
                "baseline alpha_inv (fixed varphi0)",
                "self-consistent alpha_inv (k=2)",
                "k-sensitivity table (k grid as in paper)",
                "defect-expansion diagnostics: derived candidate delta2 factors + required delta2 term to match CODATA exactly (target for next derivation)",
                "effective k needed to match CODATA (diagnostic)",
            ],
            formulas=[
                "CFE: alpha^3 - 2 c3^3 alpha^2 - 8 b1 c3^6 ln(1/varphi)=0",
                "backreaction: varphi(alpha)=varphi_tree + delta_top exp(-k alpha)",
                "next correction template: + delta2 exp(-2k alpha) (no claim; used as a debug target)",
                "ppm := 1e6*(pred-ref)/ref",
            ],
            validation=[
                "reproduces the paper’s baseline and self-consistent alpha_inv benchmarks (within tolerance)",
                "produces a stable monotone k-sensitivity table",
            ],
            determinism="Deterministic (root finding + fixed-point iteration).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        # Policy:
        # - Primary comparison observable: ᾱ^(5)(MZ) in MSbar (RG-clean, perturbative).
        # - Secondary/diagnostic: CODATA α(0) in the Thomson limit (on-shell / IR quantity).
        ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference_minimal.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        alpha_inv_codata_ref = mp.mpf(str(ref["observables"]["alpha_inv_codata_2022"]["mean"]))
        alpha_inv_codata_sigma = mp.mpf(str(ref["observables"]["alpha_inv_codata_2022"]["sigma"]))
        alpha_bar5_inv_mz_ref = mp.mpf(str(ref["observables"]["alpha_bar5_inv_MZ"]["mean"]))
        alpha_bar5_inv_mz_sigma = mp.mpf(str(ref["observables"]["alpha_bar5_inv_MZ"]["sigma"]))

        # External SM/QED running inputs for the secondary check.
        pdg_path = Path(__file__).resolve().parent.parent / "data" / "alpha_running_pdg.json"
        pdg = json.loads(pdg_path.read_text(encoding="utf-8"))
        inputs = alpha_running_inputs_from_pdg(pdg=pdg)

        # Baseline: solve CFE with fixed varphi0 := varphi_tree + delta_top (no self-consistency iteration).
        varphi0_fixed = c.varphi0
        alpha_baseline = _solve_cfe_for(c3=c.c3, b1=c.b1, varphi0=varphi0_fixed)
        alpha_inv_baseline = mp.mpf(1) / alpha_baseline

        # Self-consistent: k=2 as in paper double-cover lemma.
        sol_k2 = _fixed_point_alpha(c3=c.c3, b1=c.b1, varphi_tree=c.varphi0_tree, delta_top=c.delta_top, k=mp.mpf(2))

        # Defect-expansion diagnostics:
        #
        # If the single-defect sector contributes δ_top e^{-2α}, then the non-interacting 2-defect sector
        # contributes ~ (1/2!) (δ_top e^{-2α})^2 = (1/2) δ_top^2 e^{-4α}.
        #
        # We report both:
        # - "noninteracting_indistinguishable" candidate: delta2 = (1/2) delta_top^2
        # - legacy "ordered occupancy" candidate:        delta2 = 1 * delta_top^2
        delta2_two_defect_half = mp.mpf("0.5") * (c.delta_top**2)
        sol_k2_two_defect_half = _fixed_point_alpha(
            c3=c.c3,
            b1=c.b1,
            varphi_tree=c.varphi0_tree,
            delta_top=c.delta_top,
            k=mp.mpf(2),
            delta2=delta2_two_defect_half,
        )

        # Legacy parameter-free candidate (kept for backwards compatibility with earlier suite outputs).
        sol_k2_two_defect = _fixed_point_alpha(
            c3=c.c3,
            b1=c.b1,
            varphi_tree=c.varphi0_tree,
            delta_top=c.delta_top,
            k=mp.mpf(2),
            delta2=c.delta_top**2,
        )

        # k-sensitivity grid (paper-style diagnostics)
        k_grid = [mp.mpf(0), mp.mpf(1), mp.mpf("1.5"), mp.mpf(2), mp.mpf("2.5"), mp.mpf(3)]
        table: list[dict[str, Any]] = []
        for k in k_grid:
            sol = _fixed_point_alpha(c3=c.c3, b1=c.b1, varphi_tree=c.varphi0_tree, delta_top=c.delta_top, k=k)
            table.append(
                {
                    "k": str(k),
                    "alpha_inv": sol.alpha_inv,
                    "ppm_vs_codata": _ppm(sol.alpha_inv, alpha_inv_codata_ref),
                    "iterations": sol.iterations,
                    "converged": sol.converged,
                }
            )

        # Diagnostic: find effective k that matches CODATA (no claim; indicates how far "k=2" would need to shift).
        # Use a short bracket search around k=2.
        def alpha_inv_of_k(kv: mp.mpf) -> mp.mpf:
            return _fixed_point_alpha(c3=c.c3, b1=c.b1, varphi_tree=c.varphi0_tree, delta_top=c.delta_top, k=kv).alpha_inv

        def f_k(kv: mp.mpf) -> mp.mpf:
            return alpha_inv_of_k(kv) - alpha_inv_codata_ref

        # bracket near 2 (paper table indicates f(2) slightly negative and f(1.5) positive)
        k_lo = mp.mpf("1.5")
        k_hi = mp.mpf("2.0")
        f_lo = f_k(k_lo)
        f_hi = f_k(k_hi)
        k_match = None
        if f_lo == 0:
            k_match = k_lo
        elif f_hi == 0:
            k_match = k_hi
        elif f_lo * f_hi < 0:
            # bisection for determinism
            for _ in range(80):
                km = (k_lo + k_hi) / 2
                fm = f_k(km)
                if abs(fm) < mp.mpf("1e-30"):
                    k_match = km
                    break
                if f_lo * fm < 0:
                    k_hi = km
                    f_hi = fm
                else:
                    k_lo = km
                    f_lo = fm
            k_match = k_match or (k_lo + k_hi) / 2

        # Debug target: include a delta2 term (exp(-2k alpha)) and solve delta2 required to match CODATA at k=2.
        # This is the minimal "next correction" amplitude needed (to be derived geometrically, not fitted).
        def alpha_inv_of_delta2(delta2: mp.mpf) -> mp.mpf:
            sol = _fixed_point_alpha(
                c3=c.c3,
                b1=c.b1,
                varphi_tree=c.varphi0_tree,
                delta_top=c.delta_top,
                k=mp.mpf(2),
                delta2=delta2,
                max_iter=80,
            )
            return sol.alpha_inv

        def f_d2(d2: mp.mpf) -> mp.mpf:
            return alpha_inv_of_delta2(d2) - alpha_inv_codata_ref

        # Bracket around 0 by trying symmetric magnitudes scaled to delta_top.
        d2_match = None
        d2_scale = c.delta_top
        d2_a = mp.mpf(0)
        f0 = f_d2(d2_a)
        if f0 == 0:
            d2_match = d2_a
        else:
            # try small bracket
            d2_lo = -mp.mpf("1e-3") * d2_scale
            d2_hi = mp.mpf("1e-3") * d2_scale
            f_lo = f_d2(d2_lo)
            f_hi = f_d2(d2_hi)
            # expand if needed
            grow = mp.mpf(10)
            for _ in range(12):
                if f_lo * f_hi < 0:
                    break
                d2_lo *= grow
                d2_hi *= grow
                f_lo = f_d2(d2_lo)
                f_hi = f_d2(d2_hi)
            if f_lo * f_hi < 0:
                for _ in range(100):
                    mid = (d2_lo + d2_hi) / 2
                    fm = f_d2(mid)
                    if abs(fm) < mp.mpf("1e-30"):
                        d2_match = mid
                        break
                    if f_lo * fm < 0:
                        d2_hi = mid
                        f_hi = fm
                    else:
                        d2_lo = mid
                        f_lo = fm
                d2_match = d2_match or (d2_lo + d2_hi) / 2

        # Benchmarks from paper v2.4 table
        paper_baseline = mp.mpf("137.03650146")
        paper_selfcons = mp.mpf("137.03599410")

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="baseline_matches_paper",
                passed=bool(abs(alpha_inv_baseline - paper_baseline) < mp.mpf("5e-7")),
                detail=f"alpha_inv_baseline={alpha_inv_baseline} vs paper={paper_baseline}",
            )
        )
        checks.append(
            Check(
                check_id="selfconsistent_k2_matches_paper",
                passed=bool(abs(sol_k2.alpha_inv - paper_selfcons) < mp.mpf("5e-7")),
                detail=f"alpha_inv_k2={sol_k2.alpha_inv} vs paper={paper_selfcons} (iters={sol_k2.iterations}, converged={sol_k2.converged})",
            )
        )
        checks.append(
            Check(
                check_id="codata_ppm_residual_reproduced",
                passed=bool(abs(_ppm(sol_k2.alpha_inv, alpha_inv_codata_ref) - mp.mpf("-0.037")) < mp.mpf("0.01")),
                detail=f"ppm(k=2 vs CODATA)={_ppm(sol_k2.alpha_inv, alpha_inv_codata_ref)} (paper ~ -0.037 ppm)",
            )
        )
        checks.append(
            Check(
                check_id="codata_primary_reference_present",
                passed=bool(mp.isfinite(alpha_inv_codata_ref) and alpha_inv_codata_ref > 0),
                detail=f"alpha_inv_CODATA(0)={alpha_inv_codata_ref}",
            )
        )
        checks.append(
            Check(
                check_id="alpha_bar5_MZ_primary_reference_present",
                passed=bool(mp.isfinite(alpha_bar5_inv_mz_ref) and alpha_bar5_inv_mz_ref > 0),
                detail=f"alpha_bar5_inv(MZ)={alpha_bar5_inv_mz_ref} ± {alpha_bar5_inv_mz_sigma}",
            )
        )

        # Secondary conversion: alpha(0) -> alpha_bar^(5)(MZ)
        # Use the "non-interacting indistinguishable" two-defect candidate as the most derivation-shaped option.
        alpha0_pred = mp.mpf(1) / sol_k2_two_defect_half.alpha_inv
        alpha_mz_pred, dal = alpha_bar5_MZ_from_alpha0(alpha0=alpha0_pred, inputs=inputs)
        alpha_bar5_inv_mz_pred = mp.mpf(1) / alpha_mz_pred
        alpha_bar5_z = (
            (alpha_bar5_inv_mz_pred - alpha_bar5_inv_mz_ref) / alpha_bar5_inv_mz_sigma
            if alpha_bar5_inv_mz_sigma != 0
            else mp.mpf("nan")
        )
        checks.append(
            Check(
                check_id="alpha_bar5_MZ_within_5sigma",
                passed=bool(mp.isfinite(alpha_bar5_z) and abs(alpha_bar5_z) < mp.mpf(5)),
                detail=f"alpha_bar5_inv(MZ) pred={alpha_bar5_inv_mz_pred}, ref={alpha_bar5_inv_mz_ref} ± {alpha_bar5_inv_mz_sigma} => z={alpha_bar5_z}",
            )
        )

        lines: list[str] = []
        lines += [
            "α precision audit (assumption-explicit)",
            "",
            f"reference file: {ref_path}",
            f"PRIMARY reference: alpha_bar5_inv(MZ) = {alpha_bar5_inv_mz_ref} ± {alpha_bar5_inv_mz_sigma}  (MSbar-at-MZ policy)",
            f"SECONDARY reference (diagnostic): alpha_inv(CODATA 2022, alpha(0)) = {alpha_inv_codata_ref}",
            f"secondary running inputs: {pdg_path}",
            "",
            "Secondary conversion model (external SM/QED running):",
            f"- Δα_had^(5)(MZ) = {inputs.delta_alpha_had5_MZ} (PDG input)",
            f"- Δα_lept(MZ) (1-loop) = {dal['delta_alpha_lept_1loop']}",
            f"- Δα_extra_msbar(MZ) = {dal['delta_alpha_extra_msbar']} (explicit MSbar remainder input)",
            f"- Δα_total(MZ) = {dal['delta_alpha_total']}",
            "",
            "Model 1 (baseline / historical): fixed varphi0 = varphi_tree + delta_top",
            f"- varphi0_fixed = {varphi0_fixed}",
            f"- alpha_inv_baseline = {alpha_inv_baseline}",
            f"- ppm(baseline vs CODATA) = {_ppm(alpha_inv_baseline, alpha_inv_codata_ref)}",
            "",
            "Model 2 (canonical): self-consistent fixed point with k=2 (double cover)",
            f"- alpha_inv(k=2) = {sol_k2.alpha_inv}  (iters={sol_k2.iterations}, converged={sol_k2.converged})",
            f"- ppm(k=2 vs CODATA) = {_ppm(sol_k2.alpha_inv, alpha_inv_codata_ref)}",
            "",
            "Defect-expansion candidates for the next term (no new parameters):",
            f"- Model 3 (non-interacting, indistinguishable): delta2=(1/2)·delta_top^2 => alpha_inv = {sol_k2_two_defect_half.alpha_inv}  (iters={sol_k2_two_defect_half.iterations}, converged={sol_k2_two_defect_half.converged})",
            f"  ppm(Model 3 vs CODATA) = {_ppm(sol_k2_two_defect_half.alpha_inv, alpha_inv_codata_ref)}",
            f"- Model 4 (legacy/ordered occupancy): delta2=delta_top^2 => alpha_inv = {sol_k2_two_defect.alpha_inv}  (iters={sol_k2_two_defect.iterations}, converged={sol_k2_two_defect.converged})",
            f"  ppm(Model 4 vs CODATA) = {_ppm(sol_k2_two_defect.alpha_inv, alpha_inv_codata_ref)}",
            "",
            "Secondary consistency check (TFPT alpha(0) + SM/QED running):",
            f"- alpha_bar5_inv(MZ) pred = {alpha_bar5_inv_mz_pred}",
            f"- z(alpha_bar5_inv(MZ)) = {alpha_bar5_z}",
            "",
            "Exponent sensitivity table (varphi = varphi_tree + delta_top * exp(-k alpha)):",
            "k      alpha_inv             ppm_vs_CODATA      iters  converged",
        ]
        for row in table:
            lines.append(
                f"{row['k']:<4s}  {str(row['alpha_inv']):<20s}  {str(row['ppm_vs_codata']):<16s}  {row['iterations']:<5d}  {row['converged']}"
            )

        lines += [
            "",
            "Diagnostics / next-step targets:",
            f"- effective k to match CODATA (by bisection near 2): {k_match if k_match is not None else 'n/a'}",
            f"- required delta2 (at k=2) to match CODATA exactly: {d2_match if d2_match is not None else 'n/a'}",
        ]
        if k_match is not None:
            lines.append(f"  (gamma_defect required in k = 2 + gamma_defect) = {k_match - mp.mpf(2)}")
        if d2_match is not None:
            lines.append(f"  (delta2 / delta_top) = {d2_match / c.delta_top}")
            # Derived interaction factor target under the (1/2!) combinatorics baseline.
            if c.delta_top != 0:
                C_int_required = d2_match / delta2_two_defect_half
                lines.append(f"  (C_int required in delta2 = (1/2)·C_int·delta_top^2) = {C_int_required}")

        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module does NOT introduce a fit parameter into TFPT. It computes the minimal correction amplitude (delta2) that would be required, which is then a target for a geometric/QFT derivation.",
            "- As parameter-free defect-expansion candidates for the 'next term', it evaluates delta2=(1/2)·delta_top^2 (non-interacting, indistinguishable defects) and delta2=delta_top^2 (legacy/ordered occupancy).",
            "- Under the MSbar-at-MZ policy, the strict comparison observable is ᾱ^(5)(MZ); CODATA α(0) is kept as an IR/on-shell diagnostic.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"alpha_k_sensitivity_png": None, "alpha_defect_zscore_overview_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_alpha_k_sensitivity(
                out_dir=out_dir,
                k_sensitivity=table,
                k_match=k_match,
                alpha_inv_codata_ref=alpha_inv_codata_ref,
                alpha_inv_two_defect_k2=sol_k2_two_defect.alpha_inv,
            )
            plot_defect, plot_defect_warnings = _plot_alpha_defect_zscore_overview(
                out_dir=out_dir,
                alpha_inv_codata_ref=alpha_inv_codata_ref,
                alpha_inv_sigma=alpha_inv_codata_sigma,
                alpha_inv_baseline=alpha_inv_baseline,
                alpha_inv_k2=sol_k2.alpha_inv,
                alpha_inv_two_defect_half=sol_k2_two_defect_half.alpha_inv,
            )
            plot.update(plot_defect)
            warnings.extend(plot_defect_warnings)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference": {
                    "file": str(ref_path),
                    "alpha_inv_codata_2022": alpha_inv_codata_ref,
                    "alpha_bar5_inv_MZ": alpha_bar5_inv_mz_ref,
                    "alpha_bar5_inv_MZ_sigma": alpha_bar5_inv_mz_sigma,
                },
                "secondary_running_inputs": {"file": str(pdg_path)},
                "baseline": {
                    "varphi0_fixed": varphi0_fixed,
                    "alpha_inv": alpha_inv_baseline,
                    "ppm_vs_codata": _ppm(alpha_inv_baseline, alpha_inv_codata_ref),
                },
                "self_consistent": {
                    "k": mp.mpf(2),
                    "alpha_inv": sol_k2.alpha_inv,
                    "ppm_vs_codata": _ppm(sol_k2.alpha_inv, alpha_inv_codata_ref),
                    "iterations": sol_k2.iterations,
                    "converged": sol_k2.converged,
                },
                "self_consistent_two_defect": {
                    "k": mp.mpf(2),
                    "delta2": c.delta_top**2,
                    "delta2_over_delta_top": (c.delta_top**2) / c.delta_top if c.delta_top != 0 else None,
                    "alpha_inv": sol_k2_two_defect.alpha_inv,
                    "ppm_vs_codata": _ppm(sol_k2_two_defect.alpha_inv, alpha_inv_codata_ref),
                    "iterations": sol_k2_two_defect.iterations,
                    "converged": sol_k2_two_defect.converged,
                },
                "self_consistent_two_defect_half": {
                    "k": mp.mpf(2),
                    "delta2": delta2_two_defect_half,
                    "delta2_over_delta_top": (delta2_two_defect_half / c.delta_top) if c.delta_top != 0 else None,
                    "delta2_over_delta_top_sq": (delta2_two_defect_half / (c.delta_top**2)) if c.delta_top != 0 else None,
                    "alpha_inv": sol_k2_two_defect_half.alpha_inv,
                    "ppm_vs_codata": _ppm(sol_k2_two_defect_half.alpha_inv, alpha_inv_codata_ref),
                    "iterations": sol_k2_two_defect_half.iterations,
                    "converged": sol_k2_two_defect_half.converged,
                },
                "secondary_alpha_bar5_MZ": {
                    "alpha_inv_0_used": sol_k2_two_defect_half.alpha_inv,
                    "alpha_inv_0_used_model": "k=2, delta2=(1/2)delta_top^2 (non-interacting two-defect)",
                    "alpha_bar5_inv_MZ_pred": alpha_bar5_inv_mz_pred,
                    "delta_alpha_lept_1loop": dal["delta_alpha_lept_1loop"],
                    "delta_alpha_had5_MZ": inputs.delta_alpha_had5_MZ,
                    "delta_alpha_msbar_on_shell_shift_MZ": dal["delta_alpha_msbar_on_shell_shift"],
                    "delta_alpha_top_decoupling_MZ": dal["delta_alpha_top_decoupling"],
                    "delta_alpha_extra_msbar_MZ": dal["delta_alpha_extra_msbar"],
                    "delta_alpha_extra_total_MZ": dal["delta_alpha_extra_total"],
                    "delta_alpha_total_MZ": dal["delta_alpha_total"],
                },
                "k_sensitivity": table,
                "diagnostics": {
                    "k_match_codata": k_match,
                    "gamma_defect_required": (k_match - mp.mpf(2)) if k_match is not None else None,
                    "delta2_match_codata_at_k2": d2_match,
                    "delta2_over_delta_top": (d2_match / c.delta_top) if (d2_match is not None and c.delta_top != 0) else None,
                    "C_int_required_in_delta2_half_factor": (
                        (d2_match / delta2_two_defect_half) if (d2_match is not None and delta2_two_defect_half != 0) else None
                    ),
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

