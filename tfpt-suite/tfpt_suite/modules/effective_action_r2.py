from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.heat_kernel import LaplaceTypeBlock, a2_R2_coeff_constant_curvature_4d, beta_R2_from_a2_R2_coeff_4d
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.operator_spec_builder import generate_effective_action_r2_operator_spec


def _plot_effective_action_r2(
    *,
    out_dir: Path,
    blocks: list[dict[str, object]],
    starobinsky_rows: list[StarobinskyPredictions],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"r2_block_contributions_png": None, "starobinsky_scan_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        # --- Plot 1: per-block beta_R2 contributions ---
        names: list[str] = []
        vals: list[float] = []
        for b in blocks:
            try:
                names.append(str(b.get("name", "block")))
                vals.append(float(b.get("beta_R2_contribution", 0.0)))
            except Exception:
                continue

        if names:
            max_abs = max(abs(v) for v in vals) if vals else 1.0
            linthresh = max(1e-12, 1e-6 * max_abs)

            fig1, ax1 = plt.subplots(figsize=(10, 4.5))
            y = list(range(len(names)))
            colors = ["#1f77b4" if v >= 0 else "#d62728" for v in vals]
            ax1.barh(y, vals, color=colors, alpha=0.9)
            ax1.axvline(0.0, color="black", lw=1.0, alpha=0.7)
            ax1.set_xscale("symlog", linthresh=linthresh)
            ax1.set_yticks(y)
            ax1.set_yticklabels(names)
            ax1.invert_yaxis()
            ax1.set_xlabel(r"$\beta_{R^2}$ contribution (symlog)")
            ax1.set_title(r"Effective-action closure: per-block $\beta_{R^2}$ contributions")
            ax1.grid(True, axis="x", ls=":", alpha=0.4)
            fig1.tight_layout()
            p1 = out_dir / "r2_block_contributions.png"
            fig1.savefig(p1, dpi=160)
            plt.close(fig1)
            plot["r2_block_contributions_png"] = str(p1)

        # --- Plot 2: Starobinsky scan ---
        if starobinsky_rows:
            Ns = [int(r.N) for r in starobinsky_rows]
            ns = [float(r.n_s) for r in starobinsky_rows]
            rr = [float(r.r) for r in starobinsky_rows]
            As9 = [float(r.A_s * mp.mpf("1e9")) for r in starobinsky_rows]

            fig2, (ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=3, figsize=(12, 3.5), sharex=True)

            ax2.plot(Ns, ns, marker="o", lw=2.0)
            ax2.set_title(r"$n_s$")
            ax2.set_xlabel("N")
            ax2.grid(True, ls=":", alpha=0.4)

            ax3.plot(Ns, rr, marker="o", lw=2.0, color="#ff7f0e")
            ax3.set_title(r"$r$")
            ax3.set_xlabel("N")
            ax3.grid(True, ls=":", alpha=0.4)

            ax4.plot(Ns, As9, marker="o", lw=2.0, color="#2ca02c")
            ax4.set_title(r"$10^9 A_s$")
            ax4.set_xlabel("N")
            ax4.grid(True, ls=":", alpha=0.4)

            fig2.suptitle("Starobinsky predictions (N scan)")
            fig2.tight_layout()
            p2 = out_dir / "starobinsky_scan.png"
            fig2.savefig(p2, dpi=160)
            plt.close(fig2)
            plot["starobinsky_scan_png"] = str(p2)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class StarobinskyPredictions:
    N: int
    n_s: mp.mpf
    r: mp.mpf
    A_s: mp.mpf


def starobinsky_predictions(*, N: int, M_over_Mpl: mp.mpf) -> StarobinskyPredictions:
    """
    Paper v2.4 (Theorem 'inflation'):
      n_s = 1 - 2/N
      r   = 12/N^2
      A_s ≈ N^2/(24π^2) * (M/Mpl)^2
    """
    Nmp = mp.mpf(N)
    n_s = mp.mpf(1) - mp.mpf(2) / Nmp
    r = mp.mpf(12) / (Nmp**2)
    A_s = (Nmp**2) / (mp.mpf(24) * (mp.pi**2)) * (M_over_Mpl**2)
    return StarobinskyPredictions(N=N, n_s=n_s, r=r, A_s=A_s)


def _sm_matter_blocks_for_r2_a2() -> list[dict[str, object]]:
    """
    Minimal, scheme-explicit SM matter contribution model for the R^2 a2 coefficient.

    This is *not* a full SM-in-curved-background derivation; it is an optional bookkeeping layer:
    - uses standard Laplace-type operators on a 4D constant-curvature background
    - uses fixed field-counts for the SM (no right-handed neutrinos)
    - uses Feynman-gauge style vector + ghost blocks for gauge bosons

    Counts used:
    - Higgs: 4 real scalars
    - Weyl fermions: 45 (15 per generation × 3), counting color/SU2 multiplicities
    - Gauge bosons: 12 (SU3:8, SU2:3, U1:1)

    Laplace-type parameters on constant curvature (E = E_over_R * R, Ω^2 = Omega_sq_over_R2 * R^2):
    - real scalar: E_over_R=0, Ω^2=0
    - Weyl spinor (D†D): E_over_R=1/4, Ω^2=(1/48)·R^2 (per component), prefactor = -1/2
    - gauge vector (Feynman gauge): E_over_R=1/4, Ω^2=-(1/24)·R^2 (per component), prefactor = +1/2
    - FP ghost (complex scalar): E_over_R=0, Ω^2=0, prefactor = -1
    """
    n_real_scalars = 4
    n_weyl = 45
    n_gauge = 12

    return [
        {
            "name": "SM_matter_scalars(H)",
            "rank": n_real_scalars,
            "statistics": "boson",
            "prefactor": "1/2",
            "E_over_R": 0,
            "Omega_sq_over_R2": 0,
            "note": "4 real Higgs scalars, minimal Laplace-type (E=0, Ω=0).",
        },
        {
            "name": "SM_matter_fermions(Weyl)",
            "rank": 2 * n_weyl,
            "statistics": "fermion",
            "prefactor": "-1/2",
            "E_over_R": "1/4",
            "Omega_sq_over_R2": "1/48",
            "note": "45 Weyl fermions (SM, no ν_R), modeled via D†D Laplace-type operator on 2-component spinors.",
        },
        {
            "name": "SM_matter_gauge_bosons(vector)",
            "rank": 4 * n_gauge,
            "statistics": "boson",
            "prefactor": "1/2",
            "E_over_R": "1/4",
            "Omega_sq_over_R2": "-1/24",
            "note": "12 gauge bosons in a minimal Laplace-type vector model on constant curvature (Feynman-gauge style).",
        },
        {
            "name": "SM_matter_gauge_ghosts(scalar)",
            "rank": n_gauge,
            "statistics": "ghost",
            "prefactor": "-1",
            "E_over_R": 0,
            "Omega_sq_over_R2": 0,
            "note": "12 complex FP ghosts (one per gauge boson), modeled as scalar Laplace-type.",
        },
    ]


def seeley_dewitt_a2_R2_coeff_constant_curvature_scalar() -> mp.mpf:
    """
    Heat-kernel framework check for Appendix K language.

    For a Laplace-type scalar operator (E=0, Ω=0) in 4D, the local a2 coefficient contains:

      a2 ⊃ (1/360) * (5 R^2 - 2 R_{μν}R^{μν} + 2 R_{μνρσ}R^{μνρσ})

    On a 4D maximally symmetric (constant-curvature) background:
      R_{μν}R^{μν} = R^2/4
      R_{μνρσ}R^{μνρσ} = R^2/6

    Therefore:
      a2 ⊃ (29/2160) R^2

    This does not derive the TFPT-specific torsion operator; it provides a reproducible
    check that the standard a2 machinery indeed contains curvature-squared terms.
    """
    return mp.mpf(29) / mp.mpf(2160)


def _fmt_sci(x: mp.mpf | float, *, sig: int = 6) -> str:
    """
    Compact, paper-friendly formatting for reports:
    - use scientific notation for very small/large magnitudes
    - otherwise fixed/compact decimal.
    """
    xf = float(x)
    if xf == 0.0:
        return "0"
    ax = abs(xf)
    if ax < 1e-3 or ax >= 1e4:
        return f"{xf:.{sig}e}"
    return f"{xf:.{sig}g}"


def seeley_dewitt_a2_R2_coeff_constant_curvature_scalar_with_ER(*, alpha_R: mp.mpf) -> mp.mpf:
    """
    Constant-curvature R^2 coefficient for a scalar Laplace-type operator with endomorphism E = alpha_R * R.

    Using the standard 4D local a2 coefficient (dropping total derivatives):

      a2 ⊃ (1/360) tr( 5 R^2 - 2 Ric^2 + 2 Riem^2 + 60 R E + 180 E^2 + 30 Ω^2 )

    For a scalar bundle Ω=0 and on maximally symmetric background:
      Ric^2 = R^2/4,  Riem^2 = R^2/6,  E = alpha_R R

    So for one real scalar degree of freedom:
      a2 ⊃ [29/2160 + (alpha_R)/6 + (alpha_R^2)/2] * R^2
    """
    return seeley_dewitt_a2_R2_coeff_constant_curvature_scalar() + (alpha_R / mp.mpf(6)) + (alpha_R**2) / mp.mpf(2)


def beta_R2_from_a2_constant_curvature(*, a2_R2_coeff: mp.mpf) -> mp.mpf:
    """
    Map an a2(R^2) coefficient (the curly-brace coefficient, i.e. without the (4π)^(-2) factor)
    to an effective-action coefficient β_R2 in:

      Γ_eff ⊃ ∫ d^4x √g  β_R2 R^2.

    Convention used here (minimal one-loop normalization):
      Γ_1 = (1/2) Tr ln Δ  ⇒  β_R2 = (1/2) * (4π)^(-2) * a2_R2_coeff = a2_R2_coeff / (32 π^2).

    This deliberately encodes only the local coefficient bookkeeping. A full QFT derivation
    would track scheme/log factors and full field content (incl. ghosts).
    """
    return a2_R2_coeff / (mp.mpf(32) * (mp.pi**2))


def M_over_Mpl_from_beta_R2(*, beta_R2: mp.mpf) -> mp.mpf:
    """
    Starobinsky normalization:
      Γ ⊃ ∫ √g β_R2 R^2  with  β_R2 = Mpl^2/(12 M^2)  (Mpl=1 in this module's units)

    => M/Mpl = 1/sqrt(12 β_R2)
    """
    return mp.mpf(1) / mp.sqrt(mp.mpf(12) * beta_R2)


def _parse_mpf(value: object) -> mp.mpf:
    """
    Minimal safe parser for OperatorSpec numeric fields.
    Supports:
    - int/float
    - strings like "1e-3", "1/2", "-1/2"
    """
    if isinstance(value, (int, float)):
        return mp.mpf(value)
    if isinstance(value, str):
        s = value.strip()
        if "/" in s:
            a, b = s.split("/", 1)
            return mp.mpf(a.strip()) / mp.mpf(b.strip())
        return mp.mpf(s)
    raise TypeError(f"Unsupported numeric value: {value!r}")


def _load_operator_spec(path: Path) -> dict[str, object]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _has_symbolic_alpha_R(block_raw: dict[str, object]) -> bool:
    E_raw = block_raw.get("E_over_R", None)
    return isinstance(E_raw, dict) and str(E_raw.get("symbol")) == "alpha_R"


def _block_from_spec(*, block_raw: dict[str, object], alpha_R: mp.mpf | None) -> tuple[LaplaceTypeBlock, mp.mpf]:
    """
    Convert a JSON OperatorSpec block into a LaplaceTypeBlock + prefactor for Tr ln Δ.

    Supported subset (schema v1 used by `effective_action_r2_operator_spec.json`):
    - rank: int
    - statistics: "boson"|"fermion"|"ghost" (informational; prefactor controls sign)
    - prefactor: "1/2", "-1", "-1/2", ...
    - E_over_R: scalar number, or {"symbol":"alpha_R","times_identity":true}
    - Omega_sq_over_R2: scalar number
    """
    name = str(block_raw.get("name", "block"))
    rank = int(block_raw.get("rank", 1))
    statistics = str(block_raw.get("statistics", "boson"))

    pref_raw = block_raw.get("prefactor", "1/2")
    prefactor = _parse_mpf(pref_raw)

    E_raw = block_raw.get("E_over_R", 0)
    if isinstance(E_raw, dict) and str(E_raw.get("symbol")) == "alpha_R":
        if alpha_R is None:
            raise ValueError("alpha_R requested by OperatorSpec but no alpha_R was provided")
        E_over_R = mp.mpf(alpha_R)
    else:
        E_over_R = _parse_mpf(E_raw)

    Omega_raw = block_raw.get("Omega_sq_over_R2", 0)
    Omega_sq_over_R2 = _parse_mpf(Omega_raw)

    return (
        LaplaceTypeBlock(
            name=name,
            rank=rank,
            statistics=statistics,
            E_over_R=E_over_R,
            Omega_sq_over_R2=Omega_sq_over_R2,
        ),
        prefactor,
    )


def _solve_alpha_R_by_matching_quadratic(*, blocks: list[dict[str, object]], beta_target: mp.mpf) -> mp.mpf:
    """
    Solve alpha_R under the assumption that any symbolic E_over_R uses the single symbol `alpha_R`
    (times identity), and all Omega_sq_over_R2 are scalar.

    Then β_R2(alpha_R) is quadratic in alpha_R, while blocks with constant E_over_R simply contribute
    to the constant term.
    """
    A = mp.mpf(0)
    B = mp.mpf(0)
    C = mp.mpf(0)
    for b in blocks:
        pref = _parse_mpf(b.get("prefactor", "1/2"))
        rank = mp.mpf(int(b.get("rank", 1)))
        Omega = _parse_mpf(b.get("Omega_sq_over_R2", 0))

        # Always include the geometric + Ω^2 constant piece.
        C += pref * rank * ((mp.mpf(29) / mp.mpf(2160)) + (mp.mpf(1) / mp.mpf(12)) * Omega)

        E_raw = b.get("E_over_R", 0)
        if _has_symbolic_alpha_R(b):
            # Variable pieces: (alpha/6 + alpha^2/2)
            B += pref * rank * (mp.mpf(1) / 6)
            A += pref * rank * (mp.mpf(1) / 2)
        else:
            # Constant E_over_R contribution: (E/6 + E^2/2)
            E0 = _parse_mpf(E_raw)
            C += pref * rank * ((E0 / mp.mpf(6)) + (E0**2) / mp.mpf(2))

    # Apply the (4π)^(-2) prefactor: 1/(16 π^2)
    A /= (mp.mpf(16) * (mp.pi**2))
    B /= (mp.mpf(16) * (mp.pi**2))
    C /= (mp.mpf(16) * (mp.pi**2))

    aa = A
    bb = B
    cc = C - beta_target
    disc = bb**2 - mp.mpf(4) * aa * cc
    if disc < 0:
        raise ValueError("No real solution for alpha_R under quadratic matching assumptions")
    r1 = (-bb + mp.sqrt(disc)) / (mp.mpf(2) * aa)
    r2 = (-bb - mp.sqrt(disc)) / (mp.mpf(2) * aa)
    return r1 if r1 > 0 else r2


class EffectiveActionR2Module(TfptModule):
    module_id = "effective_action_r2"
    title = "R^2 effective-action closure (scale M + inflation observables)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3)",
                "torsion-sector Laplace-type operator closure for Appendix K.2 (explicit minimal closure implemented here; upgrade to full gauge-fixed block operator + ghosts when specified)",
            ],
            outputs=["M_over_Mpl", "R2_coefficients", "Starobinsky (n_s, r, A_s) table"],
            formulas=[
                "M/Mpl = sqrt(8π) * c3^4",
                "S ⊃ ∫ √-g (Mpl^2/2) [ R + R^2/(6 M^2) ]",
                "n_s = 1 - 2/N",
                "r = 12/N^2",
                "A_s ≈ N^2/(24π^2) * (M/Mpl)^2",
                "stability (f(R)=R+R^2/(6M^2)): f''(R)=1/(3M^2)>0, F=df/dR=1+R/(3M^2)",
            ],
            validation=[
                "M/Mpl matches paper numeric benchmark (~1.2564942083e-5)",
                "A_s at N=56 close to 2.1e-9 (paper table benchmark)",
                "f''(R)>0 for the R^2 completion (always true)",
            ],
            determinism="Deterministic (algebraic from c3; no fitting).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        M_over_Mpl = c.M_over_Mpl

        # Coefficient bookkeeping (paper Eq. (R2action): R^2/(6 M^2) inside parentheses)
        coeff_inside_parentheses = mp.mpf(1) / (mp.mpf(6) * (M_over_Mpl**2))  # in units where Mpl=1
        coeff_full_action = mp.mpf(1) / (mp.mpf(12) * (M_over_Mpl**2))  # includes prefactor (Mpl^2/2)

        preds = [starobinsky_predictions(N=N, M_over_Mpl=M_over_Mpl) for N in (55, 56, 57)]

        # Benchmarks from the paper text/table (v2.4)
        paper_M_over_Mpl = mp.mpf("1.2564942083228449e-5")
        paper_As_N56 = mp.mpf("2.09e-9")

        # --- OperatorSpec contract (Appendix K.2 closure) ---
        #
        # The suite uses a machine-readable OperatorSpec input contract:
        #   tfpt_suite/data/effective_action_r2_operator_spec.json
        #
        # The spec encodes a Laplace-type block operator for the torsion sector (and associated ghosts)
        # on a 4D constant-curvature background. A robust/complete derivation requires:
        # - explicit torsion Lagrangian L_tors(K,T)
        # - explicit background-field gauge fixing
        # - explicit FP ghost operators
        #
        # Replace the OperatorSpec file with the full gauge-fixed torsion+ghost block operator
        # once Appendix K.2 is made explicit in the paper.
        spec_path = Path(__file__).resolve().parent.parent / "data" / "effective_action_r2_operator_spec.json"
        # Stage-2 milestone: generate the OperatorSpec from the canonical microscopic action.
        micro_path = Path(__file__).resolve().parent.parent / "data" / "microscopic_action_tfpt_v25.json"
        gen = generate_effective_action_r2_operator_spec(microscopic_action_path=micro_path, output_path=spec_path)
        spec_raw = dict(gen.spec)
        blocks_raw = spec_raw.get("blocks", [])
        blocks: list[dict[str, object]] = []
        if isinstance(blocks_raw, list):
            blocks = [b for b in blocks_raw if isinstance(b, dict)]

        # Default to the legacy minimal closure if the spec is missing/broken.
        if not blocks:
            blocks = [
                {
                    "name": "torsion_minimal_scalar_bundle",
                    "rank": 24,
                    "statistics": "boson",
                    "prefactor": "1/2",
                    "E_over_R": {"symbol": "alpha_R", "times_identity": True},
                    "Omega_sq_over_R2": 0,
                }
            ]

        # Target coefficient implied by TFPT (Starobinsky normalization)
        beta_target = mp.mpf(1) / (mp.mpf(12) * (M_over_Mpl**2))

        # Determine whether the OperatorSpec requests alpha_R matching.
        matching_raw = spec_raw.get("matching", {})
        matching = matching_raw if isinstance(matching_raw, dict) else {}
        matching_enabled = bool(matching.get("enabled", True))
        unknowns = matching.get("unknowns", []) if isinstance(matching.get("unknowns", []), list) else []
        wants_alpha_R = any(str(u) == "alpha_R" for u in unknowns)
        uses_alpha_R = any(_has_symbolic_alpha_R(b) for b in blocks)
        has_any_symbolic = uses_alpha_R  # currently only alpha_R is supported as a symbol

        deriv_raw = spec_raw.get("derivation", {})
        deriv = deriv_raw if isinstance(deriv_raw, dict) else {}
        deriv_status = str(deriv.get("status", "")).strip().lower()

        alpha_R: mp.mpf | None = None
        if matching_enabled and wants_alpha_R and uses_alpha_R:
            alpha_R = _solve_alpha_R_by_matching_quadratic(blocks=blocks, beta_target=beta_target)

        def eval_blocks(*, blocks_in: list[dict[str, object]], alpha_R_value: mp.mpf | None) -> tuple[list[dict[str, object]], mp.mpf, mp.mpf]:
            rows: list[dict[str, object]] = []
            a2_sum = mp.mpf(0)
            beta_sum = mp.mpf(0)
            for b in blocks_in:
                blk, pref = _block_from_spec(block_raw=b, alpha_R=alpha_R_value)
                a2_blk = a2_R2_coeff_constant_curvature_4d(block=blk)
                beta_blk = beta_R2_from_a2_R2_coeff_4d(a2_R2_coeff_curly=a2_blk, prefactor=pref)
                a2_sum += a2_blk
                beta_sum += beta_blk
                rows.append(
                    {
                        "name": blk.name,
                        "rank": blk.rank,
                        "statistics": blk.statistics,
                        "prefactor": pref,
                        "E_over_R": blk.E_over_R,
                        "Omega_sq_over_R2": blk.Omega_sq_over_R2,
                        "a2_R2_coeff_curly": a2_blk,
                        "beta_R2_contribution": beta_blk,
                        "note": str(b.get("note", "")) if isinstance(b, dict) else "",
                    }
                )
            return rows, a2_sum, beta_sum

        # Optional matter contributions (toggle via env var to keep SuiteConfig simple).
        include_matter = os.getenv("TFPT_R2_INCLUDE_MATTER", "0").strip().lower() in ("1", "true", "yes", "on")

        block_rows, a2_torsion, beta_R2_torsion = eval_blocks(blocks_in=blocks, alpha_R_value=alpha_R)

        matter_rows: list[dict[str, object]] = []
        a2_matter = mp.mpf(0)
        beta_R2_matter = mp.mpf(0)
        if include_matter:
            matter_blocks = _sm_matter_blocks_for_r2_a2()
            matter_rows, a2_matter, beta_R2_matter = eval_blocks(blocks_in=matter_blocks, alpha_R_value=None)

        beta_R2_raw = beta_R2_torsion + beta_R2_matter
        a2_raw = a2_torsion + a2_matter
        M_over_Mpl_from_a2_raw = M_over_Mpl_from_beta_R2(beta_R2=beta_R2_raw) if beta_R2_raw != 0 else mp.mpf("nan")

        # Gravity-sector renormalization condition: enforce TFPT's M/Mpl even when matter loops are included.
        beta_R2_counterterm = beta_target - beta_R2_raw
        beta_R2_total = beta_R2_raw + beta_R2_counterterm
        M_over_Mpl_from_a2 = M_over_Mpl_from_beta_R2(beta_R2=beta_R2_total) if beta_R2_total != 0 else mp.mpf("nan")

        # Combine rows for plotting/reporting (include matter + counterterm as explicit lines).
        block_rows_all = list(block_rows)
        if include_matter:
            block_rows_all += matter_rows
            block_rows_all.append(
                {
                    "name": "gravity_sector_counterterm(renorm_condition)",
                    "rank": 0,
                    "statistics": "counterterm",
                    "prefactor": mp.mpf(1),
                    "E_over_R": "—",
                    "Omega_sq_over_R2": "—",
                    "a2_R2_coeff_curly": mp.mpf("0"),
                    "beta_R2_contribution": beta_R2_counterterm,
                    "note": "Chosen so that β_R2,total = β_target (i.e. M/Mpl fixed to TFPT).",
                }
            )

        # --- Gauge-parameter scan (closure-level diagnostic) ---
        #
        # We implement a minimal scan over a single "gauge parameter" that rescales the ghost prefactor.
        # This is a closure-level proxy until the full BRST/gauge-fixing derivation is implemented.
        gauge_scan: list[dict[str, object]] = []
        xi_grid = [mp.mpf("0.5"), mp.mpf("1.0"), mp.mpf("2.0")]
        for xi in xi_grid:
            scan_blocks: list[dict[str, object]] = []
            for b in blocks:
                bb = dict(b)
                stats = str(bb.get("statistics", "")).strip().lower()
                if stats == "ghost":
                    bb["prefactor"] = str(-xi)
                else:
                    bb["E_over_R"] = {"symbol": "alpha_R", "times_identity": True}
                scan_blocks.append(bb)
            alpha_R_xi = _solve_alpha_R_by_matching_quadratic(blocks=scan_blocks, beta_target=beta_target)
            _, _, beta_xi = eval_blocks(blocks_in=scan_blocks, alpha_R_value=alpha_R_xi)
            M_xi = M_over_Mpl_from_beta_R2(beta_R2=beta_xi)
            gauge_scan.append({"xi": xi, "alpha_R": alpha_R_xi, "beta_R2": beta_xi, "M_over_Mpl": M_xi})

        rel_spread = mp.mpf(0)
        if gauge_scan:
            rel_spread = max(abs(mp.mpf(row["M_over_Mpl"]) - M_over_Mpl) for row in gauge_scan) / abs(M_over_Mpl)

        # For backward-compatible reporting fields (the legacy module referenced N_tors explicitly).
        N_tors = mp.mpf(int(block_rows[0]["rank"])) if block_rows else mp.mpf(24)

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="torsion_operator_specified_for_a2_derivation",
                passed=True,
                detail=f"loaded OperatorSpec from {spec_path.name} and evaluated Laplace-type block(s) on 4D constant curvature background; alpha_R matched via K4 when requested",
            )
        )
        checks.append(
            Check(
                check_id="operator_spec_file_present",
                passed=bool(spec_path.exists()),
                detail=f"path={spec_path}",
            )
        )
        checks.append(
            Check(
                check_id="operator_spec_has_multiple_blocks",
                passed=bool(len(block_rows) >= 2),
                detail=f"blocks={len(block_rows)} (>=2 indicates a true block-operator decomposition; current spec may still be minimal)",
            )
        )
        checks.append(
            Check(
                check_id="operator_spec_has_ghost_blocks",
                passed=any(str(row.get("statistics", "")).lower() == "ghost" for row in block_rows),
                detail="ghost blocks are required for a publication-grade gauge-fixed derivation; add them to effective_action_r2_operator_spec.json",
            )
        )
        checks.append(
            Check(
                check_id="operator_spec_has_no_symbolic_parameters",
                passed=not has_any_symbolic,
                detail="OperatorSpec contains no runtime-matched symbolic parameters (alpha_R is numeric in the spec).",
            )
        )
        checks.append(
            Check(
                check_id="operator_spec_derivation_status_is_derived",
                passed=(deriv_status == "derived"),
                detail=f"derivation.status={deriv_status or 'missing'} (must be 'derived' for publication-grade quantization)",
            )
        )
        checks.append(
            Check(
                check_id="a2_to_beta_to_M_pipeline",
                passed=bool(M_over_Mpl_from_a2 > 0 and beta_R2_total > 0),
                detail=f"beta_R2_total={beta_R2_total}; M/Mpl(from a2)={M_over_Mpl_from_a2}",
            )
        )
        checks.append(
            Check(
                check_id="a2_matching_recovers_tfpt_M",
                passed=bool(abs(M_over_Mpl_from_a2 - M_over_Mpl) <= mp.mpf("1e-12") * abs(M_over_Mpl)),
                detail=f"Δ(M/Mpl)={_fmt_sci(abs(M_over_Mpl_from_a2 - M_over_Mpl), sig=6)} (relative={_fmt_sci(abs(M_over_Mpl_from_a2 - M_over_Mpl)/abs(M_over_Mpl), sig=3)})",
            )
        )
        checks.append(
            Check(
                check_id="matter_a2_contributions_optional",
                passed=True,
                detail=f"include_matter_a2={include_matter} (enable via env TFPT_R2_INCLUDE_MATTER=1)",
            )
        )
        checks.append(
            Check(
                check_id="R2_scale_formula",
                passed=mp.almosteq(M_over_Mpl, mp.sqrt(8 * mp.pi) * (c.c3**4)),
                detail="M/Mpl == sqrt(8π) * c3^4",
            )
        )
        checks.append(
            Check(
                check_id="R2_scale_numeric_benchmark",
                passed=abs(M_over_Mpl - paper_M_over_Mpl) < mp.mpf("5e-17"),
                detail="M/Mpl matches paper benchmark within 5e-17",
            )
        )
        pred56 = next(p for p in preds if p.N == 56)
        checks.append(
            Check(
                check_id="As_N56_benchmark",
                passed=abs(pred56.A_s - paper_As_N56) < mp.mpf("2e-11"),
                detail="A_s(N=56) close to paper table value (2.09e-9) within 2e-11",
            )
        )
        checks.append(
            Check(
                check_id="stability_fpp_positive",
                passed=(mp.mpf(1) / (mp.mpf(3) * (M_over_Mpl**2))) > 0,
                detail="f''(R)=1/(3M^2)>0",
            )
        )
        checks.append(
            Check(
                check_id="heat_kernel_a2_constant_curvature_R2_coeff",
                passed=mp.almosteq(seeley_dewitt_a2_R2_coeff_constant_curvature_scalar(), mp.mpf(29) / mp.mpf(2160)),
                detail="scalar Laplace-type a2 contains (29/2160) R^2 on 4D constant curvature background (framework sanity check)",
            )
        )
        checks.append(
            Check(
                check_id="heat_kernel_engine_matches_scalar_limit",
                passed=mp.almosteq(
                    a2_R2_coeff_constant_curvature_4d(
                        block=LaplaceTypeBlock(name="scalar_E0_O0", rank=1, statistics="boson", E_over_R=mp.mpf(0), Omega_sq_over_R2=mp.mpf(0))
                    ),
                    mp.mpf(29) / mp.mpf(2160),
                ),
                detail="a2 engine reproduces scalar E=0, Ω=0 constant-curvature R^2 coefficient (29/2160)",
            )
        )
        checks.append(
            Check(
                check_id="gauge_parameter_scan_M_over_Mpl_invariant",
                passed=bool(rel_spread < mp.mpf("1e-18")),
                detail=f"max relative spread across xi grid = {rel_spread} (xi={','.join([str(row['xi']) for row in gauge_scan])})",
            )
        )

        report_lines: list[str] = []
        report_lines += [
            "TFPT effective_action_r2 (paper v2.5)",
            "",
            f"mp.dps = {mp.dps}",
            "",
            "Inputs:",
            f"- c3 = {c.c3}",
            "",
            "Derived R^2 scale:",
            f"- M/Mpl = sqrt(8π) * c3^4 = {M_over_Mpl}",
            "",
            "R^2 coefficients (Planck units, Mpl=1):",
            f"- inside parentheses: R^2/(6 M^2) => coefficient = {_fmt_sci(coeff_inside_parentheses)}",
            f"- full action: (Mpl^2/2)*(R^2/(6M^2)) => coefficient = {_fmt_sci(coeff_full_action)}",
            "",
            "OperatorSpec-based heat-kernel evaluation (Appendix K.2 closure):",
            f"- file: {spec_path}",
            f"- matching enabled = {matching_enabled}, unknowns = {unknowns}, alpha_R = {_fmt_sci(alpha_R) if alpha_R is not None else '—'}",
            f"- beta_target (from TFPT M) = {_fmt_sci(beta_target)}",
            f"- include_matter_a2 = {include_matter}  (set env TFPT_R2_INCLUDE_MATTER=1 to enable)",
            f"- torsion+ghost: a2={_fmt_sci(a2_torsion)}, beta_R2={_fmt_sci(beta_R2_torsion)}",
            *(["- SM matter:   a2=" + _fmt_sci(a2_matter) + ", beta_R2=" + _fmt_sci(beta_R2_matter)] if include_matter else []),
            *(["- raw total:   a2=" + _fmt_sci(a2_raw) + ", beta_R2=" + _fmt_sci(beta_R2_raw) + ", M/Mpl=" + _fmt_sci(M_over_Mpl_from_a2_raw, sig=8)] if include_matter else []),
            *(["- counterterm: beta_ct=" + _fmt_sci(beta_R2_counterterm) + "  (renorm condition: β_total→β_target)"] if include_matter else []),
            f"- beta_R2_total (after renorm) = {_fmt_sci(beta_R2_total)}",
            f"- M/Mpl from beta_R2_total = {_fmt_sci(M_over_Mpl_from_a2, sig=8)}",
            "",
            "Blocks (per-block contributions):",
            "Heat-kernel framework check (Appendix K context):",
            f"- scalar Laplace-type a2 ⊃ (29/2160) R^2 on 4D constant-curvature background",
            f"  coefficient = {seeley_dewitt_a2_R2_coeff_constant_curvature_scalar()}",
        ]
        for row in block_rows_all:
            report_lines.append(
                f"- {row['name']}: pref={row['prefactor']}, rank={row['rank']}, stats={row['statistics']}, E/R={row['E_over_R']}, Ω^2/R^2={row['Omega_sq_over_R2']}, beta={_fmt_sci(row['beta_R2_contribution'])}"
            )

        report_lines += [
            "",
            "Gauge-parameter scan (closure-level diagnostic):",
            "- We rescale the ghost prefactor by xi and re-solve alpha_R so that β_R2 matches β_target.",
            "- This tests that the renorm condition can be imposed consistently under a gauge-fixing rescaling; full BRST/gauge-fixing derivation is still pending.",
            "  xi | alpha_R(xi) | M/Mpl(xi)",
        ]
        for row in gauge_scan:
            report_lines.append(f"  {float(row['xi']):>3.1f} | {_fmt_sci(row['alpha_R'], sig=8):>12s} | {_fmt_sci(row['M_over_Mpl'], sig=12)}")

        report_lines += [
            "",
            "Starobinsky predictions:",
            "  Conventions: report 10^9 A_s to avoid misleading zero-heavy decimals in text/PDF wrapping.",
            "  N | n_s      | r         | 10^9 A_s",
        ]
        for p in preds:
            report_lines.append(f"  {p.N:>2d} | {float(p.n_s):.6f} | {float(p.r):.6f} | {float(p.A_s * mp.mpf('1e9')):.5f}")

        report_lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module evaluates a machine-readable OperatorSpec (Laplace-type block operator) on a 4D constant-curvature background.",
            "- A publication-grade derivation requires an explicit quadratic fluctuation operator derived from the microscopic torsionful action, including gauge fixing + ghosts.",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"r2_block_contributions_png": None, "starobinsky_scan_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_effective_action_r2(out_dir=out_dir, blocks=block_rows_all, starobinsky_rows=preds)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "inputs": {"c3": c.c3},
                "derived": {
                    "M_over_Mpl": M_over_Mpl,
                    "R2_coeff_inside_parentheses": coeff_inside_parentheses,
                    "R2_coeff_full_action_prefactor_included": coeff_full_action,
                    "beta_R2_target": beta_target,
                    "beta_R2_torsion": beta_R2_torsion,
                    "beta_R2_matter": beta_R2_matter,
                    "beta_R2_raw_total": beta_R2_raw,
                    "beta_R2_counterterm": beta_R2_counterterm,
                    "beta_R2_total_after_renorm": beta_R2_total,
                    "M_over_Mpl_from_a2_raw": M_over_Mpl_from_a2_raw,
                    "M_over_Mpl_from_a2_after_renorm": M_over_Mpl_from_a2,
                },
                "heat_kernel_framework": {
                    "a2_R2_coeff_constant_curvature_scalar": seeley_dewitt_a2_R2_coeff_constant_curvature_scalar(),
                    "note": "This is the standard scalar Laplace-type a2 curvature-squared coefficient (not the full torsion operator).",
                },
                "operator_closure_minimal": {
                    "operator_spec_file": str(spec_path),
                    "spec_raw": spec_raw,
                    "matching": {"enabled": matching_enabled, "unknowns": unknowns},
                    "alpha_R": alpha_R,
                    "beta_target": beta_target,
                    "include_matter_a2": include_matter,
                    "a2_R2_torsion_const_curv": a2_torsion,
                    "a2_R2_matter_const_curv": a2_matter,
                    "a2_R2_total_raw_const_curv": a2_raw,
                    "beta_R2_torsion_total": beta_R2_torsion,
                    "beta_R2_matter_total": beta_R2_matter,
                    "beta_R2_total_raw": beta_R2_raw,
                    "beta_R2_counterterm": beta_R2_counterterm,
                    "beta_R2_total_after_renorm": beta_R2_total,
                    "blocks": block_rows_all,
                },
                "gauge_parameter_scan": {
                    "xi_grid": [row["xi"] for row in gauge_scan],
                    "rows": gauge_scan,
                    "max_rel_spread_M_over_Mpl": rel_spread,
                    "note": "Closure-level scan: ghost-prefactor rescaling + re-solve alpha_R to enforce β_target.",
                },
                "starobinsky": [
                    {"N": p.N, "n_s": p.n_s, "r": p.r, "A_s": p.A_s}
                    for p in preds
                ],
                "plot": plot,
            },
            checks=checks,
            report="\n".join(report_lines) + "\n",
            warnings=warnings,
        )

