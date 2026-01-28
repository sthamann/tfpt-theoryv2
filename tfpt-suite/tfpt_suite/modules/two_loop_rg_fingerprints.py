from __future__ import annotations

import csv
import hashlib
import importlib
import json
import math
from dataclasses import dataclass
from math import log10
from pathlib import Path
import sys
from typing import Optional

import numpy as np

from tfpt_suite.conventions import g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


@dataclass(frozen=True)
class GaugeRunningTable:
    source_path: str
    mu_GeV: np.ndarray
    log10_mu: np.ndarray
    alpha3: np.ndarray
    alpha2: np.ndarray
    alpha1_gut: np.ndarray
    alphaY: np.ndarray


def _workspace_root() -> Path:
    # .../wolfram_latex_attachments/tfpt-suite/tfpt_suite/modules/<file>.py -> parents[3] is workspace root
    return Path(__file__).resolve().parents[3]


def _default_pyrate3_pythonoutput_dir() -> Path:
    """
    Default location of the PyR@TE3-generated PythonOutput package we use for the numeric two-loop run.
    """
    sm_dir = (
        _workspace_root()
        / "Pyrate3"
        / "pyrate"
        / "results"
        / "SM_TFPT_2loop_v25"
        / "SM_TFPT_2Loop_v25"
        / "PythonOutput"
    )
    if sm_dir.is_dir():
        return sm_dir
    return (
        _workspace_root()
        / "Pyrate3"
        / "pyrate"
        / "results"
        / "E8Cascade2LoopGravityV2"
        / "PythonOutput"
    )


def _default_model_module_for_dir(pythonoutput_dir: Path) -> str:
    if pythonoutput_dir.name == "PythonOutput" and pythonoutput_dir.parent.name == "SM_TFPT_2Loop_v25":
        return "SM_TFPT_2Loop_v25"
    return "E8Cascade2LoopGravityV2"


def _load_cfg() -> dict[str, object]:
    cfg_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "two_loop_rg_fingerprints.json"
    if not cfg_path.exists():
        return {"_cfg_path": str(cfg_path), "_cfg_loaded": False}
    raw = json.loads(cfg_path.read_text(encoding="utf-8"))
    raw["_cfg_path"] = str(cfg_path)
    raw["_cfg_loaded"] = True
    return raw


def _require_expected_model_in_used(*, expected: str, used_dir: Path, used_module: str) -> None:
    exp = str(expected).strip()
    if not exp:
        raise ValueError("model_name_expected must be a non-empty string")
    hay = (str(used_dir) + " " + str(used_module)).lower()
    if exp.lower() not in hay:
        raise RuntimeError(
            "Model mismatch (fail-fast): "
            f"model_name_expected={exp!r} not found in pythonoutput_dir={str(used_dir)!r} "
            f"or model_module={str(used_module)!r}"
        )


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _load_sm_inputs_mz() -> SmMzInputs:
    sm_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "sm_inputs_mz.json"
    raw = json.loads(sm_path.read_text(encoding="utf-8"))
    return SmMzInputs(
        mu_GeV=float(raw["mu_GeV"]),
        alpha_em_inv=float(raw["alpha_em_inv"]),
        sin2_thetaW=float(raw["sin2_thetaW"]),
        alpha_s=float(raw["alpha_s"]),
    )


def _generate_gauge_table_from_pyrate3_pythonoutput(
    *,
    pythonoutput_dir: Path,
    model_module: str,
    output_csv: Path,
    npoints: int,
    gravity_alpha3_enabled: bool,
    gravity_kappa_vector: tuple[float, float, float],
    gravity_c3: float,
    gravity_alpha_definition_for_u1: str,
) -> dict[str, object]:
    """
    Numerically solve the PyR@TE3 RGEs (PythonOutput) and write a small canonical gauge table.

    Output schema matches the loader in this module:
      mu_GeV, log10_mu, alphaY, alpha1_GUT, alpha2, alpha3
    """
    if not pythonoutput_dir.is_dir():
        raise FileNotFoundError(f"PyR@TE3 PythonOutput directory not found: {pythonoutput_dir}")
    if npoints < 2:
        raise ValueError("npoints must be >= 2")

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    # Import generated model module from the given PythonOutput directory.
    sys.path.insert(0, str(pythonoutput_dir))
    try:
        mod = importlib.import_module(model_module)
    finally:
        try:
            sys.path.remove(str(pythonoutput_dir))
        except ValueError:
            pass

    sm_inp = _load_sm_inputs_mz()
    mu0 = float(sm_inp.mu_GeV)
    if mu0 <= 0:
        raise ValueError("Invalid SM input scale mu_GeV")

    # PyR@TE convention: t = log10(mu/GeV). Start exactly at μ=MZ (paper wording: PDG boundary at MZ).
    tmin = float(log10(mu0))
    tmax = 19.0
    t0 = float(tmin)
    rge = mod.RGEsolver("rge", tmin=tmin, tmax=tmax, initialScale=t0)

    rge.loops = {
        "GaugeCouplings": 2,
        "Yukawas": 2,
        "QuarticTerms": 2,
        "TrilinearTerms": 2,
        "ScalarMasses": 2,
        "Vevs": 2,
    }

    # Optional TFPT runner patch: Δ(dα_i)/d ln μ = κ_i * c3 * α_i^3, implemented as a correction to dg_i/dt (t=log10 μ).
    # We use αY for U(1) consistently with PyR@TE's convention g1 ≡ gY.
    if gravity_alpha3_enabled:
        alpha_def = str(gravity_alpha_definition_for_u1).strip().lower()
        if alpha_def not in ("alphay", "alpha_y", "alpha-y"):
            raise ValueError(f"Unsupported alpha_definition_for_U1={gravity_alpha_definition_for_u1!r}; expected 'alphaY'")
        k1, k2, k3 = [float(x) for x in gravity_kappa_vector]
        c3 = float(gravity_c3)
        ln10 = float(np.log(10.0))
        pi = float(np.pi)
        orig_beta = rge.betaFunction

        def patched_beta(t: float, couplingsArray):  # type: ignore[no-untyped-def]
            dy = orig_beta(t, couplingsArray)
            # Couplings are in the solver's native convention: g1 ≡ gY, g2, g3.
            g1_sm, g2, g3 = rge.extractCouplings(couplingsArray, "GaugeCouplings")
            for i, (g, kk) in enumerate(((g1_sm, k1), (g2, k2), (g3, k3))):
                if kk == 0.0 or c3 == 0.0:
                    continue
                gi = float(g)
                # dα/d ln μ = κ c3 α^3  ⇒  dg/dt = ln(10) * κ c3 * g^5 / (32 π^2)
                dy[i] = float(dy[i]) + ln10 * kk * c3 * (gi**5) / (32.0 * (pi**2))
            return dy

        rge.betaFunction = patched_beta  # type: ignore[assignment]

    # Initial conditions at μ=MZ from the shared SM input table (same source as CKM module).
    #
    # Convention:
    # - In this PyR@TE PythonOutput, we treat rge.g1 as the SM hypercharge coupling g' (= g_Y),
    #   and we *report* both α_Y = g'^2/(4π) and α_1(GUT) = (5/3) α_Y.
    g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(sm_inp)  # g1 is GUT-normalized here
    gut_norm = g1_gut_over_gY()
    gY = gY_from_g1_gut(g1_gut_mz)

    rge.g1.initialValue = float(gY)
    rge.g2.initialValue = float(g2_mz)
    rge.g3.initialValue = float(g3_mz)

    rge.Yu.initialValue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.94]]
    rge.Yd.initialValue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.017]]
    rge.Ye.initialValue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.010]]
    if hasattr(rge, "yN"):
        rge.yN.initialValue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    if hasattr(rge, "lambda_"):
        rge.lambda_.initialValue = 0.13
    if hasattr(rge, "lPhi"):
        rge.lPhi.initialValue = 0.0
    if hasattr(rge, "lHphi"):
        rge.lHphi.initialValue = 0.0
    if hasattr(rge, "cR3"):
        rge.cR3.initialValue = 0.0
    if hasattr(rge, "mu2"):
        rge.mu2.initialValue = 0.0
    if hasattr(rge, "MPhi"):
        rge.MPhi.initialValue = 0.0
    if hasattr(rge, "vSM"):
        rge.vSM.initialValue = 246.0
    if hasattr(rge, "vPQ"):
        rge.vPQ.initialValue = 0.0

    rge.solve(Npoints=int(npoints))

    t = np.array(rge.tList, dtype=float)
    mu = np.array(10 ** t, dtype=float)

    g1_sm = np.array(rge.solutions["g1"], dtype=float)
    g2 = np.array(rge.solutions["g2"], dtype=float)
    g3 = np.array(rge.solutions["g3"], dtype=float)

    gut_norm = g1_gut_over_gY()
    g1_gut_arr = gut_norm * g1_sm

    alphaY = (g1_sm**2) / (4.0 * np.pi)
    alpha1_gut = (g1_gut_arr**2) / (4.0 * np.pi)
    alpha2 = (g2**2) / (4.0 * np.pi)
    alpha3 = (g3**2) / (4.0 * np.pi)

    # Optional additional couplings (future RG pipeline upgrades: 2-loop Yukawas + λ_H).
    yt = None
    yb = None
    ytau = None
    lamH = None
    try:
        if "Yu" in rge.solutions:
            Yu = np.array(rge.solutions["Yu"], dtype=float)
            yt = Yu[:, 2, 2]
        if "Yd" in rge.solutions:
            Yd = np.array(rge.solutions["Yd"], dtype=float)
            yb = Yd[:, 2, 2]
        if "Ye" in rge.solutions:
            Ye = np.array(rge.solutions["Ye"], dtype=float)
            ytau = Ye[:, 2, 2]
        if "lambda_" in rge.solutions:
            lamH = np.array(rge.solutions["lambda_"], dtype=float)
    except Exception:
        yt = None
        yb = None
        ytau = None
        lamH = None

    with output_csv.open("w", encoding="utf-8", newline="") as f:
        extra_fields: list[str] = []
        if yt is not None:
            extra_fields.append("yt")
        if yb is not None:
            extra_fields.append("yb")
        if ytau is not None:
            extra_fields.append("ytau")
        if lamH is not None:
            extra_fields.append("lambdaH")

        writer = csv.DictWriter(
            f,
            fieldnames=["mu_GeV", "log10_mu", "alphaY", "alpha1_GUT", "alpha2", "alpha3"] + extra_fields,
        )
        writer.writeheader()
        for i in range(int(t.size)):
            row = {
                "mu_GeV": float(mu[i]),
                "log10_mu": float(t[i]),
                "alphaY": float(alphaY[i]),
                "alpha1_GUT": float(alpha1_gut[i]),
                "alpha2": float(alpha2[i]),
                "alpha3": float(alpha3[i]),
            }
            if yt is not None:
                row["yt"] = float(yt[i])
            if yb is not None:
                row["yb"] = float(yb[i])
            if ytau is not None:
                row["ytau"] = float(ytau[i])
            if lamH is not None:
                row["lambdaH"] = float(lamH[i])
            writer.writerow(row)

    return {
        "pythonoutput_dir": _relpath(pythonoutput_dir),
        "model_module": model_module,
        "gravity_alpha3_patch": {
            "enabled": bool(gravity_alpha3_enabled),
            "kappa_vector": [float(x) for x in gravity_kappa_vector],
            "c3": float(gravity_c3),
            "alpha_definition_for_U1": str(gravity_alpha_definition_for_u1),
        },
        "tmin_log10": tmin,
        "tmax_log10": tmax,
        "initialScale_log10": t0,
        "npoints": int(npoints),
        "initial_conditions": {
            "mu0_GeV": float(mu0),
            "gY": float(rge.g1.initialValue),
            "g1_gut": float(g1_gut_mz),
            "g2": float(rge.g2.initialValue),
            "g3": float(rge.g3.initialValue),
            "alpha_s_MZ": float(sm_inp.alpha_s),
            "Yu33": 0.94,
            "Yd33": 0.017,
            "Ye33": 0.010,
            "lambda": 0.13,
        },
        "output_csv": _relpath(output_csv),
    }


def _load_pyrate_gauge_table(path: Path) -> GaugeRunningTable:
    mu: list[float] = []
    log10_mu: list[float] = []
    alpha3: list[float] = []
    alpha2: list[float] = []
    alpha1_gut: list[float] = []
    alphaY: list[float] = []

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mu.append(float(row["mu_GeV"]))
            log10_mu.append(float(row["log10_mu"]))
            alpha3.append(float(row["alpha3"]))
            alpha2.append(float(row["alpha2"]))
            alpha1_gut.append(float(row["alpha1_GUT"]))
            if "alphaY" in row and row["alphaY"] not in (None, ""):
                alphaY.append(float(row["alphaY"]))

    mu_arr = np.array(mu, dtype=float)
    log10_arr = np.array(log10_mu, dtype=float)
    alpha3_arr = np.array(alpha3, dtype=float)
    alpha2_arr = np.array(alpha2, dtype=float)
    alpha1_arr = np.array(alpha1_gut, dtype=float)
    if len(alphaY) == len(alpha1_gut):
        alphaY_arr = np.array(alphaY, dtype=float)
    else:
        # Backwards compatibility: older CSV schema stored only alpha1_GUT, and alphaY = (3/5) alpha1_GUT.
        alphaY_arr = (3.0 / 5.0) * alpha1_arr

    if not (mu_arr.size >= 2 and np.all(np.diff(log10_arr) > 0)):
        raise ValueError(f"Unexpected ordering / size in RG table: {path}")

    return GaugeRunningTable(
        source_path=_relpath(path),
        mu_GeV=mu_arr,
        log10_mu=log10_arr,
        alpha3=alpha3_arr,
        alpha2=alpha2_arr,
        alpha1_gut=alpha1_arr,
        alphaY=alphaY_arr,
    )


def _interp_in_log10(table: GaugeRunningTable, *, mu_GeV: float, y: np.ndarray) -> float:
    if mu_GeV <= 0:
        raise ValueError("mu_GeV must be positive")
    x = float(log10(mu_GeV))
    x0 = float(table.log10_mu[0])
    x1 = float(table.log10_mu[-1])
    if x < x0 or x > x1:
        raise ValueError(f"mu_GeV={mu_GeV} out of table range: 10^{x0}..10^{x1} GeV")
    return float(np.interp(x, table.log10_mu, y))


@dataclass(frozen=True)
class Crossing:
    target: float
    mu_GeV: float
    alpha_at_mu: float
    rel_dev: float
    bracket: Optional[tuple[int, int]]


def _find_crossing(table: GaugeRunningTable, *, target: float, y: np.ndarray) -> Crossing:
    # Prefer an actual sign-change bracket (monotone crossings), otherwise fall back to best gridpoint.
    diff = y - float(target)
    bracket: Optional[tuple[int, int]] = None
    for i in range(int(diff.size - 1)):
        a = float(diff[i])
        b = float(diff[i + 1])
        if a == 0.0:
            mu = float(table.mu_GeV[i])
            return Crossing(target=target, mu_GeV=mu, alpha_at_mu=float(y[i]), rel_dev=0.0, bracket=(i, i))
        if (a > 0 and b < 0) or (a < 0 and b > 0):
            bracket = (i, i + 1)
            break

    if bracket is None:
        j = int(np.argmin(np.abs(diff)))
        mu = float(table.mu_GeV[j])
        alpha = float(y[j])
        rel = float(abs(alpha - target) / abs(target)) if target != 0 else float("nan")
        return Crossing(target=target, mu_GeV=mu, alpha_at_mu=alpha, rel_dev=rel, bracket=None)

    i0, i1 = bracket
    x0 = float(table.log10_mu[i0])
    x1 = float(table.log10_mu[i1])
    y0 = float(y[i0])
    y1 = float(y[i1])
    if y1 == y0:
        x_star = 0.5 * (x0 + x1)
    else:
        x_star = x0 + (target - y0) * (x1 - x0) / (y1 - y0)
    mu_star = float(10 ** x_star)
    alpha_star = float(np.interp(x_star, table.log10_mu, y))
    rel = float(abs(alpha_star - target) / abs(target)) if target != 0 else float("nan")
    return Crossing(target=target, mu_GeV=mu_star, alpha_at_mu=alpha_star, rel_dev=rel, bracket=bracket)


def _plot_two_loop_overview(*, table: GaugeRunningTable, out_dir: Path, const: TfptConstants) -> list[str]:
    """
    Write one or more PNG plots into out_dir.
    Returns a list of warning strings (empty if successful).
    """
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        mu = table.mu_GeV
        a1 = table.alpha1_gut
        a2 = table.alpha2
        a3 = table.alpha3

        # Plot 1: alpha_i running
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_xscale("log")
        ax.plot(mu, a1, label=r"$\alpha_1=(5/3)\alpha_Y$ (GUT norm.)", linewidth=2.0)
        ax.plot(mu, a2, label=r"$\alpha_2$", linewidth=2.0)
        ax.plot(mu, a3, label=r"$\alpha_3$", linewidth=2.0)
        ax.axvline(1e6, linestyle="--", linewidth=1.5, alpha=0.7, color="black")
        ax.axhline(float(const.varphi0), linestyle=":", linewidth=1.5, alpha=0.8, color="black", label=r"$\varphi_0$")
        ax.axhline(float(const.c3), linestyle=":", linewidth=1.5, alpha=0.8, color="gray", label=r"$c_3$")
        ax.set_xlabel(r"$\mu$ [GeV]")
        ax.set_ylabel(r"$\alpha_i(\mu)$")
        ax.set_title("Two-loop gauge running (PyR@TE3 numeric solve)")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(out_dir / "two_loop_gauge_running.png", dpi=200)
        plt.close(fig)

        # Plot 2: alpha3 fingerprint zoom
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        ax2.set_xscale("log")
        ax2.plot(mu, a3, label=r"$\alpha_3(\mu)$", linewidth=2.5)
        ax2.axvline(1e6, linestyle="--", linewidth=1.5, alpha=0.7, color="black", label="1 PeV")
        ax2.axvline(2.5e8, linestyle="--", linewidth=1.5, alpha=0.7, color="gray", label=r"$\mu_{c_3}$ (benchmark)")
        ax2.axhline(float(const.varphi0), linestyle=":", linewidth=2.0, alpha=0.8, color="black", label=r"$\varphi_0$")
        ax2.axhline(float(const.c3), linestyle=":", linewidth=2.0, alpha=0.8, color="gray", label=r"$c_3$")
        ax2.set_xlim(max(mu.min(), 1e4), min(mu.max(), 1e10))
        ax2.set_xlabel(r"$\mu$ [GeV]")
        ax2.set_ylabel(r"$\alpha_3(\mu)$")
        ax2.set_title(r"RG fingerprints (zoom): $\alpha_3$ vs $(\varphi_0,c_3)$")
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc="best")
        fig2.tight_layout()
        fig2.savefig(out_dir / "two_loop_rg_fingerprints_zoom.png", dpi=200)
        plt.close(fig2)

        # Plot 3: alpha^{-1} unification overview with mu* annotation
        inv1 = 1.0 / a1
        inv2 = 1.0 / a2
        inv3 = 1.0 / a3
        log_mu = np.log10(mu)
        diff12 = np.abs(inv1 - inv2)
        idx = int(np.nanargmin(diff12)) if diff12.size else 0
        mu_star = float(mu[idx]) if mu.size else float("nan")
        log_mu_star = float(log_mu[idx]) if log_mu.size else float("nan")
        inv12_avg = float((inv1[idx] + inv2[idx]) / 2.0) if inv1.size else float("nan")
        mismatch = float(abs(inv3[idx] - inv12_avg)) if inv3.size else float("nan")

        fig3, ax3 = plt.subplots(figsize=(10, 6))
        ax3.plot(log_mu, inv1, label=r"$\alpha_1^{-1}$", linewidth=2.0)
        ax3.plot(log_mu, inv2, label=r"$\alpha_2^{-1}$", linewidth=2.0)
        ax3.plot(log_mu, inv3, label=r"$\alpha_3^{-1}$", linewidth=2.0)
        if math.isfinite(log_mu_star):
            ax3.axvline(log_mu_star, linestyle="--", linewidth=1.5, alpha=0.7, color="black", label=f"$\\mu_*$≈{mu_star:.2e} GeV")
            ax3.annotate(
                f"mismatch≈{mismatch:.2g}",
                xy=(log_mu_star, inv12_avg),
                xytext=(log_mu_star + 0.2, inv12_avg),
                arrowprops={"arrowstyle": "->", "alpha": 0.6},
            )
        ax3.set_xlabel(r"$\log_{10}(\mu/\mathrm{GeV})$")
        ax3.set_ylabel(r"$\alpha_i^{-1}(\mu)$")
        ax3.set_title("Gauge-coupling unification overview")
        ax3.grid(True, alpha=0.3)
        ax3.legend(loc="best")
        fig3.tight_layout()
        fig3.savefig(out_dir / "gauge_unification_running.png", dpi=200)
        plt.close(fig3)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")
    return warnings


class TwoLoopRgFingerprintsModule(TfptModule):
    module_id = "two_loop_rg_fingerprints"
    title = "Two-loop RG fingerprints (α3 scale matches to TFPT invariants)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Two-loop gauge running table (generated from PyR@TE3 PythonOutput and cached into this module output directory as `gauge_couplings.csv`)",
                "PyR@TE3 PythonOutput package selected via `tfpt_suite/data/two_loop_rg_fingerprints.json` (fail-fast, no silent fallback)",
                "SM boundary conditions at μ=MZ: `tfpt-suite/tfpt_suite/data/sm_inputs_mz.json`",
                "TFPT invariants from paper v2.5: c3=1/(8π), varphi0=1/(6π)+3/(256π^4)",
            ],
            outputs=[
                "α3(1 PeV) and its relative deviation to varphi0",
                "scales μ where α3 crosses varphi0 and c3 (interpolated in log10 μ)",
            ],
            formulas=[
                "c3 = 1/(8π)",
                "varphi0 = 1/(6π) + 3/(256π^4)",
                "relative deviation: |α3(μ)-target|/target",
                "interpolation/crossing: linear in log10(μ) vs α3 using the tabulated points",
            ],
            validation=[
                "table ordering: log10(μ) strictly increasing",
                "α3(1 PeV) is within 2% of varphi0 (paper claim-map falsification target)",
                "targets (varphi0, c3) lie within the table’s α3 range (so a crossing exists)",
            ],
            determinism="Deterministic given PyR@TE3 PythonOutput + fixed numeric integration grid + TFPT constants.",
        )

    def run(self, config) -> ModuleResult:
        const = TfptConstants.compute()
        c3 = float(const.c3)
        varphi0 = float(const.varphi0)

        out_dir = self.output_dir(config)
        generated_csv = out_dir / "gauge_couplings.csv"
        generated_csv_grav = out_dir / "gauge_couplings_gravity.csv"

        sm_inp = _load_sm_inputs_mz()
        mu0 = float(sm_inp.mu_GeV)
        tmin_target = float(log10(mu0))

        cfg = _load_cfg()
        cfg_path = Path(str(cfg.get("_cfg_path", ""))) if str(cfg.get("_cfg_path", "")).strip() else None
        cfg_loaded = bool(cfg.get("_cfg_loaded", False))

        # Preferred config path (schema_version>=2): select via central PyR@TE registry.
        pythonoutput_kind = str(cfg.get("pyrate_pythonoutput_kind", "")).strip()
        if pythonoutput_kind:
            mdl = get_pyrate_pythonoutput(pythonoutput_kind)
            pythonoutput_dir = Path(mdl.pythonoutput_dir)
            model_name_expected = str(mdl.model_name_expected).strip()
            # PyR@TE PythonOutput convention: module file matches the model name.
            model_module = model_name_expected
            yaml_source_path = Path(mdl.yaml_source) if mdl.yaml_source is not None else None
        else:
            # Legacy schema (v1) fallback: explicit directory + module fields, or deterministic default.
            model_name_expected = str(cfg.get("model_name_expected", "")).strip()
            if not model_name_expected:
                # Backwards-compatible deterministic fallback: follow the legacy directory preference.
                pythonoutput_dir = _default_pyrate3_pythonoutput_dir()
                model_module = _default_model_module_for_dir(pythonoutput_dir)
                model_name_expected = model_module
            else:
                pythonoutput_dir = _workspace_root() / str(cfg.get("pythonoutput_dir", "")).strip()
                model_module = str(cfg.get("model_module", "")).strip() or model_name_expected
            yaml_source_path = None
            if "yaml_source" in cfg and str(cfg.get("yaml_source", "")).strip():
                yaml_source_path = _workspace_root() / str(cfg.get("yaml_source", "")).strip()

        _require_expected_model_in_used(expected=model_name_expected, used_dir=pythonoutput_dir, used_module=model_module)

        pythonoutput_module_file = pythonoutput_dir / f"{model_module}.py"
        if not pythonoutput_module_file.exists():
            raise FileNotFoundError(f"Expected PyR@TE3 module file not found: {pythonoutput_module_file}")

        # Optional gravity α^3 patch (explicit; runner-level, not in PyR@TE model YAML).
        grav_cfg = cfg.get("gravity_alpha3_patch", {}) if isinstance(cfg.get("gravity_alpha3_patch", {}), dict) else {}
        gravity_enabled = bool(grav_cfg.get("enabled", False))
        alpha_definition_for_u1 = str(grav_cfg.get("alpha_definition_for_U1", "alphaY")).strip()
        kappa_vec_raw = grav_cfg.get("kappa_vector", [0.0, 0.0, 0.0])
        if not (isinstance(kappa_vec_raw, list) and len(kappa_vec_raw) == 3):
            raise ValueError("gravity_alpha3_patch.kappa_vector must be a length-3 list")
        gravity_kappa_vector = (float(kappa_vec_raw[0]), float(kappa_vec_raw[1]), float(kappa_vec_raw[2]))

        generation_info: dict[str, object] | None = None
        generation_info_grav: dict[str, object] | None = None
        table: GaugeRunningTable | None = None
        table_grav: GaugeRunningTable | None = None

        def regenerate() -> None:
            nonlocal generation_info, table
            generation_info = _generate_gauge_table_from_pyrate3_pythonoutput(
                pythonoutput_dir=pythonoutput_dir,
                model_module=model_module,
                output_csv=generated_csv,
                npoints=171,  # fixed for determinism + speed
                gravity_alpha3_enabled=False,
                gravity_kappa_vector=(0.0, 0.0, 0.0),
                gravity_c3=c3,
                gravity_alpha_definition_for_u1=alpha_definition_for_u1,
            )
            table = None

        def regenerate_gravity() -> None:
            nonlocal generation_info_grav, table_grav
            generation_info_grav = _generate_gauge_table_from_pyrate3_pythonoutput(
                pythonoutput_dir=pythonoutput_dir,
                model_module=model_module,
                output_csv=generated_csv_grav,
                npoints=171,
                gravity_alpha3_enabled=True,
                gravity_kappa_vector=gravity_kappa_vector,
                gravity_c3=c3,
                gravity_alpha_definition_for_u1=alpha_definition_for_u1,
            )
            table_grav = None

        if config.overwrite or (not generated_csv.exists()):
            regenerate()
        else:
            # Guard against stale cached tables from older conventions (e.g. starting at 10^2 GeV).
            try:
                table = _load_pyrate_gauge_table(generated_csv)
                # Guard against stale cached tables from a *different model*.
                prev_results_path = out_dir / "results.json"
                if prev_results_path.exists():
                    prev = json.loads(prev_results_path.read_text(encoding="utf-8"))
                    prev_fp = (prev.get("results", {}) or {}).get("model_fingerprint", {}) if isinstance(prev, dict) else {}
                    prev_used = str(prev_fp.get("pythonoutput_module_sha256", "")).strip()
                    cur_used = _sha256_file(pythonoutput_module_file)
                    if prev_used and prev_used != cur_used:
                        regenerate()
                if abs(float(table.log10_mu[0]) - tmin_target) > 1e-3:
                    regenerate()
            except Exception:
                regenerate()

        if table is None:
            table = _load_pyrate_gauge_table(generated_csv)

        if gravity_enabled:
            if config.overwrite or (not generated_csv_grav.exists()):
                regenerate_gravity()
            else:
                try:
                    table_grav = _load_pyrate_gauge_table(generated_csv_grav)
                    if abs(float(table_grav.log10_mu[0]) - tmin_target) > 1e-3:
                        regenerate_gravity()
                except Exception:
                    regenerate_gravity()
            if table_grav is None:
                table_grav = _load_pyrate_gauge_table(generated_csv_grav)

        # basic integrity checks
        alpha3_decreasing = bool(np.all(np.diff(table.alpha3) < 0))
        alpha3_min = float(np.min(table.alpha3))
        alpha3_max = float(np.max(table.alpha3))

        # fingerprint evaluation at a fixed reference scale (1 PeV = 10^6 GeV)
        mu_ref = 1e6
        alpha3_1pev = _interp_in_log10(table, mu_GeV=mu_ref, y=table.alpha3)
        rel_1pev_vs_varphi0 = float(abs(alpha3_1pev - varphi0) / abs(varphi0))

        cross_varphi0 = _find_crossing(table, target=varphi0, y=table.alpha3)
        cross_c3 = _find_crossing(table, target=c3, y=table.alpha3)

        checks: list[Check] = []
        warnings: list[str] = []
        checks.append(
            Check(
                check_id="table_alpha3_monotone_decreasing",
                passed=alpha3_decreasing,
                detail=f"α3 decreases with μ across table (min={alpha3_min:.6g}, max={alpha3_max:.6g})",
            )
        )
        checks.append(
            Check(
                check_id="alpha3_target_in_range_varphi0",
                passed=bool(alpha3_min <= varphi0 <= alpha3_max),
                detail=f"varphi0={varphi0:.6g} within α3 range [{alpha3_min:.6g}, {alpha3_max:.6g}]",
            )
        )
        checks.append(
            Check(
                check_id="alpha3_target_in_range_c3",
                passed=bool(alpha3_min <= c3 <= alpha3_max),
                detail=f"c3={c3:.6g} within α3 range [{alpha3_min:.6g}, {alpha3_max:.6g}]",
            )
        )
        checks.append(
            Check(
                check_id="alpha3_1PeV_close_to_varphi0",
                passed=bool(np.isfinite(rel_1pev_vs_varphi0) and rel_1pev_vs_varphi0 < 0.02),
                detail=f"α3(1 PeV)={alpha3_1pev:.6g}, varphi0={varphi0:.6g}, rel dev={rel_1pev_vs_varphi0:.3%}",
            )
        )

        if config.plot:
            warnings.extend(_plot_two_loop_overview(table=table, out_dir=out_dir, const=const))

        csv_sha256 = _sha256_file(generated_csv) if generated_csv.exists() else None
        csv_grav_sha256 = _sha256_file(generated_csv_grav) if (gravity_enabled and generated_csv_grav.exists()) else None
        pythonoutput_sha256 = _sha256_file(pythonoutput_module_file)
        yaml_source_sha256 = _sha256_file(yaml_source_path) if (yaml_source_path is not None and yaml_source_path.exists()) else None

        gravity_diff: dict[str, object] | None = None
        if gravity_enabled and table_grav is not None:
            # Compare gravity-on to baseline on the baseline μ-grid (interpolated in log10 μ).
            a3_grav_on_base = np.interp(table.log10_mu, table_grav.log10_mu, table_grav.alpha3)
            a2_grav_on_base = np.interp(table.log10_mu, table_grav.log10_mu, table_grav.alpha2)
            a1_grav_on_base = np.interp(table.log10_mu, table_grav.log10_mu, table_grav.alpha1_gut)

            rel_a3 = np.abs((a3_grav_on_base - table.alpha3) / table.alpha3)
            imax = int(np.nanargmax(rel_a3)) if rel_a3.size else 0
            mu_at_max = float(table.mu_GeV[imax]) if table.mu_GeV.size else float("nan")
            max_rel = float(rel_a3[imax]) if rel_a3.size else float("nan")

            # Snapshot diffs at a few reference scales (if within range).
            ref_mus = [1.0e6, 1.0e10, 1.0e16]
            points: list[dict[str, float]] = []
            for mu_ref in ref_mus:
                try:
                    a1b = _interp_in_log10(table, mu_GeV=mu_ref, y=table.alpha1_gut)
                    a2b = _interp_in_log10(table, mu_GeV=mu_ref, y=table.alpha2)
                    a3b = _interp_in_log10(table, mu_GeV=mu_ref, y=table.alpha3)
                    a1g = _interp_in_log10(table_grav, mu_GeV=mu_ref, y=table_grav.alpha1_gut)
                    a2g = _interp_in_log10(table_grav, mu_GeV=mu_ref, y=table_grav.alpha2)
                    a3g = _interp_in_log10(table_grav, mu_GeV=mu_ref, y=table_grav.alpha3)
                    points.append(
                        {
                            "mu_GeV": float(mu_ref),
                            "delta_alpha1_gut": float(a1g - a1b),
                            "delta_alpha2": float(a2g - a2b),
                            "delta_alpha3": float(a3g - a3b),
                            "rel_delta_alpha3": float(abs(a3g - a3b) / abs(a3b)) if a3b != 0 else float("nan"),
                        }
                    )
                except Exception:
                    continue

            gravity_diff = {
                "mu_GeV_at_max_rel_delta_alpha3": mu_at_max,
                "max_rel_delta_alpha3": max_rel,
                "points": points,
            }

        conventions_line = "Conventions: hypercharge Q=T3+Y (SM); α3=g3^2/(4π), αY=g′^2/(4π), α1(GUT)=(5/3)αY."
        g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(sm_inp)
        gut_norm = g1_gut_over_gY()
        gY_mz = gY_from_g1_gut(g1_gut_mz)

        ratio = float(g1_gut_mz / gY_mz) if gY_mz != 0 else float("nan")
        checks.append(
            Check(
                check_id="g1_gut_over_gY_convention",
                passed=bool(np.isfinite(ratio) and abs(ratio - gut_norm) < 1e-12),
                detail=f"g1_GUT/gY={ratio:.15g} vs sqrt(5/3)={gut_norm:.15g}",
            )
        )
        # Global sanity: α1(GUT) must equal (5/3) αY for every row in the cached table.
        alpha_ratio_target = 5.0 / 3.0
        alpha_ratio = table.alpha1_gut / table.alphaY
        max_alpha_ratio_dev = float(np.nanmax(np.abs(alpha_ratio - alpha_ratio_target))) if alpha_ratio.size else float("nan")
        checks.append(
            Check(
                check_id="alpha1_gut_over_alphaY_ratio_global",
                passed=bool(np.isfinite(max_alpha_ratio_dev) and max_alpha_ratio_dev < 1e-12),
                detail=f"max |α1_GUT/αY - 5/3| across table = {max_alpha_ratio_dev:.3e}",
            )
        )

        csv_header = ""
        try:
            if generated_csv.exists():
                csv_header = generated_csv.read_text(encoding="utf-8").splitlines()[0].strip()
        except Exception:
            csv_header = ""

        report = "\n".join(
            [
                "TFPT two-loop RG fingerprints (α3 scale matching)",
                "",
                conventions_line,
                "",
                f"Input boundary: μ0 = {mu0:.6g} GeV (MZ), log10(μ0/GeV) = {tmin_target:.6g}, αs(MZ) = {float(sm_inp.alpha_s):.6g}",
                f"Initial couplings @ μ0: g′={gY_mz:.6g}, g2={g2_mz:.6g}, g3={g3_mz:.6g}  (g1_GUT={g1_gut_mz:.6g})",
                "",
                f"Source CSV: {table.source_path}" + (f"  (sha256={csv_sha256})" if csv_sha256 else ""),
                *(["CSV columns: " + csv_header] if csv_header else []),
                f"Model expected: {model_name_expected}",
                f"PyR@TE3 PythonOutput: {_relpath(pythonoutput_dir)}",
                f"PyR@TE3 module file: {_relpath(pythonoutput_module_file)}" + (f"  (sha256={pythonoutput_sha256})" if pythonoutput_sha256 else ""),
                *( [f"PyR@TE3 YAML source: {_relpath(yaml_source_path)}" + (f"  (sha256={yaml_source_sha256})" if yaml_source_sha256 else "")] if yaml_source_path is not None else [] ),
                *( [f"Module config: {_relpath(cfg_path)}"] if (cfg_loaded and cfg_path is not None) else [] ),
                *(["", f"PyR@TE solve grid: Npoints={generation_info['npoints']}, tmax={generation_info['tmax_log10']:.6g}"] if generation_info is not None else []),
                "",
                "Gravity α^3 patch (runner-level, optional):",
                f"- enabled = {gravity_enabled}",
                f"- kappa_vector = {list(gravity_kappa_vector)}",
                f"- c3 = {c3:.12g} (TFPT constant)",
                f"- alpha_definition_for_U1 = {alpha_definition_for_u1}",
                *(
                    [
                        f"- gravity CSV: {_relpath(generated_csv_grav)}"
                        + (f"  (sha256={csv_grav_sha256})" if csv_grav_sha256 else "")
                    ]
                    if gravity_enabled
                    else []
                ),
                *(
                    [
                        f"- max rel Δα3 across table: {float(gravity_diff['max_rel_delta_alpha3']):.3e} at μ={float(gravity_diff['mu_GeV_at_max_rel_delta_alpha3']):.6g} GeV",
                        *[
                            f"- Δα3(μ={p['mu_GeV']:.3g} GeV) = {p['delta_alpha3']:.3e}  (rel {p['rel_delta_alpha3']:.3e})"
                            for p in (gravity_diff.get('points', []) if isinstance(gravity_diff, dict) else [])
                            if isinstance(p, dict) and 'mu_GeV' in p
                        ],
                    ]
                    if gravity_enabled and gravity_diff is not None
                    else []
                ),
                "",
                "TFPT targets (paper v2.5):",
                f"- c3      = {c3:.12g}",
                f"- varphi0 = {varphi0:.12g}",
                "",
                "Fingerprint evaluations:",
                f"- α3(1 PeV = 1e6 GeV) = {alpha3_1pev:.12g}",
                f"  rel dev vs varphi0  = {rel_1pev_vs_varphi0:.6%}",
                "",
                "Interpolated crossing scales (linear in log10 μ):",
                f"- α3(μ)=varphi0: μ ≈ {cross_varphi0.mu_GeV:.6g} GeV  (α3={cross_varphi0.alpha_at_mu:.6g}, rel dev={cross_varphi0.rel_dev:.3e})",
                f"- α3(μ)=c3:      μ ≈ {cross_c3.mu_GeV:.6g} GeV  (α3={cross_c3.alpha_at_mu:.6g}, rel dev={cross_c3.rel_dev:.3e})",
                "",
                "Checks:",
                *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
                "",
            ]
        )

        return ModuleResult(
            results={
                "model_fingerprint": {
                    "pyrate_pythonoutput_kind": pythonoutput_kind if pythonoutput_kind else None,
                    "model_name_expected": model_name_expected,
                    "pythonoutput_module_file": _relpath(pythonoutput_module_file),
                    "pythonoutput_module_sha256": pythonoutput_sha256,
                    "yaml_source": _relpath(yaml_source_path) if yaml_source_path is not None else None,
                    "yaml_source_sha256": yaml_source_sha256,
                    "config_file": _relpath(cfg_path) if (cfg_loaded and cfg_path is not None) else None,
                },
                "gravity_alpha3_patch": {
                    "enabled": gravity_enabled,
                    "kappa_vector": list(gravity_kappa_vector),
                    "c3": c3,
                    "alpha_definition_for_U1": alpha_definition_for_u1,
                    "source_csv_gravity": _relpath(generated_csv_grav) if gravity_enabled else None,
                    "csv_gravity_sha256": csv_grav_sha256,
                    "diff_summary": gravity_diff,
                },
                "source_csv": table.source_path,
                "generation": generation_info,
                "generation_gravity": generation_info_grav,
                "reproducibility": {
                    "conventions": conventions_line,
                    "sm_boundary": {
                        "mu0_GeV": mu0,
                        "log10_mu0": tmin_target,
                        "alpha_s_MZ": float(sm_inp.alpha_s),
                    },
                    "artifacts": {
                        "csv_sha256": csv_sha256,
                        "csv_gravity_sha256": csv_grav_sha256,
                        "pyrate_pythonoutput_kind": pythonoutput_kind if pythonoutput_kind else None,
                        "pythonoutput_dir": _relpath(pythonoutput_dir),
                        "pythonoutput_module_file": _relpath(pythonoutput_module_file),
                        "pythonoutput_module_sha256": pythonoutput_sha256,
                        "yaml_source": _relpath(yaml_source_path) if yaml_source_path is not None else None,
                        "yaml_source_sha256": yaml_source_sha256,
                        "config_file": _relpath(cfg_path) if (cfg_loaded and cfg_path is not None) else None,
                    },
                },
                "targets": {"c3": c3, "varphi0": varphi0},
                "alpha3": {
                    "range": [alpha3_min, alpha3_max],
                    "is_monotone_decreasing": alpha3_decreasing,
                    "alpha3_at_1PeV": alpha3_1pev,
                    "rel_dev_1PeV_vs_varphi0": rel_1pev_vs_varphi0,
                    "crossing_varphi0": {
                        "mu_GeV": cross_varphi0.mu_GeV,
                        "alpha3": cross_varphi0.alpha_at_mu,
                        "rel_dev": cross_varphi0.rel_dev,
                        "bracket_idx": list(cross_varphi0.bracket) if cross_varphi0.bracket is not None else None,
                    },
                    "crossing_c3": {
                        "mu_GeV": cross_c3.mu_GeV,
                        "alpha3": cross_c3.alpha_at_mu,
                        "rel_dev": cross_c3.rel_dev,
                        "bracket_idx": list(cross_c3.bracket) if cross_c3.bracket is not None else None,
                    },
                },
            },
            checks=checks,
            report=report,
            warnings=warnings,
        )

