from __future__ import annotations

import gzip
import json
import lzma
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

import numpy as np


DatasetKind = Literal["multivariate_gaussian", "one_sided_upper"]


@dataclass(frozen=True)
class DatasetResult:
    dataset_id: str
    kind: str
    labels: list[str]
    chi2: float
    dof: int
    loglike: float
    details: dict[str, Any]


def _as_float(x: object) -> float:
    try:
        return float(x)  # type: ignore[arg-type]
    except Exception:
        return float("nan")


def _safe_logdet(mat: np.ndarray) -> float:
    sign, logdet = np.linalg.slogdet(mat)
    if sign <= 0:
        return float("nan")
    return float(logdet)


def _chi2_sf(chi2: float, dof: int) -> float:
    """
    Chi-square survival function. We keep this dependency-light by using SciPy only if present.
    If SciPy isn't present, return NaN (module can still report chi2 deterministically).
    """
    try:
        from scipy.stats import chi2 as chi2_dist  # type: ignore

        return float(chi2_dist.sf(float(chi2), int(dof)))
    except Exception:
        return float("nan")


def load_likelihood_spec(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _open_text_path(path: Path):
    if path.suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def load_grid_table(
    path: Path,
    *,
    param_columns: Sequence[str],
    chi2_column: str = "chi2",
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(f"grid file not found: {path}")
    header: list[str] | None = None
    rows: list[list[float]] = []
    with _open_text_path(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if header is None and any(not _is_float_token(p) for p in parts):
                header = [p.strip() for p in parts]
                continue
            try:
                rows.append([float(p) for p in parts])
            except Exception:
                continue
    if not rows:
        raise ValueError(f"grid file {path} has no numeric rows")

    data = np.array(rows, dtype=float)
    if header:
        col_index = {name: idx for idx, name in enumerate(header)}
        if chi2_column not in col_index:
            raise ValueError(f"grid file {path} missing chi2 column '{chi2_column}'")
        for col in param_columns:
            if col not in col_index:
                raise ValueError(f"grid file {path} missing parameter column '{col}'")
        param_idx = [col_index[c] for c in param_columns]
        chi2_idx = col_index[chi2_column]
    else:
        expected_cols = len(param_columns) + 1
        if data.shape[1] != expected_cols:
            raise ValueError(f"grid file {path} expected {expected_cols} columns, got {data.shape[1]}")
        param_idx = list(range(len(param_columns)))
        chi2_idx = data.shape[1] - 1

    points = data[:, param_idx]
    chi2 = data[:, chi2_idx]
    meta = {
        "rows": int(data.shape[0]),
        "param_columns": list(param_columns),
        "chi2_column": chi2_column,
        "has_header": bool(header),
        "path": str(path),
    }
    return points, chi2, meta


def interpolate_grid_chi2(
    *,
    points: np.ndarray,
    chi2: np.ndarray,
    target: Sequence[float],
    method: str = "inverse_distance",
    k_nearest: int = 8,
    scales: Sequence[float] | None = None,
) -> tuple[float, dict[str, Any]]:
    if points.ndim != 2:
        raise ValueError("grid points must be 2D array")
    if len(target) != points.shape[1]:
        raise ValueError(f"target dimension {len(target)} != grid dimension {points.shape[1]}")
    if points.shape[0] == 0:
        raise ValueError("grid points empty")

    tgt = np.array(target, dtype=float)
    scl = np.ones(points.shape[1], dtype=float) if scales is None else np.array(scales, dtype=float)
    if np.any(scl <= 0) or scl.shape[0] != points.shape[1]:
        raise ValueError("invalid scales for grid interpolation")

    diffs = (points - tgt) / scl
    d2 = np.sum(diffs * diffs, axis=1)
    if method == "nearest" or k_nearest <= 1 or points.shape[0] == 1:
        idx = int(np.argmin(d2))
        return float(chi2[idx]), {"method": "nearest", "k_used": 1, "d2_min": float(d2[idx])}

    k_use = min(int(k_nearest), int(points.shape[0]))
    idx = np.argpartition(d2, k_use - 1)[:k_use]
    d2_sel = d2[idx]
    if np.any(d2_sel == 0.0):
        idx0 = int(idx[int(np.argmin(d2_sel))])
        return float(chi2[idx0]), {"method": "nearest_exact", "k_used": 1, "d2_min": float(d2_sel.min())}
    weights = 1.0 / np.maximum(d2_sel, 1.0e-12)
    chi2_val = float(np.sum(weights * chi2[idx]) / np.sum(weights))
    meta = {"method": "inverse_distance", "k_used": int(k_use), "d2_min": float(d2_sel.min())}
    return chi2_val, meta


def _is_float_token(token: str) -> bool:
    try:
        float(token)
        return True
    except Exception:
        return False


def evaluate_dataset(
    *,
    dataset: Mapping[str, Any],
    predictions: Mapping[str, float],
    include_norm_constant: bool = True,
) -> DatasetResult:
    kind = str(dataset.get("kind", "")).strip()
    dataset_id = str(dataset.get("dataset_id", "")).strip() or "unnamed_dataset"

    if kind == "multivariate_gaussian":
        labels = [str(x) for x in (dataset.get("labels", []) or [])]
        means = dataset.get("means", [])
        cov = dataset.get("covariance", [])
        floors = dataset.get("theory_floors", {}) if isinstance(dataset.get("theory_floors", {}), dict) else {}

        y = np.array([_as_float(predictions.get(k, float("nan"))) for k in labels], dtype=float)
        mu = np.array([_as_float(m) for m in means], dtype=float)
        C = np.array(cov, dtype=float)

        if C.shape != (len(labels), len(labels)):
            raise ValueError(f"{dataset_id}: covariance has wrong shape {C.shape}, expected {(len(labels), len(labels))}")

        # Apply theory floors as diagonal additions (σ_floor^2).
        if floors:
            diag = np.zeros(len(labels), dtype=float)
            for i, k in enumerate(labels):
                sig = _as_float(floors.get(k, 0.0))
                if math.isfinite(sig) and sig > 0:
                    diag[i] = sig * sig
            C = C + np.diag(diag)

        r = y - mu
        C_inv = np.linalg.inv(C)
        chi2 = float(r.T @ C_inv @ r)
        dof = int(len(labels))
        logdet = _safe_logdet(C)
        if include_norm_constant and math.isfinite(logdet):
            loglike = float(-0.5 * (chi2 + logdet + dof * math.log(2.0 * math.pi)))
        else:
            loglike = float(-0.5 * chi2)

        p = _chi2_sf(chi2=chi2, dof=dof)
        return DatasetResult(
            dataset_id=dataset_id,
            kind=kind,
            labels=labels,
            chi2=chi2,
            dof=dof,
            loglike=loglike,
            details={
                "predictions": {k: float(predictions.get(k, float("nan"))) for k in labels},
                "means": [float(x) for x in mu.tolist()],
                "covariance": C.tolist(),
                "p_value": p,
                "include_norm_constant": bool(include_norm_constant),
            },
        )

    if kind == "one_sided_upper":
        label = str(dataset.get("label", "")).strip()
        upper = _as_float(dataset.get("upper", float("nan")))
        pred = _as_float(predictions.get(label, float("nan")))
        if not math.isfinite(pred) or not math.isfinite(upper) or upper <= 0:
            raise ValueError(f"{dataset_id}: invalid pred/upper for one-sided bound")
        chi2 = 0.0 if pred <= upper else float(((pred - upper) / upper) ** 2)
        return DatasetResult(
            dataset_id=dataset_id,
            kind=kind,
            labels=[label],
            chi2=float(chi2),
            dof=1,
            loglike=float(-0.5 * chi2),
            details={"pred": float(pred), "upper": float(upper), "cl": _as_float(dataset.get("cl", float("nan"))), "likelihood": str(dataset.get("likelihood", ""))},
        )

    raise ValueError(f"{dataset_id}: unsupported dataset kind: {kind}")


def evaluate_likelihood_spec(
    *,
    spec: Mapping[str, Any],
    predictions: Mapping[str, float],
    include_norm_constant: bool = True,
) -> tuple[list[DatasetResult], float]:
    datasets = spec.get("datasets", [])
    if not isinstance(datasets, list):
        raise ValueError("likelihood spec: datasets must be a list")

    results: list[DatasetResult] = []
    loglike_total = 0.0
    for d in datasets:
        if not isinstance(d, dict):
            continue
        enabled = d.get("enabled", True)
        if enabled is False:
            continue
        res = evaluate_dataset(dataset=d, predictions=predictions, include_norm_constant=include_norm_constant)
        results.append(res)
        loglike_total += float(res.loglike)
    return results, float(loglike_total)


def predictions_from_global_consistency_results(payload: Mapping[str, Any]) -> dict[str, float]:
    """
    Extracts the canonical prediction keys used by TFPT scorecards from a
    `global_consistency_test/results.json` payload.
    """
    preds: dict[str, float] = {}
    try:
        res = payload.get("results", {}) if isinstance(payload.get("results", {}), dict) else {}
        p = res.get("predictions", {}) if isinstance(res.get("predictions", {}), dict) else {}
        terms = res.get("terms", []) if isinstance(res.get("terms", []), list) else []

        def _term_pred(name: str) -> float | None:
            if not isinstance(terms, list):
                return None
            for row in terms:
                if isinstance(row, dict) and row.get("name") == name:
                    value = row.get("pred", row.get("value", float("nan")))
                    v = _as_float(value)
                    return v if math.isfinite(v) else None
            return None

        for k in ["alpha_inv_0", "alpha_bar5_inv_MZ", "beta_deg", "cabibbo_lambda", "n_s", "r", "A_s"]:
            value = None
            if isinstance(p, dict) and k in p:
                v = _as_float(p[k])
                value = v if math.isfinite(v) else None
            if value is None:
                value = _term_pred(k)
            if value is not None:
                preds[k] = float(value)
    except Exception:
        pass
    return preds

