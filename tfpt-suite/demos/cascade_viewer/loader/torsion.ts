/**
 * Torsion Bounds Loader
 * =====================
 * Loads torsion bounds mapping and TFPT predictions.
 * 
 * Data sources:
 *   - out/torsion_bounds_mapping/results.json
 *   - tfpt_suite/data/torsion_regimes.json
 * 
 * Paper reference: Section 6 (Torsion), Appendix M (SME mapping)
 * 
 * Key concepts:
 *   - S_μ: TFPT axial torsion vector
 *   - b_μ = k · S_μ (SME mapping, k = 3/4 for minimal coupling)
 *   - ξ = c₃/φ₀ (normalization factor)
 */

import * as fs from 'fs';

/** A single torsion bound entry */
export interface TorsionBound {
  label: string;
  coefficient: string;
  abs_max_input_GeV: number;
  inferred_abs_max_S_mu_GeV: number;
  mapping_used: string;
  applies_to: string;
}

/** Torsion regime definition */
export interface TorsionRegime {
  id: string;
  label: string;
  model: string;
  c_factor?: number;
  H0_km_s_Mpc?: number;
}

/** TFPT torsion prediction */
export interface TorsionPrediction {
  S_mu_abs_GeV: number;
  regime: string;
  regime_details: TorsionRegime;
  theorem_basis: string;
}

/** Mapping constants */
export interface TorsionMapping {
  k_default: number;
  k_policy: string;
  xi_tree: number;
  xi: number;
  convention: string;
}

/** Full torsion data */
export interface TorsionData {
  /** Inferred S_μ bounds from SME data */
  inferred_bounds: TorsionBound[];
  /** TFPT prediction for current regime */
  prediction: TorsionPrediction;
  /** Mapping constants */
  mapping: TorsionMapping;
  /** Ratio of prediction to bounds (all should be << 1) */
  prediction_to_bound_ratios: { label: string; ratio: number }[];
  /** Available regimes */
  regimes: TorsionRegime[];
  /** Data file path */
  data_file: string;
  /** Checks */
  checks: { check_id: string; passed: boolean; detail: string }[];
}

/**
 * Load torsion bounds and predictions.
 * This function does NOT compute anything - it only reads and structures.
 */
export function loadTorsion(boundsMappingPath: string, regimesPath: string): TorsionData {
  const raw = JSON.parse(fs.readFileSync(boundsMappingPath, 'utf8'));
  const results = raw.results;
  
  // Load regimes
  const regimesRaw = JSON.parse(fs.readFileSync(regimesPath, 'utf8'));
  const regimes: TorsionRegime[] = regimesRaw.regimes.map((r: any) => ({
    id: r.id,
    label: r.label,
    model: r.model,
    c_factor: r.c_factor,
    H0_km_s_Mpc: r.H0_km_s_Mpc
  }));
  
  // Parse inferred bounds
  const inferred_bounds: TorsionBound[] = results.inferred_bounds.map((b: any) => ({
    label: b.label,
    coefficient: b.coefficient,
    abs_max_input_GeV: b.abs_max_input_GeV,
    inferred_abs_max_S_mu_GeV: b.inferred_abs_max_S_mu_GeV,
    mapping_used: b.mapping_used,
    applies_to: b.applies_to
  }));
  
  // Parse prediction
  const pred = results.tfpt_prediction;
  const prediction: TorsionPrediction = {
    S_mu_abs_GeV: pred.S_mu_abs_GeV,
    regime: pred.regime,
    regime_details: pred.regime_details,
    theorem_basis: pred.theorem_basis
  };
  
  // Parse mapping
  const mapRaw = results.mapping;
  const mapping: TorsionMapping = {
    k_default: mapRaw.k_default,
    k_policy: mapRaw.k_policy,
    xi_tree: parseFloat(mapRaw.xi_tree),
    xi: parseFloat(mapRaw.xi),
    convention: mapRaw.mapping_raw?.convention ?? ''
  };
  
  // Compute prediction/bound ratios
  const prediction_to_bound_ratios = inferred_bounds.map(b => ({
    label: b.label,
    ratio: prediction.S_mu_abs_GeV / b.inferred_abs_max_S_mu_GeV
  }));
  
  // Extract checks
  const checks = raw.checks.map((c: any) => ({
    check_id: c.check_id,
    passed: c.passed,
    detail: c.detail
  }));
  
  return {
    inferred_bounds,
    prediction,
    mapping,
    prediction_to_bound_ratios,
    regimes,
    data_file: results.data_file,
    checks
  };
}
