/**
 * Alpha Precision / Homeostasis Loader
 * =====================================
 * Loads α(0) self-consistency data and k-sensitivity from alpha_precision_audit.
 * 
 * Data source:
 *   - out/alpha_precision_audit/results.json
 * 
 * Paper reference: Section 3.2 (CFE cubic equation), Eq. (3.7)
 * 
 * The CFE (Cubic Fixed-point Equation):
 *   α³ - 2c₃³α² - 8b₁c₃⁶ ln(1/φ) = 0
 * 
 * With backreaction:
 *   φ(α) = φ_tree + δ_top exp(-kα)
 */

import * as fs from 'fs';

/** A single k-sensitivity data point */
export interface KSensitivityPoint {
  k: number;
  alpha_inv: string;
  alpha_inv_numeric: number;
  ppm_vs_codata: number;
  converged: boolean;
  iterations: number;
}

/** Self-consistent solution data */
export interface SelfConsistentSolution {
  k: number;
  alpha_inv: string;
  alpha_inv_numeric: number;
  ppm_vs_codata: number;
  converged: boolean;
  iterations: number;
}

/** Defect expansion diagnostics */
export interface DefectDiagnostics {
  delta2_match_codata_at_k2: number;
  delta2_over_delta_top: number;
  k_match_codata: number;
  gamma_defect_required: number;
}

/** Reference values */
export interface AlphaReference {
  alpha_inv_codata_2022: number;
  alpha_bar5_inv_MZ: number;
  alpha_bar5_inv_MZ_sigma: number;
}

/** Defect partition g=5 result (theoryv3) */
export interface DefectPartitionG5 {
  g_value: number;
  delta2: number;
  delta2_over_delta_top2: number;
  model_id: string;
  alpha_inv_pred: number;
  alpha_inv_ref: number;
  alpha_inv_sigma: number;
  z_score: number;
  negative_control: { g: number; z: number }[];
}

/** Full alpha precision data */
export interface AlphaData {
  /** Baseline (fixed φ₀) */
  baseline: {
    alpha_inv: string;
    alpha_inv_numeric: number;
    ppm_vs_codata: number;
  };
  /** Self-consistent solution (standard k=2) */
  self_consistent: SelfConsistentSolution;
  /** Two-defect partition solution */
  self_consistent_two_defect_half: SelfConsistentSolution | null;
  /** k-sensitivity table */
  k_sensitivity: KSensitivityPoint[];
  /** Diagnostics */
  diagnostics: DefectDiagnostics;
  /** Reference values */
  reference: AlphaReference;
  /** Checks passed */
  checks: { check_id: string; passed: boolean; detail: string }[];
  /** NEW: g=5 defect partition (theoryv3) */
  defect_g5: DefectPartitionG5 | null;
}

/**
 * Load alpha precision audit data.
 * This function does NOT compute anything - it only reads and structures.
 */
export function loadAlpha(resultsPath: string, defectG5Path?: string): AlphaData {
  const raw = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
  const results = raw.results;
  
  // Parse baseline
  const baseline = {
    alpha_inv: results.baseline.alpha_inv,
    alpha_inv_numeric: parseFloat(results.baseline.alpha_inv),
    ppm_vs_codata: parseFloat(results.baseline.ppm_vs_codata)
  };
  
  // Parse self-consistent solution
  const sc = results.self_consistent;
  const self_consistent: SelfConsistentSolution = {
    k: parseFloat(sc.k),
    alpha_inv: sc.alpha_inv,
    alpha_inv_numeric: parseFloat(sc.alpha_inv),
    ppm_vs_codata: parseFloat(sc.ppm_vs_codata),
    converged: sc.converged,
    iterations: sc.iterations
  };
  
  // Parse two-defect solution if available
  let self_consistent_two_defect_half: SelfConsistentSolution | null = null;
  if (results.self_consistent_two_defect_half) {
    const td = results.self_consistent_two_defect_half;
    self_consistent_two_defect_half = {
      k: parseFloat(td.k),
      alpha_inv: td.alpha_inv,
      alpha_inv_numeric: parseFloat(td.alpha_inv),
      ppm_vs_codata: parseFloat(td.ppm_vs_codata),
      converged: td.converged,
      iterations: td.iterations
    };
  }
  
  // Parse k-sensitivity table
  const k_sensitivity: KSensitivityPoint[] = results.k_sensitivity.map((pt: any) => ({
    k: parseFloat(pt.k),
    alpha_inv: pt.alpha_inv,
    alpha_inv_numeric: parseFloat(pt.alpha_inv),
    ppm_vs_codata: parseFloat(pt.ppm_vs_codata),
    converged: pt.converged,
    iterations: pt.iterations
  }));
  
  // Parse diagnostics
  const diag = results.diagnostics;
  const diagnostics: DefectDiagnostics = {
    delta2_match_codata_at_k2: parseFloat(diag.delta2_match_codata_at_k2),
    delta2_over_delta_top: parseFloat(diag.delta2_over_delta_top),
    k_match_codata: parseFloat(diag.k_match_codata),
    gamma_defect_required: parseFloat(diag.gamma_defect_required)
  };
  
  // Parse reference
  const ref = results.reference;
  const reference: AlphaReference = {
    alpha_inv_codata_2022: parseFloat(ref.alpha_inv_codata_2022),
    alpha_bar5_inv_MZ: parseFloat(ref.alpha_bar5_inv_MZ),
    alpha_bar5_inv_MZ_sigma: parseFloat(ref.alpha_bar5_inv_MZ_sigma)
  };
  
  // Extract checks
  const checks = raw.checks.map((c: any) => ({
    check_id: c.check_id,
    passed: c.passed,
    detail: c.detail
  }));
  
  // Load g=5 defect partition if path provided
  let defect_g5: DefectPartitionG5 | null = null;
  if (defectG5Path) {
    try {
      const g5Raw = JSON.parse(fs.readFileSync(defectG5Path, 'utf8'));
      const g5Results = g5Raw.results;
      defect_g5 = {
        g_value: g5Results.delta2.g_value,
        delta2: parseFloat(g5Results.delta2.delta2),
        delta2_over_delta_top2: parseFloat(g5Results.delta2.delta2_over_delta_top2),
        model_id: g5Results.delta2.model_id,
        alpha_inv_pred: parseFloat(g5Results.alpha_inv_0.pred),
        alpha_inv_ref: parseFloat(g5Results.alpha_inv_0.ref),
        alpha_inv_sigma: parseFloat(g5Results.alpha_inv_0.sigma),
        z_score: parseFloat(g5Results.alpha_inv_0.z),
        negative_control: g5Results.negative_control.z_by_g.map((item: any) => ({
          g: item.g,
          z: item.z
        }))
      };
    } catch (e) {
      console.warn('Could not load defect_g5 data:', e);
    }
  }
  
  return {
    baseline,
    self_consistent,
    self_consistent_two_defect_half,
    k_sensitivity,
    diagnostics,
    reference,
    checks,
    defect_g5
  };
}
