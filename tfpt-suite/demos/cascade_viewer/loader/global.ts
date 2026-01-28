/**
 * Global Consistency Loader
 * =========================
 * Loads global χ² scorecard from global_consistency_test.
 * 
 * Data sources:
 *   - out/global_consistency_test/results.json
 *   - tfpt_suite/data/global_reference.json
 * 
 * Paper reference: Table 1 (Observable predictions), Section 7 (Global fit)
 * 
 * Key metrics:
 *   - χ² per observable
 *   - z-score (number of σ from reference)
 *   - ppm deviation
 */

import * as fs from 'fs';

/** A single χ² term / observable */
export interface Chi2Term {
  name: string;
  pred: number;
  mean: number;
  sigma_exp: number;
  sigma_theory: number;
  sigma_total: number;
  z_total: number;
  chi2_total: number;
  ppm: number;
  source: string;
}

/** Global totals */
export interface GlobalTotals {
  chi2_total: number;
  dof: number;
  p_value: number;
  chi2_excluding_alpha: number;
  dof_excluding_alpha: number;
  p_value_excluding_alpha: number;
}

/** TFPT predictions summary */
export interface PredictionsSummary {
  alpha_inv_0: number;
  n_s: number;
  A_s: number;
  r: number;
  cabibbo_lambda: number;
  beta_deg: number;
}

/** Full global consistency data */
export interface GlobalData {
  /** All χ² terms */
  terms: Chi2Term[];
  /** Watchlist (terms with |z| > threshold) */
  watchlist: Chi2Term[];
  /** Global totals */
  totals: GlobalTotals;
  /** TFPT predictions */
  predictions: PredictionsSummary;
  /** Mode (engineering/physics) */
  mode: string;
  /** Checks */
  checks: { check_id: string; passed: boolean; detail: string }[];
}

/**
 * Load global consistency test data.
 * This function does NOT compute anything - it only reads and structures.
 */
export function loadGlobal(resultsPath: string): GlobalData {
  const raw = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
  const results = raw.results;
  
  // Parse χ² terms
  const terms: Chi2Term[] = results.terms.map((t: any) => ({
    name: t.name,
    pred: parseFloat(t.pred),
    mean: parseFloat(t.mean),
    sigma_exp: parseFloat(t.sigma_exp),
    sigma_theory: parseFloat(t.sigma_theory),
    sigma_total: parseFloat(t.sigma_total),
    z_total: parseFloat(t.z_total),
    chi2_total: parseFloat(t.chi2_total),
    ppm: parseFloat(t.ppm),
    source: t.source
  }));
  
  // Parse watchlist
  const watchlist: Chi2Term[] = (results.watchlist_engineering ?? []).map((t: any) => ({
    name: t.name,
    pred: parseFloat(t.pred),
    mean: parseFloat(t.mean),
    sigma_exp: parseFloat(t.sigma_exp),
    sigma_theory: parseFloat(t.sigma_theory),
    sigma_total: parseFloat(t.sigma_total),
    z_total: parseFloat(t.z_total),
    chi2_total: parseFloat(t.chi2_total),
    ppm: parseFloat(t.ppm),
    source: t.source
  }));
  
  // Parse totals
  const totalsRaw = results.totals;
  const totals: GlobalTotals = {
    chi2_total: parseFloat(totalsRaw.chi2_total_engineering),
    dof: totalsRaw.dof_core,
    p_value: parseFloat(totalsRaw.p_core_engineering),
    chi2_excluding_alpha: parseFloat(totalsRaw.chi2_total_engineering_excluding_alpha),
    dof_excluding_alpha: totalsRaw.dof_core_excluding_alpha,
    p_value_excluding_alpha: parseFloat(totalsRaw.p_core_engineering_excluding_alpha)
  };
  
  // Parse predictions
  const pred = results.predictions;
  const predictions: PredictionsSummary = {
    alpha_inv_0: parseFloat(pred.alpha_inv_0 ?? pred.alpha_inv_0_two_defect),
    n_s: parseFloat(pred.n_s),
    A_s: parseFloat(pred.A_s),
    r: parseFloat(pred.r),
    cabibbo_lambda: parseFloat(pred.cabibbo_lambda),
    beta_deg: parseFloat(pred.beta_deg)
  };
  
  // Extract checks
  const checks = raw.checks.map((c: any) => ({
    check_id: c.check_id,
    passed: c.passed,
    detail: c.detail
  }));
  
  return {
    terms,
    watchlist,
    totals,
    predictions,
    mode: results.mode,
    checks
  };
}
