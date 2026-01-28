/**
 * Spectrum / Spectral Flow Loader
 * ================================
 * Loads spectral flow data from APS η-gluing and operator eigenvalues.
 * 
 * Data sources:
 *   - out/aps_eta_gluing/results.json (spectral flow table)
 *   - out/effective_action_r2/results.json (operator blocks, E/R eigenvalues)
 *   - out/bounce_perturbations/results.json (transfer functions)
 * 
 * Paper reference: Appendix K (APS gluing), Appendix L (bounce perturbations)
 */

import * as fs from 'fs';

/** A single spectral flow entry */
export interface SpectralFlowEntry {
  m: number;
  spin: 'periodic' | 'antiperiodic';
  analytic_sf: number;
  numeric_sf: number;
  winding_det_u: number;
  eta_start_numeric: number;
  eta_end_numeric: number;
  theta_start: number;
  theta_end: number;
}

/** Operator block from effective action */
export interface OperatorBlock {
  name: string;
  rank: number;
  statistics: 'boson' | 'ghost';
  E_over_R: number;
  Omega_sq_over_R2: number;
  a2_R2_coeff_curly: number;
  beta_R2_contribution: number;
  prefactor: number;
}

/** Bounce perturbation transfer function data */
export interface TransferFunctionData {
  k_grid: number[];
  T_scalar: number[];
  T_tensor: number[];
  k_bounce_s_est: number;
  k_bounce_t_est: number;
}

/** Full spectrum data */
export interface SpectrumData {
  /** APS η-gluing spectral flow table */
  spectral_flow: SpectralFlowEntry[];
  /** Minimal seam term Δ_Γ = 2π for m=1 */
  seam_term_delta_gamma_2pi: number;
  /** Operator blocks from R² effective action */
  operator_blocks: OperatorBlock[];
  /** Total β_R² from operator closure */
  beta_R2_total: number;
  /** M/M_Pl derived from operator closure */
  M_over_Mpl: number;
  /** Transfer function data */
  transfer_functions: TransferFunctionData;
  /** Checks */
  checks: { check_id: string; passed: boolean; detail: string }[];
}

/**
 * Load APS η-gluing spectral flow data
 */
function loadApsEtaGluing(resultsPath: string): { spectralFlow: SpectralFlowEntry[], checks: any[] } {
  const raw = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
  const results = raw.results;
  
  const spectralFlow: SpectralFlowEntry[] = results.table.map((entry: any) => ({
    m: entry.m,
    spin: entry.spin as 'periodic' | 'antiperiodic',
    analytic_sf: entry.analytic_sf,
    numeric_sf: entry.numeric_sf,
    winding_det_u: entry.winding_det_u,
    eta_start_numeric: entry.eta_start_numeric,
    eta_end_numeric: entry.eta_end_numeric,
    theta_start: entry.theta_start,
    theta_end: entry.theta_end
  }));
  
  return {
    spectralFlow,
    checks: raw.checks
  };
}

/**
 * Load effective action R² operator blocks
 */
function loadEffectiveActionR2(resultsPath: string): { 
  blocks: OperatorBlock[], 
  beta_R2_total: number, 
  M_over_Mpl: number 
} {
  const raw = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
  const closure = raw.results.operator_closure_minimal;
  
  const blocks: OperatorBlock[] = closure.blocks.map((b: any) => ({
    name: b.name,
    rank: b.rank,
    statistics: b.statistics as 'boson' | 'ghost',
    E_over_R: parseFloat(b.E_over_R),
    Omega_sq_over_R2: parseFloat(b.Omega_sq_over_R2),
    a2_R2_coeff_curly: parseFloat(b.a2_R2_coeff_curly),
    beta_R2_contribution: parseFloat(b.beta_R2_contribution),
    prefactor: parseFloat(b.prefactor)
  }));
  
  return {
    blocks,
    beta_R2_total: parseFloat(closure.beta_R2_total_after_renorm),
    M_over_Mpl: parseFloat(raw.results.derived.M_over_Mpl)
  };
}

/**
 * Load bounce perturbation transfer functions
 */
function loadBouncePerturbations(resultsPath: string): TransferFunctionData {
  const raw = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
  const results = raw.results;
  
  return {
    k_grid: results.k_grid,
    T_scalar: results.T_scalar,
    T_tensor: results.T_tensor,
    k_bounce_s_est: results.diagnostics.k_bounce_s_est_raw,
    k_bounce_t_est: results.diagnostics.k_bounce_t_est_raw
  };
}

/**
 * Load all spectrum data from TFPT suite artifacts.
 * This function does NOT compute anything - it only reads and structures.
 */
export function loadSpectrum(
  apsEtaPath: string, 
  effectiveR2Path: string,
  bouncePath: string
): SpectrumData {
  const apsData = loadApsEtaGluing(apsEtaPath);
  const r2Data = loadEffectiveActionR2(effectiveR2Path);
  const bounceData = loadBouncePerturbations(bouncePath);
  
  // Compute Δ_Γ = 2π for minimal class (m=1)
  const m1_periodic = apsData.spectralFlow.find(e => e.m === 1 && e.spin === 'periodic');
  const seam_term = m1_periodic ? m1_periodic.theta_end - m1_periodic.theta_start : 0;
  
  return {
    spectral_flow: apsData.spectralFlow,
    seam_term_delta_gamma_2pi: seam_term,
    operator_blocks: r2Data.blocks,
    beta_R2_total: r2Data.beta_R2_total,
    M_over_Mpl: r2Data.M_over_Mpl,
    transfer_functions: bounceData,
    checks: apsData.checks
  };
}
