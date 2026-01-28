/**
 * Cascade / Ladder Loader
 * =======================
 * Loads cascade stages from RG thresholds and two-loop fingerprints.
 * 
 * Data sources:
 *   - tfpt_suite/data/rge_thresholds_v25.json (threshold scales)
 *   - out/two_loop_rg_fingerprints/results.json (α3 running, crossings)
 * 
 * Paper reference: Section 4 (E8 Cascade), Table 2 (Threshold scales)
 */

import * as fs from 'fs';
import * as path from 'path';

/** A single cascade stage with symmetry breaking info */
export interface CascadeStage {
  level: number;
  name: string;
  scale_GeV: number;
  log10_scale: number;
  group_before: string;
  group_after: string;
  description: string;
}

/** Crossing point where α3 equals a TFPT invariant */
export interface AlphaCrossing {
  target_name: string;
  target_value: number;
  mu_GeV: number;
  log10_mu: number;
  alpha3_at_crossing: number;
}

/** Full cascade data structure */
export interface CascadeData {
  stages: CascadeStage[];
  crossings: AlphaCrossing[];
  alpha3_range: [number, number];
  is_monotone_decreasing: boolean;
}

/**
 * Load threshold scales from rge_thresholds_v25.json
 */
function loadThresholds(thresholdsPath: string): Map<string, number> {
  const raw = JSON.parse(fs.readFileSync(thresholdsPath, 'utf8'));
  const thresholds = new Map<string, number>();
  
  for (const [name, value] of Object.entries(raw.thresholds_GeV)) {
    thresholds.set(name, value as number);
  }
  
  return thresholds;
}

/**
 * Load two-loop RG fingerprints results
 */
function loadTwoLoopResults(resultsPath: string): any {
  return JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
}

/**
 * Load cascade data from TFPT suite artifacts.
 * This function does NOT compute anything - it only reads and structures.
 */
export function loadCascade(thresholdsPath: string, twoLoopPath: string): CascadeData {
  const thresholds = loadThresholds(thresholdsPath);
  const twoLoop = loadTwoLoopResults(twoLoopPath);
  
  // Build cascade stages from thresholds
  // Order: E8 → G8 → SM × U(1)_PQ → SM
  const stages: CascadeStage[] = [];
  let level = 0;
  
  // E8 breaking stages based on the thresholds
  const stageDefinitions = [
    {
      name: "MNR3",
      group_before: "E8 → SO(10) × U(1)³",
      group_after: "E8 → SO(10) × U(1)² (ν_R3 active)",
      description: "Third right-handed neutrino threshold"
    },
    {
      name: "MNR2", 
      group_before: "E8 → SO(10) × U(1)² (ν_R3)",
      group_after: "E8 → SO(10) × U(1) (ν_R2,3 active)",
      description: "Second right-handed neutrino threshold"
    },
    {
      name: "MNR1",
      group_before: "E8 → SO(10) × U(1) (ν_R2,3)",
      group_after: "SM × U(1)_PQ (all ν_R active)",
      description: "First right-handed neutrino threshold"
    },
    {
      name: "MPhi",
      group_before: "SM × U(1)_PQ",
      group_after: "SM × U(1)_PQ (Φ active)",
      description: "PQ scalar threshold (f_a scale)"
    },
    {
      name: "MG8",
      group_before: "SM × U(1)_PQ × G8",
      group_after: "SM × U(1)_PQ",
      description: "Color octet scalar threshold"
    },
    {
      name: "MSigma",
      group_before: "SM × Σ sector",
      group_after: "SM",
      description: "Σ multiplet threshold"
    }
  ];
  
  // Sort by scale (descending)
  const sortedStages = stageDefinitions
    .filter(s => thresholds.has(s.name))
    .sort((a, b) => (thresholds.get(b.name) ?? 0) - (thresholds.get(a.name) ?? 0));
  
  for (const stageDef of sortedStages) {
    const scale = thresholds.get(stageDef.name)!;
    stages.push({
      level: level++,
      name: stageDef.name,
      scale_GeV: scale,
      log10_scale: Math.log10(scale),
      group_before: stageDef.group_before,
      group_after: stageDef.group_after,
      description: stageDef.description
    });
  }
  
  // Extract crossings from two-loop results
  const crossings: AlphaCrossing[] = [];
  
  const alpha3Results = twoLoop.results?.alpha3;
  if (alpha3Results) {
    // varphi0 crossing
    if (alpha3Results.crossing_varphi0) {
      const cv = alpha3Results.crossing_varphi0;
      crossings.push({
        target_name: "φ₀",
        target_value: cv.alpha3,
        mu_GeV: cv.mu_GeV,
        log10_mu: Math.log10(cv.mu_GeV),
        alpha3_at_crossing: cv.alpha3
      });
    }
    
    // c3 crossing
    if (alpha3Results.crossing_c3) {
      const cc = alpha3Results.crossing_c3;
      crossings.push({
        target_name: "c₃",
        target_value: cc.alpha3,
        mu_GeV: cc.mu_GeV,
        log10_mu: Math.log10(cc.mu_GeV),
        alpha3_at_crossing: cc.alpha3
      });
    }
  }
  
  return {
    stages,
    crossings,
    alpha3_range: alpha3Results?.range ?? [0, 0],
    is_monotone_decreasing: alpha3Results?.is_monotone_decreasing ?? false
  };
}
