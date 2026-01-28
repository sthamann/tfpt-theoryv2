/**
 * TFPT Store - Centralized State Container
 * =========================================
 * 
 * All views read from this store. No view touches files directly.
 * The store is populated once at startup by index.ts.
 * 
 * Philosophy:
 *   - Loader reads, Store normalizes, Views show
 *   - Single source of truth for the UI
 *   - If data is missing, the store crashes hard (intentional)
 */

import {
  CascadeData,
  AlphaData,
  SpectrumData,
  TorsionData,
  GlobalData
} from '../loader';

/** Currently active level/stage for synchronized views */
export let ACTIVE_LEVEL = 0;

/** Listeners for ACTIVE_LEVEL changes */
type LevelChangeListener = (level: number) => void;
const levelListeners: LevelChangeListener[] = [];

/**
 * Set the active cascade level and notify all listeners
 */
export function setActiveLevel(level: number): void {
  ACTIVE_LEVEL = level;
  for (const listener of levelListeners) {
    listener(level);
  }
}

/**
 * Subscribe to active level changes
 */
export function onLevelChange(listener: LevelChangeListener): () => void {
  levelListeners.push(listener);
  return () => {
    const idx = levelListeners.indexOf(listener);
    if (idx >= 0) levelListeners.splice(idx, 1);
  };
}

/** The complete TFPT state */
export interface TFPTState {
  /** Cascade/ladder data */
  cascade: CascadeData | null;
  /** Alpha precision data */
  alpha: AlphaData | null;
  /** Spectrum/spectral flow data */
  spectrum: SpectrumData | null;
  /** Torsion bounds data */
  torsion: TorsionData | null;
  /** Global consistency data */
  global: GlobalData | null;
  /** Timestamp when data was loaded */
  loaded_at: Date | null;
  /** Source paths used */
  source_paths: {
    cascade_thresholds?: string;
    cascade_two_loop?: string;
    alpha_precision?: string;
    spectrum_aps?: string;
    spectrum_r2?: string;
    spectrum_bounce?: string;
    torsion_bounds?: string;
    torsion_regimes?: string;
    global_consistency?: string;
  };
}

/**
 * The global TFPT store instance.
 * All views read from this.
 */
export const TFPT_STORE: TFPTState = {
  cascade: null,
  alpha: null,
  spectrum: null,
  torsion: null,
  global: null,
  loaded_at: null,
  source_paths: {}
};

/**
 * Check if the store is fully populated
 */
export function isStoreReady(): boolean {
  return (
    TFPT_STORE.cascade !== null &&
    TFPT_STORE.alpha !== null &&
    TFPT_STORE.spectrum !== null &&
    TFPT_STORE.torsion !== null &&
    TFPT_STORE.global !== null
  );
}

/**
 * Get a summary of loaded data for diagnostics
 */
export function getStoreSummary(): object {
  return {
    ready: isStoreReady(),
    loaded_at: TFPT_STORE.loaded_at?.toISOString() ?? null,
    cascade_stages: TFPT_STORE.cascade?.stages.length ?? 0,
    cascade_crossings: TFPT_STORE.cascade?.crossings.length ?? 0,
    alpha_k_points: TFPT_STORE.alpha?.k_sensitivity.length ?? 0,
    spectral_flow_entries: TFPT_STORE.spectrum?.spectral_flow.length ?? 0,
    operator_blocks: TFPT_STORE.spectrum?.operator_blocks.length ?? 0,
    torsion_bounds: TFPT_STORE.torsion?.inferred_bounds.length ?? 0,
    global_terms: TFPT_STORE.global?.terms.length ?? 0,
    source_paths: TFPT_STORE.source_paths
  };
}

/**
 * Serialize store to JSON for API responses
 */
export function serializeStore(): string {
  return JSON.stringify({
    cascade: TFPT_STORE.cascade,
    alpha: TFPT_STORE.alpha,
    spectrum: TFPT_STORE.spectrum,
    torsion: TFPT_STORE.torsion,
    global: TFPT_STORE.global,
    active_level: ACTIVE_LEVEL,
    loaded_at: TFPT_STORE.loaded_at?.toISOString() ?? null
  }, null, 2);
}
