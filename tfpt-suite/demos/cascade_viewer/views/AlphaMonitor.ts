/**
 * Alpha Monitor - CFE Fixpoint Visualization
 * ==========================================
 * 
 * Shows the α(0) self-consistency iteration:
 *   - k-sensitivity curve (α vs k)
 *   - ppm deviation from CODATA
 *   - Convergence status
 * 
 * Data source: TFPT_STORE.alpha
 * Paper reference: Section 3.2, Eq. (3.7) - CFE cubic equation
 * 
 * CFE: α³ - 2c₃³α² - 8b₁c₃⁶ ln(1/φ) = 0
 * Backreaction: φ(α) = φ_tree + δ_top exp(-kα)
 */

import { TFPT_STORE } from '../state/tfpt_store';

const SVG_WIDTH = 520;
const SVG_HEIGHT = 520;
const PLOT_MARGIN = { top: 60, right: 40, bottom: 70, left: 80 };

/**
 * Map data value to SVG coordinate
 */
function mapX(k: number, kMin: number, kMax: number): number {
  const plotWidth = SVG_WIDTH - PLOT_MARGIN.left - PLOT_MARGIN.right;
  return PLOT_MARGIN.left + ((k - kMin) / (kMax - kMin)) * plotWidth;
}

function mapY(ppm: number, ppmMin: number, ppmMax: number): number {
  const plotHeight = SVG_HEIGHT - PLOT_MARGIN.top - PLOT_MARGIN.bottom;
  return PLOT_MARGIN.top + ((ppmMax - ppm) / (ppmMax - ppmMin)) * plotHeight;
}

/**
 * Render the alpha monitor as SVG
 */
export function renderAlphaMonitor(): string {
  const alpha = TFPT_STORE.alpha;
  if (!alpha) {
    return `<svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}">
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#f00">
        Error: Alpha data not loaded
      </text>
    </svg>`;
  }
  
  const kSensitivity = alpha.k_sensitivity;
  
  // Calculate ranges with some padding
  const kValues = kSensitivity.map(p => p.k);
  const ppmValues = kSensitivity.map(p => p.ppm_vs_codata);
  
  const kMin = Math.min(...kValues) - 0.3;
  const kMax = Math.max(...kValues) + 0.3;
  const ppmMin = Math.min(...ppmValues, -2) - 0.5;
  const ppmMax = Math.max(...ppmValues, 4) + 0.5;
  
  // Build data points path
  let pathD = '';
  for (let i = 0; i < kSensitivity.length; i++) {
    const pt = kSensitivity[i];
    const x = mapX(pt.k, kMin, kMax);
    const y = mapY(pt.ppm_vs_codata, ppmMin, ppmMax);
    pathD += i === 0 ? `M ${x} ${y}` : ` L ${x} ${y}`;
  }
  
  // Build data point markers
  let dataPoints = '';
  for (const pt of kSensitivity) {
    const x = mapX(pt.k, kMin, kMax);
    const y = mapY(pt.ppm_vs_codata, ppmMin, ppmMax);
    const isK2 = Math.abs(pt.k - 2.0) < 0.01;
    const color = isK2 ? '#4CAF50' : '#2196F3';
    const radius = isK2 ? 10 : 6;
    
    dataPoints += `
      <circle cx="${x}" cy="${y}" r="${radius}" fill="${color}" stroke="#fff" stroke-width="2">
        <title>k = ${pt.k.toFixed(1)}
α⁻¹ = ${pt.alpha_inv_numeric.toFixed(8)}
Δppm = ${pt.ppm_vs_codata.toFixed(4)}</title>
      </circle>
    `;
  }
  
  // Zero line
  const zeroY = mapY(0, ppmMin, ppmMax);
  const zeroLine = `<line x1="${PLOT_MARGIN.left}" y1="${zeroY}" x2="${SVG_WIDTH - PLOT_MARGIN.right}" y2="${zeroY}" 
                          stroke="#4CAF50" stroke-width="2" stroke-dasharray="8,4"/>`;
  
  // CODATA reference band (±0.1 ppm)
  const bandTop = mapY(0.1, ppmMin, ppmMax);
  const bandBottom = mapY(-0.1, ppmMin, ppmMax);
  const codataBand = `<rect x="${PLOT_MARGIN.left}" y="${bandTop}" 
                            width="${SVG_WIDTH - PLOT_MARGIN.left - PLOT_MARGIN.right}" 
                            height="${bandBottom - bandTop}" 
                            fill="#4CAF50" opacity="0.15"/>`;
  
  // Grid lines
  let gridLines = '';
  for (let k = 0; k <= 3; k += 0.5) {
    if (k >= kMin && k <= kMax) {
      const x = mapX(k, kMin, kMax);
      gridLines += `<line x1="${x}" y1="${PLOT_MARGIN.top}" x2="${x}" y2="${SVG_HEIGHT - PLOT_MARGIN.bottom}" stroke="#eee" stroke-width="1"/>`;
      gridLines += `<text x="${x}" y="${SVG_HEIGHT - PLOT_MARGIN.bottom + 20}" text-anchor="middle" fill="#666" font-size="11">${k.toFixed(1)}</text>`;
    }
  }
  
  for (let ppm = Math.ceil(ppmMin); ppm <= Math.floor(ppmMax); ppm++) {
    const y = mapY(ppm, ppmMin, ppmMax);
    gridLines += `<line x1="${PLOT_MARGIN.left}" y1="${y}" x2="${SVG_WIDTH - PLOT_MARGIN.right}" y2="${y}" stroke="#eee" stroke-width="1"/>`;
    gridLines += `<text x="${PLOT_MARGIN.left - 10}" y="${y + 4}" text-anchor="end" fill="#666" font-size="11">${ppm}</text>`;
  }
  
  // Self-consistent point annotation
  const scX = mapX(alpha.self_consistent.k, kMin, kMax);
  const scY = mapY(alpha.self_consistent.ppm_vs_codata, ppmMin, ppmMax);
  const scAnnotation = `
    <line x1="${scX}" y1="${scY + 15}" x2="${scX}" y2="${SVG_HEIGHT - PLOT_MARGIN.bottom}" stroke="#4CAF50" stroke-width="1" stroke-dasharray="4,2"/>
    <text x="${scX}" y="${SVG_HEIGHT - PLOT_MARGIN.bottom + 35}" text-anchor="middle" fill="#2e7d32" font-size="10" font-weight="bold">
      k = 2 (paper value)
    </text>
  `;
  
  return `
    <svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}" xmlns="http://www.w3.org/2000/svg">
      <!-- Plot area background -->
      <rect x="${PLOT_MARGIN.left}" y="${PLOT_MARGIN.top}" 
            width="${SVG_WIDTH - PLOT_MARGIN.left - PLOT_MARGIN.right}" 
            height="${SVG_HEIGHT - PLOT_MARGIN.top - PLOT_MARGIN.bottom}" 
            fill="#fafafa" stroke="#ddd"/>
      
      <!-- Grid -->
      ${gridLines}
      
      <!-- CODATA band -->
      ${codataBand}
      
      <!-- Zero line -->
      ${zeroLine}
      
      <!-- Data curve -->
      <path d="${pathD}" fill="none" stroke="#2196F3" stroke-width="3"/>
      
      <!-- Data points -->
      ${dataPoints}
      
      <!-- Self-consistent annotation -->
      ${scAnnotation}
      
      <!-- Axes labels -->
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT - 15}" text-anchor="middle" fill="#333" font-size="12" font-weight="500">
        Backreaction exponent k
      </text>
      <text x="20" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#333" font-size="12" font-weight="500" transform="rotate(-90, 20, ${SVG_HEIGHT/2})">
        Δppm vs CODATA
      </text>
      
      <!-- Legend -->
      <g transform="translate(${SVG_WIDTH - 180}, ${PLOT_MARGIN.top + 10})">
        <rect width="165" height="85" fill="#fff" stroke="#ddd" rx="4"/>
        <circle cx="15" cy="18" r="5" fill="#2196F3"/>
        <text x="28" y="22" font-size="10" fill="#333">Iteration points</text>
        <circle cx="15" cy="40" r="8" fill="#4CAF50"/>
        <text x="28" y="44" font-size="10" fill="#333">k=2 (self-consistent)</text>
        <line x1="8" y1="62" x2="22" y2="62" stroke="#4CAF50" stroke-dasharray="8,4" stroke-width="2"/>
        <text x="28" y="66" font-size="10" fill="#333">CODATA reference</text>
        <rect x="8" y="74" width="14" height="6" fill="#4CAF50" opacity="0.3"/>
        <text x="28" y="80" font-size="9" fill="#666">±0.1 ppm band</text>
      </g>
      
      <!-- Status box: k=2 baseline -->
      <g transform="translate(${PLOT_MARGIN.left + 10}, ${PLOT_MARGIN.top + 10})">
        <rect width="185" height="62" fill="#fff" stroke="#ddd" rx="4"/>
        <text x="10" y="20" font-size="11" fill="#333" font-weight="500">
          α⁻¹(k=2) = ${alpha.self_consistent.alpha_inv_numeric.toFixed(8)}
        </text>
        <text x="10" y="38" font-size="11" fill="${alpha.self_consistent.converged ? '#4CAF50' : '#f00'}">
          ${alpha.self_consistent.converged ? '✓' : '✗'} Converged in ${alpha.self_consistent.iterations} iterations
        </text>
        <text x="10" y="54" font-size="10" fill="#666">
          Δppm = ${alpha.self_consistent.ppm_vs_codata.toFixed(4)} (baseline, no δ₂)
        </text>
      </g>
      
      ${renderDefectPartitionG5(alpha)}
    </svg>
  `;
}

/**
 * Render the g=5 defect partition section
 */
function renderDefectPartitionG5(alpha: any): string {
  const g5 = alpha.defect_g5;
  if (!g5) return '';
  
  const startY = SVG_HEIGHT - 115;
  
  // Build negative control comparison
  let controlBars = '';
  const maxAbsZ = Math.max(...g5.negative_control.map((c: any) => Math.abs(c.z)));
  const barWidth = 60;
  
  for (let i = 0; i < g5.negative_control.length; i++) {
    const c = g5.negative_control[i];
    const x = 320 + i * 65;
    const height = Math.min((Math.abs(c.z) / maxAbsZ) * 45, 45);
    const color = c.g === 5 ? '#4CAF50' : '#f44336';
    
    controlBars += `
      <rect x="${x}" y="${startY + 50 - height}" width="50" height="${height}" fill="${color}" opacity="0.7" rx="3">
        <title>g=${c.g}: z=${c.z.toFixed(2)}σ</title>
      </rect>
      <text x="${x + 25}" y="${startY + 68}" text-anchor="middle" fill="#333" font-size="10">g=${c.g}</text>
      <text x="${x + 25}" y="${startY + 82}" text-anchor="middle" fill="#666" font-size="9">${Math.abs(c.z).toFixed(1)}σ</text>
    `;
  }
  
  return `
    <!-- Defect partition g=5 section (theoryv3) -->
    <rect x="10" y="${startY - 15}" width="${SVG_WIDTH - 20}" height="110" fill="#e8f5e9" stroke="#4CAF50" stroke-width="1" rx="6"/>
    
    <text x="20" y="${startY + 5}" font-size="12" font-weight="600" fill="#2e7d32">
      Defect Partition: δ₂ = (g/4)·δ_top² with g = ${g5.g_value}
    </text>
    
    <text x="20" y="${startY + 25}" font-size="11" fill="#333">
      α⁻¹(0) = ${g5.alpha_inv_pred.toFixed(11)}
    </text>
    <text x="20" y="${startY + 43}" font-size="11" fill="#333">
      CODATA: ${g5.alpha_inv_ref} ± ${g5.alpha_inv_sigma}
    </text>
    <text x="20" y="${startY + 61}" font-size="12" font-weight="bold" fill="${Math.abs(g5.z_score) < 2 ? '#2e7d32' : '#FF9800'}">
      z = ${g5.z_score.toFixed(2)}σ ${Math.abs(g5.z_score) < 2 ? '✓ within 2σ!' : ''}
    </text>
    
    <text x="20" y="${startY + 82}" font-size="10" fill="#666">
      δ₂/δ_top² = ${g5.delta2_over_delta_top2} (discrete, not fitted)
    </text>
    
    <!-- Negative control -->
    <text x="320" y="${startY + 5}" font-size="10" fill="#666">Negative control |z|:</text>
    ${controlBars}
  `;
}

/**
 * Generate HTML wrapper
 */
export function getAlphaMonitorHTML(): string {
  return `
    <div id="alpha-monitor" class="view-panel">
      <div class="view-header">
        <h3>α Homeostasis Monitor</h3>
        <p class="view-description">
          <strong>What you see:</strong> The fine-structure constant α(0) computed from the CFE (Cubic Fixed-point Equation).
          <em>Top:</em> k-sensitivity curve — at k=2 (Möbius double cover), the baseline prediction is within 0.04 ppm.
          <em>Bottom (green):</em> The defect partition with g=5 adds a second-order correction δ₂ = (5/4)δ_top², 
          bringing the final prediction to within <strong>1.86σ</strong> of CODATA 2022 — with no fit parameters!
        </p>
        <span class="view-source">Data: alpha_precision_audit, theoryv3/defect_partition_g5_audit</span>
      </div>
      <div class="view-content">
        ${renderAlphaMonitor()}
      </div>
    </div>
  `;
}
