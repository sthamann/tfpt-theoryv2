/**
 * Torsion Bounds View
 * ===================
 * 
 * Visualizes:
 *   - SME-style torsion bounds (A_T, A_X, A_Y, A_Z)
 *   - TFPT prediction for S_μ
 *   - Prediction/bound ratios
 * 
 * Data source: TFPT_STORE.torsion
 * Paper reference: Section 6, Appendix M
 */

import { TFPT_STORE } from '../state/tfpt_store';

const SVG_WIDTH = 480;
const SVG_HEIGHT = 400;

/**
 * Render torsion bounds as a log-scale bar chart
 */
export function renderTorsionBounds(): string {
  const torsion = TFPT_STORE.torsion;
  if (!torsion) {
    return `<svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}">
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#f00">
        Error: Torsion data not loaded
      </text>
    </svg>`;
  }
  
  const bounds = torsion.inferred_bounds;
  const prediction = torsion.prediction;
  
  // Log scale range (in GeV)
  const allValues = [...bounds.map(b => b.inferred_abs_max_S_mu_GeV), prediction.S_mu_abs_GeV];
  const logMin = Math.floor(Math.log10(Math.min(...allValues))) - 3;
  const logMax = Math.ceil(Math.log10(Math.max(...allValues))) + 2;
  
  const plotWidth = 300;
  const plotHeight = 200;
  const startX = 90;
  const startY = 80;
  
  const mapLog = (val: number) => {
    const logVal = Math.log10(val);
    return startX + ((logVal - logMin) / (logMax - logMin)) * plotWidth;
  };
  
  // Bar chart for bounds
  const barHeight = 35;
  let bars = '';
  for (let i = 0; i < bounds.length; i++) {
    const bound = bounds[i];
    const y = startY + i * (barHeight + 12);
    const x = mapLog(bound.inferred_abs_max_S_mu_GeV);
    const barWidth = x - startX;
    
    // Color based on component
    const colors: Record<string, string> = {
      'A_T': '#E91E63',
      'A_X': '#9C27B0',
      'A_Y': '#3F51B5',
      'A_Z': '#009688'
    };
    const color = colors[bound.label] ?? '#666';
    const exp = Math.round(Math.log10(bound.inferred_abs_max_S_mu_GeV));
    
    bars += `
      <g>
        <text x="${startX - 10}" y="${y + barHeight/2 + 5}" text-anchor="end" fill="#333" font-size="13" font-weight="bold">
          ${bound.label}
        </text>
        <rect x="${startX}" y="${y}" width="${barWidth}" height="${barHeight}" fill="${color}" opacity="0.75" rx="4">
          <title>${bound.label}: |S_μ| &lt; ${bound.inferred_abs_max_S_mu_GeV.toExponential(1)} GeV
${bound.applies_to}</title>
        </rect>
        <text x="${x + 8}" y="${y + barHeight/2 + 5}" fill="#333" font-size="11">
          10^${exp} GeV
        </text>
      </g>
    `;
  }
  
  // TFPT prediction line
  const predX = mapLog(prediction.S_mu_abs_GeV);
  const predExp = Math.round(Math.log10(prediction.S_mu_abs_GeV));
  const predLine = `
    <line x1="${predX}" y1="${startY - 25}" x2="${predX}" y2="${startY + bounds.length * (barHeight + 12) + 10}" 
          stroke="#4CAF50" stroke-width="3" stroke-dasharray="10,5"/>
    <text x="${predX}" y="${startY - 35}" text-anchor="middle" fill="#2e7d32" font-size="11" font-weight="bold">
      TFPT: 10^${predExp} GeV
    </text>
  `;
  
  // Log scale axis
  let axisMarks = '';
  const axisY = startY + bounds.length * (barHeight + 12) + 15;
  for (let log = -45; log <= logMax; log += 5) {
    const x = mapLog(Math.pow(10, log));
    if (x >= startX && x <= startX + plotWidth) {
      axisMarks += `
        <line x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY + 8}" stroke="#666"/>
        <text x="${x}" y="${axisY + 22}" text-anchor="middle" fill="#666" font-size="10">
          10^${log}
        </text>
      `;
    }
  }
  
  // Ratio summary
  const maxRatio = Math.max(...torsion.prediction_to_bound_ratios.map(r => r.ratio));
  
  return `
    <svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}" xmlns="http://www.w3.org/2000/svg">
      <!-- Bounds bars -->
      ${bars}
      
      <!-- TFPT prediction -->
      ${predLine}
      
      <!-- Axis -->
      <line x1="${startX}" y1="${axisY}" x2="${startX + plotWidth}" y2="${axisY}" stroke="#666" stroke-width="1"/>
      ${axisMarks}
      <text x="${startX + plotWidth/2}" y="${axisY + 42}" text-anchor="middle" fill="#333" font-size="11">
        |S_μ| (GeV)
      </text>
      
      <!-- Info box -->
      <g transform="translate(20, ${SVG_HEIGHT - 85})">
        <rect width="${SVG_WIDTH - 40}" height="75" fill="#f5f5f5" stroke="#ddd" rx="6"/>
        <text x="15" y="22" fill="#333" font-size="11" font-weight="500">
          Regime: ${prediction.regime_details.label}
        </text>
        <text x="15" y="42" fill="#666" font-size="10">
          Mapping: b_μ = ${torsion.mapping.k_default} · S_μ (minimal coupling, ξ_tree = ${torsion.mapping.xi_tree})
        </text>
        <text x="15" y="62" fill="${maxRatio < 1e-6 ? '#2e7d32' : '#FF9800'}" font-size="11" font-weight="bold">
          max(pred/bound) = ${maxRatio.toExponential(2)} ${maxRatio < 1e-6 ? '✓ Far below all experimental bounds' : '⚠ Check regime assumptions'}
        </text>
      </g>
    </svg>
  `;
}

/**
 * Generate HTML wrapper
 */
export function getTorsionBoundsHTML(): string {
  return `
    <div id="torsion-bounds" class="view-panel">
      <div class="view-header">
        <h3>Torsion Bounds</h3>
        <p class="view-description">
          <strong>What you see:</strong> Experimental upper bounds on axial torsion S_μ from SME tests (colored bars),
          compared to the TFPT prediction (green dashed line). The prediction is ~10⁻⁴² GeV in the cosmological 
          regime — about 15 orders of magnitude below the tightest lab bounds. This is a falsifiability check.
        </p>
        <span class="view-source">Data: torsion_bounds_mapping/results.json, torsion_regimes.json</span>
      </div>
      <div class="view-content">
        ${renderTorsionBounds()}
      </div>
    </div>
  `;
}
