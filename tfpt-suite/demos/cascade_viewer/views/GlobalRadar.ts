/**
 * Global Radar - χ² Scorecard
 * ===========================
 * 
 * Visualizes:
 *   - Per-observable χ² contributions
 *   - z-scores (deviations from reference)
 *   - Global fit quality (p-value)
 * 
 * Data source: TFPT_STORE.global
 * Paper reference: Table 1, Section 7
 */

import { TFPT_STORE } from '../state/tfpt_store';

const SVG_WIDTH = 500;
const SVG_HEIGHT = 480;
const RADAR_CENTER_X = 160;
const RADAR_CENTER_Y = 200;
const RADAR_RADIUS = 130;

/**
 * Render radar chart for z-scores
 */
function renderRadarChart(): string {
  const global = TFPT_STORE.global;
  if (!global) return '';
  
  const terms = global.terms;
  const n = terms.length;
  if (n === 0) return '';
  
  const angleStep = (2 * Math.PI) / n;
  
  // Draw axes
  let axes = '';
  let labels = '';
  for (let i = 0; i < n; i++) {
    const angle = -Math.PI / 2 + i * angleStep;
    const x = RADAR_CENTER_X + RADAR_RADIUS * Math.cos(angle);
    const y = RADAR_CENTER_Y + RADAR_RADIUS * Math.sin(angle);
    const labelX = RADAR_CENTER_X + (RADAR_RADIUS + 35) * Math.cos(angle);
    const labelY = RADAR_CENTER_Y + (RADAR_RADIUS + 35) * Math.sin(angle);
    
    axes += `<line x1="${RADAR_CENTER_X}" y1="${RADAR_CENTER_Y}" x2="${x}" y2="${y}" stroke="#ddd" stroke-width="1"/>`;
    
    // Label with observable name
    const term = terms[i];
    const displayName = term.name
      .replace('alpha_bar5_inv_MZ', 'α̅⁵(MZ)')
      .replace('cabibbo_lambda', 'λ_Cabibbo')
      .replace('beta_deg', 'β_CMB')
      .replace('n_s', 'n_s')
      .replace('A_s', 'A_s');
    labels += `<text x="${labelX}" y="${labelY + 4}" text-anchor="middle" fill="#333" font-size="11" font-weight="500">${displayName}</text>`;
  }
  
  // Draw concentric circles for |z| = 1, 2, 3
  let circles = '';
  for (const zLevel of [1, 2, 3]) {
    const r = (zLevel / 3) * RADAR_RADIUS;
    circles += `<circle cx="${RADAR_CENTER_X}" cy="${RADAR_CENTER_Y}" r="${r}" fill="none" stroke="#eee" stroke-width="1"/>`;
    if (zLevel < 3) {
      circles += `<text x="${RADAR_CENTER_X + r + 5}" y="${RADAR_CENTER_Y - 5}" fill="#aaa" font-size="9">${zLevel}σ</text>`;
    }
  }
  circles += `<text x="${RADAR_CENTER_X + RADAR_RADIUS + 5}" y="${RADAR_CENTER_Y - 5}" fill="#aaa" font-size="9">3σ</text>`;
  
  // Build radar polygon from z-scores
  let points = '';
  for (let i = 0; i < n; i++) {
    const angle = -Math.PI / 2 + i * angleStep;
    const zMag = Math.min(Math.abs(terms[i].z_total), 3); // Cap at 3σ for display
    const r = (zMag / 3) * RADAR_RADIUS;
    const x = RADAR_CENTER_X + r * Math.cos(angle);
    const y = RADAR_CENTER_Y + r * Math.sin(angle);
    points += `${x},${y} `;
  }
  
  // Data points
  let dataPoints = '';
  for (let i = 0; i < n; i++) {
    const angle = -Math.PI / 2 + i * angleStep;
    const term = terms[i];
    const zMag = Math.min(Math.abs(term.z_total), 3);
    const r = (zMag / 3) * RADAR_RADIUS;
    const x = RADAR_CENTER_X + r * Math.cos(angle);
    const y = RADAR_CENTER_Y + r * Math.sin(angle);
    
    // Color based on z-score magnitude
    let color = '#4CAF50'; // Green for |z| < 1
    if (Math.abs(term.z_total) >= 1.5) color = '#FF9800'; // Orange for |z| >= 1.5
    if (Math.abs(term.z_total) >= 2.5) color = '#f44336'; // Red for |z| >= 2.5
    
    dataPoints += `
      <circle cx="${x}" cy="${y}" r="8" fill="${color}" stroke="#fff" stroke-width="2">
        <title>${term.name}
z-score: ${term.z_total.toFixed(3)}
χ² contribution: ${term.chi2_total.toFixed(4)}
Prediction: ${term.pred.toPrecision(6)}
Reference: ${term.mean} ± ${term.sigma_exp}</title>
      </circle>
    `;
  }
  
  return `
    <g class="radar-chart">
      ${circles}
      ${axes}
      <polygon points="${points}" fill="#2196F3" fill-opacity="0.2" stroke="#2196F3" stroke-width="2"/>
      ${dataPoints}
      ${labels}
    </g>
  `;
}

/**
 * Render χ² bar chart
 */
function renderChi2Bars(): string {
  const global = TFPT_STORE.global;
  if (!global) return '';
  
  const terms = global.terms;
  const barHeight = 24;
  const maxChi2 = Math.max(...terms.map(t => t.chi2_total), 1);
  const barWidth = 140;
  const startX = 340;
  const startY = 80;
  
  let bars = '';
  for (let i = 0; i < terms.length; i++) {
    const term = terms[i];
    const y = startY + i * (barHeight + 12);
    const width = (term.chi2_total / maxChi2) * barWidth;
    
    // Color based on contribution
    const frac = term.chi2_total / global.totals.chi2_total;
    const hue = 120 - frac * 120; // Green to red
    const color = `hsl(${hue}, 70%, 50%)`;
    
    const displayName = term.name
      .replace('alpha_bar5_inv_MZ', 'α̅⁵(MZ)')
      .replace('cabibbo_lambda', 'λ_C')
      .replace('beta_deg', 'β_CMB')
      .replace('n_s', 'n_s')
      .replace('A_s', 'A_s');
    
    bars += `
      <g>
        <rect x="${startX}" y="${y}" width="${width}" height="${barHeight}" fill="${color}" rx="3"/>
        <text x="${startX + width + 8}" y="${y + barHeight/2 + 5}" fill="#333" font-size="11">
          ${term.chi2_total.toFixed(2)}
        </text>
        <text x="${startX - 8}" y="${y + barHeight/2 + 5}" text-anchor="end" fill="#333" font-size="10">
          ${displayName}
        </text>
      </g>
    `;
  }
  
  return `
    <g class="chi2-bars">
      <text x="${startX + barWidth/2}" y="${startY - 25}" text-anchor="middle" fill="#333" font-size="12" font-weight="600">
        χ² Contributions
      </text>
      ${bars}
    </g>
  `;
}

/**
 * Render global radar view
 */
export function renderGlobalRadar(): string {
  const global = TFPT_STORE.global;
  if (!global) {
    return `<svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}">
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#f00">
        Error: Global data not loaded
      </text>
    </svg>`;
  }
  
  const totals = global.totals;
  const pValue = totals.p_value;
  const pValueColor = pValue > 0.05 ? '#2e7d32' : (pValue > 0.01 ? '#FF9800' : '#f44336');
  
  return `
    <svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}" xmlns="http://www.w3.org/2000/svg">
      <!-- Radar chart -->
      ${renderRadarChart()}
      
      <!-- Chi² bars -->
      ${renderChi2Bars()}
      
      <!-- Summary box -->
      <g transform="translate(20, ${SVG_HEIGHT - 100})">
        <rect width="200" height="90" fill="#f9f9f9" stroke="#ddd" rx="6"/>
        <text x="15" y="22" fill="#333" font-size="12" font-weight="600">Global Fit Summary</text>
        <text x="15" y="42" fill="#333" font-size="11">χ² = ${totals.chi2_total.toFixed(2)} / ${totals.dof} dof</text>
        <text x="15" y="60" fill="${pValueColor}" font-size="13" font-weight="bold">
          p-value = ${(pValue * 100).toFixed(1)}% ${pValue > 0.05 ? '✓' : '⚠'}
        </text>
        <text x="15" y="80" fill="#666" font-size="10">
          (excl. α: χ²=${totals.chi2_excluding_alpha.toFixed(2)}, p=${(totals.p_value_excluding_alpha * 100).toFixed(1)}%)
        </text>
      </g>
      
      <!-- Mode badge -->
      <g transform="translate(${SVG_WIDTH - 100}, ${SVG_HEIGHT - 100})">
        <rect width="80" height="30" fill="#e3f2fd" stroke="#1565c0" rx="15"/>
        <text x="40" y="20" text-anchor="middle" fill="#1565c0" font-size="11" font-weight="600">${global.mode}</text>
      </g>
    </svg>
  `;
}

/**
 * Generate HTML wrapper
 */
export function getGlobalRadarHTML(): string {
  return `
    <div id="global-radar" class="view-panel">
      <div class="view-header">
        <h3>Global Consistency Radar</h3>
        <p class="view-description">
          <strong>What you see:</strong> A multi-observable consistency check. The radar shows z-scores (deviations in σ)
          for each observable — smaller means better agreement. The bar chart shows χ² contributions.
          With p ≈ 80%, the global fit is statistically acceptable. Points within the 1σ circle are green.
        </p>
        <span class="view-source">Data: global_consistency_test/results.json</span>
      </div>
      <div class="view-content">
        ${renderGlobalRadar()}
      </div>
    </div>
  `;
}
