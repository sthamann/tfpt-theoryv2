/**
 * Ladder View - E8 Cascade Visualization
 * =======================================
 * 
 * Renders the cascade stages as a vertical ladder showing:
 *   - RG threshold scales (MSigma → MG8 → MPhi → MNR1-3)
 *   - Energy scale annotations (log₁₀μ in GeV)
 *   - α₃ crossings with TFPT invariants (φ₀, c₃)
 * 
 * Data source: TFPT_STORE.cascade
 * Paper reference: Section 4 (E8 Cascade), Table 2 (Threshold scales)
 */

import { TFPT_STORE, ACTIVE_LEVEL, setActiveLevel } from '../state/tfpt_store';
import { CascadeStage, AlphaCrossing } from '../loader';

const SVG_WIDTH = 480;
const SVG_HEIGHT = 520;
const LADDER_X = 120;
const STAGE_HEIGHT = 65;
const STAGE_RADIUS = 22;

/**
 * Generate SVG for a single cascade stage
 */
function renderStage(stage: CascadeStage, index: number, yPos: number, isActive: boolean): string {
  const fillColor = isActive ? '#4CAF50' : '#2196F3';
  const textColor = '#fff';
  
  // Scale in scientific notation
  const exp = Math.floor(stage.log10_scale);
  const mantissa = Math.pow(10, stage.log10_scale - exp);
  const scaleStr = mantissa > 1.5 ? `${mantissa.toFixed(1)}×10^${exp}` : `10^${exp.toFixed(1)}`;
  
  return `
    <g class="cascade-stage" data-level="${index}" onclick="window.setActiveLevel(${index})">
      <!-- Stage circle -->
      <circle cx="${LADDER_X}" cy="${yPos}" r="${STAGE_RADIUS}" fill="${fillColor}" stroke="#fff" stroke-width="2">
        <title>Click to select this stage</title>
      </circle>
      
      <!-- Stage name inside circle -->
      <text x="${LADDER_X}" y="${yPos + 5}" text-anchor="middle" fill="${textColor}" font-size="11" font-weight="bold">
        ${stage.name}
      </text>
      
      <!-- Scale annotation (left side) -->
      <text x="${LADDER_X - 45}" y="${yPos + 4}" text-anchor="end" fill="#555" font-size="11" font-weight="500">
        ${scaleStr} GeV
      </text>
      
      <!-- Description (right side) -->
      <text x="${LADDER_X + 45}" y="${yPos - 5}" text-anchor="start" fill="#666" font-size="10">
        ${stage.description}
      </text>
      <text x="${LADDER_X + 45}" y="${yPos + 10}" text-anchor="start" fill="#333" font-size="10" font-weight="500">
        → ${getShortGroupName(stage.group_after)}
      </text>
    </g>
  `;
}

/**
 * Get shortened group name for display
 */
function getShortGroupName(group: string): string {
  if (group.includes('SM × U(1)_PQ')) return 'SM × U(1)_PQ';
  if (group.includes('all ν_R')) return 'SM × U(1)_PQ (ν_R active)';
  if (group.includes('Φ active')) return 'SM × U(1)_PQ × Φ';
  if (group.includes('SM')) return 'SM';
  return group.split('→').pop()?.trim() ?? group;
}

/**
 * Generate SVG for a crossing marker
 */
function renderCrossing(crossing: AlphaCrossing, yPos: number): string {
  const exp = Math.floor(Math.log10(crossing.mu_GeV));
  
  return `
    <g class="crossing-marker">
      <line x1="${LADDER_X - 25}" y1="${yPos}" x2="${LADDER_X + 25}" y2="${yPos}" 
            stroke="#FF9800" stroke-width="2" stroke-dasharray="6,3"/>
      <circle cx="${LADDER_X}" cy="${yPos}" r="6" fill="#FF9800" stroke="#fff" stroke-width="1"/>
      <text x="${LADDER_X + 35}" y="${yPos + 4}" fill="#e65100" font-size="10" font-weight="bold">
        α₃ = ${crossing.target_name} (10^${exp.toFixed(1)} GeV)
      </text>
    </g>
  `;
}

/**
 * Render the complete ladder view as SVG
 */
export function renderLadderView(): string {
  const cascade = TFPT_STORE.cascade;
  if (!cascade) {
    return `<svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}">
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#f00">
        Error: Cascade data not loaded
      </text>
    </svg>`;
  }
  
  const stages = cascade.stages;
  const crossings = cascade.crossings;
  
  // Calculate positions
  const startY = 80;
  const stagePositions = stages.map((s, i) => ({
    stage: s,
    index: i,
    y: startY + i * STAGE_HEIGHT
  }));
  
  // Build ladder lines between stages
  let ladderLines = '';
  for (let i = 0; i < stagePositions.length - 1; i++) {
    const y1 = stagePositions[i].y + STAGE_RADIUS;
    const y2 = stagePositions[i + 1].y - STAGE_RADIUS;
    ladderLines += `<line x1="${LADDER_X}" y1="${y1}" x2="${LADDER_X}" y2="${y2}" stroke="#ccc" stroke-width="3"/>`;
  }
  
  // Build stage nodes
  let stageNodes = '';
  for (const { stage, index, y } of stagePositions) {
    stageNodes += renderStage(stage, index, y, index === ACTIVE_LEVEL);
  }
  
  // Build crossing markers (interpolate Y position based on log scale)
  let crossingMarkers = '';
  if (stages.length >= 2) {
    const minLog = stages[stages.length - 1].log10_scale;
    const maxLog = stages[0].log10_scale;
    const minY = startY + (stages.length - 1) * STAGE_HEIGHT;
    const maxY = startY;
    
    for (const crossing of crossings) {
      const frac = (crossing.log10_mu - minLog) / (maxLog - minLog);
      const y = minY - frac * (minY - maxY);
      if (y > startY - 20 && y < minY + 20) {
        crossingMarkers += renderCrossing(crossing, y);
      }
    }
  }
  
  const totalHeight = startY + stages.length * STAGE_HEIGHT + 40;
  
  return `
    <svg width="${SVG_WIDTH}" height="${Math.max(SVG_HEIGHT, totalHeight)}" xmlns="http://www.w3.org/2000/svg">
      <defs>
        <linearGradient id="ladderGrad" x1="0%" y1="0%" x2="0%" y2="100%">
          <stop offset="0%" style="stop-color:#1a237e;stop-opacity:0.08" />
          <stop offset="100%" style="stop-color:#1a237e;stop-opacity:0.15" />
        </linearGradient>
      </defs>
      
      <!-- Background -->
      <rect x="${LADDER_X - 35}" y="60" width="70" height="${stages.length * STAGE_HEIGHT + 20}" 
            fill="url(#ladderGrad)" rx="10"/>
      
      <!-- Ladder lines -->
      ${ladderLines}
      
      <!-- Crossing markers -->
      ${crossingMarkers}
      
      <!-- Stage nodes -->
      ${stageNodes}
      
      <!-- Footer info -->
      <text x="${SVG_WIDTH/2}" y="${totalHeight - 25}" text-anchor="middle" font-size="10" fill="#666">
        α₃ range: [${cascade.alpha3_range[0].toFixed(4)}, ${cascade.alpha3_range[1].toFixed(4)}]
      </text>
      <text x="${SVG_WIDTH/2}" y="${totalHeight - 8}" text-anchor="middle" font-size="10" fill="${cascade.is_monotone_decreasing ? '#4CAF50' : '#f00'}">
        ${cascade.is_monotone_decreasing ? '✓ Monotone decreasing (asymptotic freedom)' : '✗ Not monotone'}
      </text>
    </svg>
  `;
}

/**
 * Generate HTML wrapper with the ladder SVG
 */
export function getLadderViewHTML(): string {
  return `
    <div id="ladder-view" class="view-panel">
      <div class="view-header">
        <h3>E8 Cascade Ladder</h3>
        <p class="view-description">
          <strong>What you see:</strong> The symmetry breaking cascade from E8 down to the Standard Model.
          Each node represents a threshold scale where heavy particles decouple.
          Orange markers show where α₃(μ) crosses TFPT invariants φ₀ and c₃.
          <strong>Click a stage</strong> to highlight it across all views.
        </p>
        <span class="view-source">Data: rge_thresholds_v25.json, two_loop_rg_fingerprints/results.json</span>
      </div>
      <div class="view-content">
        ${renderLadderView()}
      </div>
    </div>
  `;
}
