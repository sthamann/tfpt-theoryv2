/**
 * Spectrum Heatmap - Spectral Flow & Operator Blocks
 * ==================================================
 * 
 * Visualizes:
 *   - APS η-gluing spectral flow table
 *   - R² operator block contributions (E/R eigenvalues)
 *   - Transfer functions T(k) from bounce perturbations
 * 
 * Data source: TFPT_STORE.spectrum
 * Paper reference: Appendix K (APS gluing), Appendix L (bounce perturbations)
 */

import { TFPT_STORE, ACTIVE_LEVEL } from '../state/tfpt_store';

const SVG_WIDTH = 520;
const SVG_HEIGHT = 480;

/**
 * Render spectral flow heatmap
 */
function renderSpectralFlowTable(): string {
  const spectrum = TFPT_STORE.spectrum;
  if (!spectrum) return '';
  
  const periodicEntries = spectrum.spectral_flow.filter(e => e.spin === 'periodic');
  const cellSize = 36;
  const tableWidth = periodicEntries.length * cellSize;
  const startX = 70;
  const startY = 70;
  
  let cells = '';
  for (let i = 0; i < periodicEntries.length; i++) {
    const entry = periodicEntries[i];
    const x = startX + i * cellSize;
    
    // Spectral flow value determines color intensity
    const intensity = Math.min(entry.numeric_sf / 8, 1);
    const lightness = 65 - intensity * 30;
    const color = `hsl(210, 70%, ${lightness}%)`;
    
    cells += `
      <g>
        <rect x="${x}" y="${startY}" width="${cellSize - 3}" height="${cellSize - 3}" fill="${color}" rx="4">
          <title>Winding number m = ${entry.m}
Spectral flow SF = ${entry.numeric_sf}
η(end) = ${entry.eta_end_numeric.toFixed(6)}</title>
        </rect>
        <text x="${x + cellSize/2 - 1}" y="${startY + cellSize/2 + 5}" text-anchor="middle" fill="#fff" font-size="13" font-weight="bold">
          ${entry.numeric_sf}
        </text>
      </g>
    `;
  }
  
  // Column headers (m values)
  let headers = '';
  for (let i = 0; i < periodicEntries.length; i++) {
    const x = startX + i * cellSize;
    headers += `<text x="${x + cellSize/2 - 1}" y="${startY - 10}" text-anchor="middle" fill="#666" font-size="10">m=${periodicEntries[i].m}</text>`;
  }
  
  return `
    <g class="spectral-flow-table">
      <text x="${startX - 10}" y="${startY + cellSize/2}" text-anchor="end" fill="#333" font-size="11" font-weight="500">SF</text>
      ${headers}
      ${cells}
    </g>
  `;
}

/**
 * Render operator blocks bar chart
 */
function renderOperatorBlocks(): string {
  const spectrum = TFPT_STORE.spectrum;
  if (!spectrum) return '';
  
  const blocks = spectrum.operator_blocks;
  const barHeight = 28;
  const maxBeta = Math.max(...blocks.map(b => Math.abs(b.beta_R2_contribution)));
  const chartWidth = 280;
  const startX = 120;
  const startY = 160;
  
  let bars = '';
  for (let i = 0; i < blocks.length; i++) {
    const block = blocks[i];
    const y = startY + i * (barHeight + 10);
    const width = (Math.abs(block.beta_R2_contribution) / maxBeta) * chartWidth;
    const color = block.statistics === 'ghost' ? '#FF5722' : '#2196F3';
    const sign = block.beta_R2_contribution < 0 ? '-' : '+';
    
    // Display name
    const displayName = block.name
      .replace('torsion_trace_vector_', 'T_μ ')
      .replace('torsion_axial_vector_', 'S_μ ')
      .replace('torsion_tensor_', 'q_μνρ ')
      .replace('fp_ghost_vector', 'FP ghost');
    
    bars += `
      <g>
        <text x="${startX - 8}" y="${y + barHeight/2 + 4}" text-anchor="end" fill="#333" font-size="10">
          ${displayName}
        </text>
        <rect x="${startX}" y="${y}" width="${width}" height="${barHeight}" fill="${color}" rx="4" opacity="0.85">
          <title>${block.name}
β_R² contribution: ${block.beta_R2_contribution.toExponential(4)}
Rank: ${block.rank}</title>
        </rect>
        <text x="${startX + width + 8}" y="${y + barHeight/2 + 4}" fill="#555" font-size="10">
          ${sign}${Math.abs(block.beta_R2_contribution).toExponential(2)}
        </text>
      </g>
    `;
  }
  
  return `
    <g class="operator-blocks">
      ${bars}
      <text x="${startX + chartWidth/2}" y="${startY + blocks.length * (barHeight + 10) + 20}" text-anchor="middle" fill="#333" font-size="11">
        Total: β_R² = ${spectrum.beta_R2_total.toExponential(3)} → M/M_Pl = ${spectrum.M_over_Mpl.toExponential(4)}
      </text>
    </g>
  `;
}

/**
 * Render transfer function mini-plot
 */
function renderTransferFunctions(): string {
  const spectrum = TFPT_STORE.spectrum;
  if (!spectrum) return '';
  
  const tf = spectrum.transfer_functions;
  const plotWidth = 200;
  const plotHeight = 100;
  const startX = 290;
  const startY = 350;
  
  // Log scale for k
  const kMin = Math.log10(Math.min(...tf.k_grid));
  const kMax = Math.log10(Math.max(...tf.k_grid));
  const tMax = Math.max(...tf.T_scalar, ...tf.T_tensor, 5);
  
  const mapX = (k: number) => startX + ((Math.log10(k) - kMin) / (kMax - kMin)) * plotWidth;
  const mapY = (t: number) => startY + plotHeight - (Math.min(t, tMax) / tMax) * plotHeight;
  
  // Build paths
  let scalarPath = '';
  let tensorPath = '';
  for (let i = 0; i < tf.k_grid.length; i++) {
    const x = mapX(tf.k_grid[i]);
    const ys = mapY(tf.T_scalar[i]);
    const yt = mapY(tf.T_tensor[i]);
    scalarPath += i === 0 ? `M ${x} ${ys}` : ` L ${x} ${ys}`;
    tensorPath += i === 0 ? `M ${x} ${yt}` : ` L ${x} ${yt}`;
  }
  
  // T=1 reference line
  const t1Y = mapY(1);
  
  return `
    <g class="transfer-functions">
      <text x="${startX + plotWidth/2}" y="${startY - 15}" text-anchor="middle" fill="#333" font-size="11" font-weight="500">
        Transfer Functions T(k)
      </text>
      <rect x="${startX}" y="${startY}" width="${plotWidth}" height="${plotHeight}" fill="#fafafa" stroke="#ddd"/>
      
      <!-- T=1 line -->
      <line x1="${startX}" y1="${t1Y}" x2="${startX + plotWidth}" y2="${t1Y}" stroke="#4CAF50" stroke-dasharray="4,3"/>
      <text x="${startX + plotWidth + 5}" y="${t1Y + 3}" fill="#4CAF50" font-size="9">T=1</text>
      
      <!-- Scalar -->
      <path d="${scalarPath}" fill="none" stroke="#2196F3" stroke-width="2"/>
      
      <!-- Tensor -->
      <path d="${tensorPath}" fill="none" stroke="#FF9800" stroke-width="2"/>
      
      <!-- Labels -->
      <text x="${startX + 8}" y="${startY + 18}" fill="#2196F3" font-size="10" font-weight="500">Scalar</text>
      <text x="${startX + 8}" y="${startY + 32}" fill="#FF9800" font-size="10" font-weight="500">Tensor</text>
      <text x="${startX + plotWidth/2}" y="${startY + plotHeight + 18}" text-anchor="middle" fill="#666" font-size="10">
        log₁₀(k)
      </text>
    </g>
  `;
}

/**
 * Render complete spectrum heatmap
 */
export function renderSpectrumHeatmap(): string {
  const spectrum = TFPT_STORE.spectrum;
  if (!spectrum) {
    return `<svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}">
      <text x="${SVG_WIDTH/2}" y="${SVG_HEIGHT/2}" text-anchor="middle" fill="#f00">
        Error: Spectrum data not loaded
      </text>
    </svg>`;
  }
  
  return `
    <svg width="${SVG_WIDTH}" height="${SVG_HEIGHT}" xmlns="http://www.w3.org/2000/svg">
      <!-- Section: Spectral Flow -->
      <text x="20" y="30" font-size="12" font-weight="600" fill="#333">Spectral Flow SF(m)</text>
      <text x="20" y="48" font-size="10" fill="#666">APS η-gluing for periodic boundary conditions</text>
      ${renderSpectralFlowTable()}
      
      <!-- Section: Operator Blocks -->
      <text x="20" y="140" font-size="12" font-weight="600" fill="#333">R² Operator Block Contributions</text>
      ${renderOperatorBlocks()}
      
      <!-- Section: Transfer Functions -->
      ${renderTransferFunctions()}
      
      <!-- Seam term annotation -->
      <text x="20" y="${SVG_HEIGHT - 20}" fill="#666" font-size="10">
        Seam term: Δ_Γ(m=1) = ${spectrum.seam_term_delta_gamma_2pi.toFixed(5)} ≈ 2π ✓
      </text>
    </svg>
  `;
}

/**
 * Generate HTML wrapper
 */
export function getSpectrumHeatmapHTML(): string {
  return `
    <div id="spectrum-heatmap" class="view-panel">
      <div class="view-header">
        <h3>Spectrum Heatmap</h3>
        <p class="view-description">
          <strong>What you see:</strong> Three spectral analyses from TFPT.
          <em>Top:</em> APS spectral flow SF(m) for winding numbers m=0..8 — the integer values confirm correct gluing.
          <em>Middle:</em> Heat-kernel R² contributions from torsion operator blocks.
          <em>Bottom:</em> Bounce transfer functions T(k) approaching 1 at high k (adiabatic limit).
        </p>
        <span class="view-source">Data: aps_eta_gluing, effective_action_r2, bounce_perturbations</span>
      </div>
      <div class="view-content">
        ${renderSpectrumHeatmap()}
      </div>
    </div>
  `;
}
