/**
 * TFPT Cascade Viewer - Main Entry Point
 * ======================================
 * 
 * This is the ONLY file that performs I/O.
 * It loads all data into TFPT_STORE, then serves the visualization.
 * 
 * Usage:
 *   npx ts-node index.ts
 * or after compilation:
 *   node dist/index.js
 * 
 * Then open http://localhost:3333 in your browser.
 */

import * as http from 'http';
import * as fs from 'fs';
import * as path from 'path';
import * as yaml from 'js-yaml';

// Loaders
import { loadCascade } from './loader/cascade';
import { loadAlpha } from './loader/alpha';
import { loadSpectrum } from './loader/spectrum';
import { loadTorsion } from './loader/torsion';
import { loadGlobal } from './loader/global';

// State
import { TFPT_STORE, serializeStore, getStoreSummary, setActiveLevel, ACTIVE_LEVEL } from './state/tfpt_store';

// Views
import { getLadderViewHTML } from './views/LadderView';
import { getAlphaMonitorHTML } from './views/AlphaMonitor';
import { getSpectrumHeatmapHTML } from './views/SpectrumHeatmap';
import { getTorsionBoundsHTML } from './views/TorsionBounds';
import { getGlobalRadarHTML } from './views/GlobalRadar';

const PORT = 3333;
const CONFIG_PATH = path.join(__dirname, 'config', 'sources.yaml');

/** Load configuration from YAML */
function loadConfig(): any {
  const configContent = fs.readFileSync(CONFIG_PATH, 'utf8');
  return yaml.load(configContent);
}

/** Resolve path relative to tfpt-suite root */
function resolvePath(basePath: string, relativePath: string): string {
  // basePath is "../../.." from config dir
  // config dir is cascade_viewer/config, so ../../.. goes to tfpt-suite
  const configDir = path.join(__dirname, 'config');
  const tfptSuiteRoot = path.resolve(configDir, basePath);
  return path.join(tfptSuiteRoot, relativePath);
}

/** Initialize the TFPT store by loading all data */
function initializeStore(): void {
  console.log('Loading TFPT data...');
  
  const config = loadConfig();
  const basePath = config.base_path;
  
  // Record source paths
  TFPT_STORE.source_paths = {
    cascade_thresholds: config.cascade.thresholds,
    cascade_two_loop: config.cascade.two_loop_rg,
    alpha_precision: config.alpha.precision_audit,
    spectrum_aps: config.spectrum.aps_eta_gluing,
    spectrum_r2: config.spectrum.effective_action_r2,
    spectrum_bounce: config.spectrum.bounce_perturbations,
    torsion_bounds: config.torsion.bounds_mapping,
    torsion_regimes: config.torsion.regimes,
    global_consistency: config.global.consistency_test
  };
  
  // Load cascade
  console.log('  - Loading cascade data...');
  const cascadeThresholdsPath = resolvePath(basePath, config.cascade.thresholds);
  const cascadeTwoLoopPath = resolvePath(basePath, config.cascade.two_loop_rg);
  TFPT_STORE.cascade = loadCascade(cascadeThresholdsPath, cascadeTwoLoopPath);
  console.log(`    Loaded ${TFPT_STORE.cascade.stages.length} stages, ${TFPT_STORE.cascade.crossings.length} crossings`);
  
  // Load alpha
  console.log('  - Loading alpha data...');
  const alphaPrecisionPath = resolvePath(basePath, config.alpha.precision_audit);
  const alphaDefectG5Path = config.alpha.defect_partition_g5 
    ? resolvePath(basePath, config.alpha.defect_partition_g5) 
    : undefined;
  TFPT_STORE.alpha = loadAlpha(alphaPrecisionPath, alphaDefectG5Path);
  console.log(`    Loaded ${TFPT_STORE.alpha.k_sensitivity.length} k-sensitivity points`);
  if (TFPT_STORE.alpha.defect_g5) {
    console.log(`    Loaded g=${TFPT_STORE.alpha.defect_g5.g_value} defect partition (z=${TFPT_STORE.alpha.defect_g5.z_score.toFixed(2)}σ)`);
  }
  
  // Load spectrum
  console.log('  - Loading spectrum data...');
  const spectrumApsPath = resolvePath(basePath, config.spectrum.aps_eta_gluing);
  const spectrumR2Path = resolvePath(basePath, config.spectrum.effective_action_r2);
  const spectrumBouncePath = resolvePath(basePath, config.spectrum.bounce_perturbations);
  TFPT_STORE.spectrum = loadSpectrum(spectrumApsPath, spectrumR2Path, spectrumBouncePath);
  console.log(`    Loaded ${TFPT_STORE.spectrum.spectral_flow.length} spectral flow entries, ${TFPT_STORE.spectrum.operator_blocks.length} operator blocks`);
  
  // Load torsion
  console.log('  - Loading torsion data...');
  const torsionBoundsPath = resolvePath(basePath, config.torsion.bounds_mapping);
  const torsionRegimesPath = resolvePath(basePath, config.torsion.regimes);
  TFPT_STORE.torsion = loadTorsion(torsionBoundsPath, torsionRegimesPath);
  console.log(`    Loaded ${TFPT_STORE.torsion.inferred_bounds.length} bounds, ${TFPT_STORE.torsion.regimes.length} regimes`);
  
  // Load global
  console.log('  - Loading global consistency data...');
  const globalConsistencyPath = resolvePath(basePath, config.global.consistency_test);
  TFPT_STORE.global = loadGlobal(globalConsistencyPath);
  console.log(`    Loaded ${TFPT_STORE.global.terms.length} χ² terms`);
  
  TFPT_STORE.loaded_at = new Date();
  console.log('All data loaded successfully!');
  console.log('Store summary:', JSON.stringify(getStoreSummary(), null, 2));
}

/** Generate the main HTML page */
function generateMainPage(): string {
  return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>TFPT Cascade Viewer</title>
  <style>
    * {
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
      background: #f0f2f5;
      color: #333;
      min-height: 100vh;
    }
    header {
      background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
      color: white;
      padding: 20px 30px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.2);
    }
    header h1 {
      font-size: 24px;
      font-weight: 600;
    }
    header p {
      font-size: 14px;
      opacity: 0.9;
      margin-top: 8px;
      max-width: 800px;
    }
    .container {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(520px, 1fr));
      gap: 24px;
      padding: 24px;
      max-width: 1900px;
      margin: 0 auto;
    }
    .view-panel {
      background: white;
      border-radius: 12px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      overflow: visible;
    }
    .view-header {
      background: #f5f5f5;
      padding: 14px 18px;
      border-bottom: 1px solid #eee;
    }
    .view-header h3 {
      font-size: 16px;
      font-weight: 600;
      color: #333;
      margin-bottom: 6px;
    }
    .view-description {
      font-size: 12px;
      color: #666;
      line-height: 1.5;
      margin-bottom: 8px;
    }
    .view-source {
      font-size: 10px;
      color: #999;
      font-family: monospace;
    }
    .view-content {
      padding: 20px;
      display: flex;
      justify-content: center;
      align-items: flex-start;
      min-height: 420px;
      overflow: visible;
    }
    .view-content svg {
      display: block;
    }
    footer {
      text-align: center;
      padding: 24px;
      color: #666;
      font-size: 12px;
      border-top: 1px solid #ddd;
      margin-top: 20px;
      background: white;
    }
    footer a {
      color: #1a237e;
      text-decoration: none;
    }
    footer a:hover {
      text-decoration: underline;
    }
    .status-bar {
      background: #e3f2fd;
      padding: 12px 24px;
      font-size: 12px;
      color: #1565c0;
      display: flex;
      justify-content: space-between;
      align-items: center;
      border-bottom: 1px solid #bbdefb;
    }
    .status-bar .active-level {
      background: #1565c0;
      color: white;
      padding: 4px 12px;
      border-radius: 12px;
      font-weight: 600;
    }
    @media (max-width: 1100px) {
      .container {
        grid-template-columns: 1fr;
      }
    }
    .cascade-stage { cursor: pointer; }
    .cascade-stage:hover circle { filter: brightness(1.1); }
  </style>
</head>
<body>
  <header>
    <h1>TFPT Cascade Viewer</h1>
    <p>Interactive visualization of the E8→SM symmetry breaking cascade. All data is read directly from tfpt-suite analysis outputs — no values are computed or duplicated here.</p>
  </header>
  
  <div class="status-bar">
    <span>Data loaded: ${TFPT_STORE.loaded_at?.toISOString().replace('T', ' ').substring(0, 19)} UTC</span>
    <span>Selected cascade level: <span class="active-level">n = ${ACTIVE_LEVEL}</span></span>
  </div>
  
  <div class="container">
    ${getLadderViewHTML()}
    ${getAlphaMonitorHTML()}
    ${getSpectrumHeatmapHTML()}
    ${getTorsionBoundsHTML()}
    ${getGlobalRadarHTML()}
  </div>
  
  <footer>
    <p><strong>TFPT Suite Demo</strong> — Pure read-only visualization of analysis artifacts.</p>
    <p style="margin-top: 8px;">
      Paper reference: TFPT v2.5 — Section 3 (α homeostasis), Section 4 (E8 cascade), Section 6 (Torsion bounds), Appendix K-L (Spectral analysis)
    </p>
  </footer>
  
  <script>
    // Global function for view synchronization
    window.setActiveLevel = function(level) {
      console.log('Setting active level to:', level);
      var statusEl = document.querySelector('.active-level');
      if (statusEl) statusEl.textContent = 'n = ' + level + ' (updating...)';
      
      fetch('/api/set-level/' + level, { method: 'POST' })
        .then(function(response) { return response.json(); })
        .then(function(data) {
          console.log('Level set:', data);
          window.location.reload();
        })
        .catch(function(err) {
          console.error('Error:', err);
          alert('Failed to set level: ' + err);
        });
    };
    
    // Event delegation for cascade stage clicks
    document.addEventListener('click', function(e) {
      var target = e.target;
      // Walk up to find cascade-stage group
      while (target && target !== document) {
        if (target.classList && target.classList.contains('cascade-stage')) {
          var level = target.getAttribute('data-level');
          if (level !== null) {
            console.log('Cascade stage clicked, level:', level);
            window.setActiveLevel(parseInt(level, 10));
          }
          return;
        }
        // Also check if clicked on a child element inside cascade-stage
        if (target.parentElement && target.parentElement.classList && 
            target.parentElement.classList.contains('cascade-stage')) {
          var level = target.parentElement.getAttribute('data-level');
          if (level !== null) {
            console.log('Cascade stage child clicked, level:', level);
            window.setActiveLevel(parseInt(level, 10));
          }
          return;
        }
        target = target.parentElement;
      }
    });
    
    console.log('TFPT Cascade Viewer loaded. Click on cascade stages to change active level.');
  </script>
</body>
</html>`;
}

/** Create and start the HTTP server */
function startServer(): void {
  const server = http.createServer((req, res) => {
    const url = req.url ?? '/';
    
    // API: Get store data
    if (url === '/api/store') {
      res.writeHead(200, { 'Content-Type': 'application/json' });
      res.end(serializeStore());
      return;
    }
    
    // API: Set active level
    if (url.startsWith('/api/set-level/') && req.method === 'POST') {
      const level = parseInt(url.split('/').pop() ?? '0', 10);
      setActiveLevel(level);
      res.writeHead(200, { 'Content-Type': 'application/json' });
      res.end(JSON.stringify({ ok: true, level }));
      return;
    }
    
    // API: Get summary
    if (url === '/api/summary') {
      res.writeHead(200, { 'Content-Type': 'application/json' });
      res.end(JSON.stringify(getStoreSummary(), null, 2));
      return;
    }
    
    // Main page
    if (url === '/' || url === '/index.html') {
      res.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8' });
      res.end(generateMainPage());
      return;
    }
    
    // 404
    res.writeHead(404, { 'Content-Type': 'text/plain' });
    res.end('Not Found');
  });
  
  server.listen(PORT, () => {
    console.log('');
    console.log('========================================');
    console.log('  TFPT Cascade Viewer');
    console.log(`  Server running at http://localhost:${PORT}`);
    console.log('========================================');
    console.log('');
    console.log('API endpoints:');
    console.log('  GET  /api/store   - Full store data (JSON)');
    console.log('  GET  /api/summary - Store summary');
    console.log('  POST /api/set-level/:n - Set active cascade level');
    console.log('');
  });
}

// Main
try {
  initializeStore();
  startServer();
} catch (error) {
  console.error('Failed to initialize TFPT Cascade Viewer:', error);
  process.exit(1);
}
