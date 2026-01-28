"""
Chain search module for finding optimal monotonic chains in E8 nilpotent orbit lattice.

This module implements a beam search algorithm to explore the combinatorial space
of valid Bala-Carter chains, evaluating them with fit-free theoretical criteria.
"""

import numpy as np
import pandas as pd
import re
from collections import Counter
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


def tokenize_label(label: str) -> List[str]:
    """
    Tokenize orbit label into constituent parts.
    
    Examples:
    - "E6(a3)+A1" -> ["E6(a3)", "A1"]
    - "2A3" -> ["A3", "A3"]
    - "A4+2A1" -> ["A4", "A1", "A1"]
    """
    parts = []
    # Split by + and process each part
    for part in label.replace(' ', '').split('+'):
        # Match patterns like 2A3, A4, D5(a1), E6(a3), etc.
        m = re.match(r'(\d+)?([ADE]\d+(?:\([a-z]\d+\))?|E8(?:\([a-z]\d+\))?)', part.strip())
        if m:
            count = int(m.group(1) or 1)
            base = m.group(2)
            parts.extend([base] * count)
    return parts


def jaccard_distance_multiset(tokens_a: List[str], tokens_b: List[str]) -> float:
    """
    Calculate Jaccard distance between two multisets of tokens.
    
    Returns value in [0, 1] where 0 means identical, 1 means completely different.
    """
    if not tokens_a and not tokens_b:
        return 0.0
    
    ca = Counter(tokens_a)
    cb = Counter(tokens_b)
    
    # Intersection and union of multisets
    intersection = sum((ca & cb).values())
    union = sum((ca | cb).values())
    
    if union == 0:
        return 0.0
    
    return 1.0 - (intersection / union)


def build_hasse_graph(df: pd.DataFrame, 
                      label_threshold: float = 1.0,  # Default: no filter
                      max_height_diff: int = 999,
                      allowed_steps: Tuple[int, ...] = (2, 4)) -> Tuple[Dict, Dict]:
    """
    Build approximate Hasse diagram from orbit data.
    
    Parameters
    ----------
    df : pd.DataFrame
        Orbit data with columns: Label, Dim, Height
    label_threshold : float
        Maximum allowed label distance for edge creation
    max_height_diff : int
        Maximum allowed height difference for edges
        
    Returns
    -------
    layers : dict
        Orbits grouped by D value
    edges : dict
        Adjacency list of edges with metadata
    """
    # Prepare data
    df = df.copy()
    df['D'] = 248 - df['Dim']
    df['tokens'] = df['Label'].apply(tokenize_label)
    
    # Group by D value (layers in the Hasse diagram)
    layers = {}
    for D, group in df.groupby('D'):
        layers[int(D)] = group
    
    # Build edges (directed, downward in D)
    edges = {i: [] for i in df.index}
    total_edges = 0
    
    for D in sorted(layers.keys(), reverse=True):
        layer_edges = 0
        # Look for connections to next layers based on allowed steps
        for step in allowed_steps:
            target_D = D - step
            if target_D not in layers:
                continue
            
            for i, u in layers[D].iterrows():
                for j, v in layers[target_D].iterrows():
                    # Calculate metrics (no hard filtering, only features)
                    height_diff = abs(v['Height'] - u['Height'])
                    label_dist = jaccard_distance_multiset(u['tokens'], v['tokens'])
                    
                    # Only apply VERY loose filters if specified
                    # These should primarily be in the cost function, not as gates
                    if max_height_diff < 999 and height_diff > max_height_diff:
                        continue
                    if label_threshold < 1.0 and label_dist > label_threshold:
                        continue
                    
                    # Add edge with ALL metadata for cost evaluation
                    edges[i].append((j, {
                        'delta_D': step,
                        'delta_height': int(v['Height'] - u['Height']),
                        'label_dist': label_dist,
                        'height_diff': height_diff
                    }))
                    layer_edges += 1
                    total_edges += 1
        
        # Per-layer diagnostics
        if layer_edges > 0:
            logger.debug(f"Layer D={D}: {len(layers[D])} nodes, {layer_edges} outgoing edges")
    
    # Summary
    logger.info(f"Built Hasse graph: {len(layers)} layers, {total_edges} edges")
    
    # Sanity check
    if total_edges < len(layers) - 1:
        logger.warning(f"Graph is too sparse! Only {total_edges} edges for {len(layers)} layers.")
        logger.warning("Consider relaxing label_threshold or max_height_diff constraints.")
    
    return layers, edges


def calculate_delta3(x: np.ndarray) -> np.ndarray:
    """
    Calculate third forward difference: Δ³x_n.
    
    Δ³x_n = x_{n+3} - 3*x_{n+2} + 3*x_{n+1} - x_n
    """
    if len(x) < 4:
        return np.array([])
    return x[3:] - 3*x[2:-1] + 3*x[1:-2] - x[:-3]


def evaluate_chain_metrics(lnD: List[float], 
                           heights: List[int], 
                           edge_meta: List[Dict]) -> Dict[str, float]:
    """
    Calculate fit-free metrics for chain evaluation.
    
    Returns
    -------
    dict
        Metrics including:
        - cv_d3: Coefficient of variation of Δ³ln(D)
        - smooth_s: RMS of second differences of step sizes
        - sum_delta_height: Sum of absolute height changes
        - label_penalty: Cumulative label distance
        - n_jumps: Number of D jumps of size 4
    """
    lnD = np.array(lnD)
    heights = np.array(heights)
    
    # Metric 1: CV of third differences of ln(D)
    d3 = calculate_delta3(lnD)
    if len(d3) > 0 and np.mean(d3) != 0:
        cv_d3 = np.std(d3) / abs(np.mean(d3))
    else:
        cv_d3 = np.inf
    
    # Metric 2: Smoothness of step sizes
    if len(lnD) >= 2:
        s = np.diff(lnD)  # Step sizes
        if len(s) >= 2:
            d2_s = np.diff(s)  # Second differences
            smooth_s = np.sqrt(np.mean(d2_s**2))
        else:
            smooth_s = 0.0
    else:
        smooth_s = np.inf
    
    # Metric 3: Sum of absolute height changes
    if len(heights) >= 2:
        sum_delta_height = np.sum(np.abs(np.diff(heights)))
    else:
        sum_delta_height = 0.0
    
    # Metric 4: Cumulative label distance
    label_penalty = sum(e.get('label_dist', 0) for e in edge_meta)
    
    # Metric 5: Number of large jumps (ΔD = 4)
    n_jumps = sum(1 for e in edge_meta if e.get('delta_D', 0) == 4)
    
    return {
        'cv_d3': float(cv_d3),
        'smooth_s': float(smooth_s),
        'sum_delta_height': float(sum_delta_height),
        'label_penalty': float(label_penalty),
        'n_jumps': int(n_jumps)
    }


def combine_cost(metrics: Dict[str, float], weights: Dict[str, float]) -> float:
    """
    Combine metrics into single cost value using weights.
    """
    return (weights['w_cubic'] * metrics['cv_d3'] +
            weights['w_smooth'] * metrics['smooth_s'] +
            weights['w_height'] * metrics['sum_delta_height'] +
            weights['w_label'] * metrics['label_penalty'] +
            weights['w_jump'] * metrics['n_jumps'])


def beam_search_chains(df: pd.DataFrame,
                       beam_width: int = 32,
                       top_k: int = 10,
                       weights: Optional[Dict[str, float]] = None,
                       label_threshold: float = 1.0,  # Default: no label filtering
                       max_height_diff: int = 999,     # Default: no height filtering
                       allowed_steps: Tuple[int, ...] = (2, 4),  # Allowed ΔD steps
                       start_D: int = 60,
                       target_D: int = 8,
                       verbose: bool = True) -> List[Dict]:
    """
    Find optimal chains using beam search.
    
    Parameters
    ----------
    df : pd.DataFrame
        Orbit data
    beam_width : int
        Maximum number of paths to keep at each step
    top_k : int
        Number of best chains to return
    weights : dict
        Weights for cost function
    label_threshold : float
        Maximum label distance for edges
    start_D : int
        Starting D value (should have A4+A1)
    target_D : int
        Target D value
        
    Returns
    -------
    list
        Top k chains with their metrics and paths
    """
    # Default weights
    if weights is None:
        weights = {
            'w_cubic': 1.0,    # Most important: cubic pattern in ln(D)
            'w_smooth': 0.6,   # Smoothness of steps
            'w_height': 0.3,   # Height continuity
            'w_label': 0.2,    # Label similarity
            'w_jump': 0.5      # Penalty for large jumps
        }
    
    # Build graph
    if verbose:
        logger.info(f"Building Hasse graph with label_threshold={label_threshold}, max_height_diff={max_height_diff}, allowed_steps={allowed_steps}")
    
    layers, edges = build_hasse_graph(df, label_threshold, max_height_diff, allowed_steps)
    
    # Check if graph is viable
    total_edges = sum(len(e) for e in edges.values())
    expected_min_edges = len(layers) - 1  # Minimum for a connected path
    
    if total_edges < expected_min_edges:
        logger.warning(f"Sparse graph: {total_edges} edges across {len(layers)} layers (expected at least {expected_min_edges}); continuing anyway.")
    
    if verbose and allowed_steps == (2,):
        logger.info(f"Strict ΔD=2 mode: expecting {(start_D - target_D)//2 + 1} steps")
    
    # Find starting points (must include A4+A1 at D=60)
    if start_D not in layers:
        logger.error(f"No orbits found at D={start_D}")
        return []
    
    start_layer = layers[start_D]
    start_indices = []
    for idx, row in start_layer.iterrows():
        if row['Label'].replace(' ', '') == 'A4+A1':
            start_indices.append(idx)
            logger.info(f"Starting from {row['Label']} (D={start_D})")
    
    if not start_indices:
        logger.warning("A4+A1 not found at D=60, using all orbits at this level")
        start_indices = list(start_layer.index)
    
    # Initialize beam with starting states
    # State: (path_indices, D_values, lnD_values, heights, edge_metadata)
    beam = []
    for idx in start_indices:
        row = df.loc[idx]
        beam.append({
            'path': [idx],
            'D_values': [int(row['D'])],
            'lnD_values': [np.log(row['D'])],
            'heights': [int(row['Height'])],
            'edge_meta': [],
            'cost': 0.0
        })
    
    # Beam search with progress tracking
    current_D = start_D
    iteration = 0
    
    if verbose:
        logger.info(f"Starting beam search from D={start_D} to D={target_D}")
    
    while current_D > target_D and beam:
        next_beam = []
        
        for state in beam:
            last_idx = state['path'][-1]
            
            # Explore all edges from current node
            for next_idx, edge_info in edges.get(last_idx, []):
                next_row = df.loc[next_idx]
                next_D = int(next_row['D'])
                
                # Skip if not going down
                if next_D >= current_D:
                    continue
                
                # Create new state
                new_state = {
                    'path': state['path'] + [next_idx],
                    'D_values': state['D_values'] + [next_D],
                    'lnD_values': state['lnD_values'] + [np.log(next_D)],
                    'heights': state['heights'] + [int(next_row['Height'])],
                    'edge_meta': state['edge_meta'] + [edge_info],
                    'cost': 0.0
                }
                
                # Calculate cost if we have enough points
                if len(new_state['lnD_values']) >= 4:
                    metrics = evaluate_chain_metrics(
                        new_state['lnD_values'],
                        new_state['heights'],
                        new_state['edge_meta']
                    )
                    new_state['cost'] = combine_cost(metrics, weights)
                    new_state['metrics'] = metrics
                
                next_beam.append(new_state)
        
        # Prune beam
        if next_beam:
            next_beam.sort(key=lambda s: s['cost'])
            beam = next_beam[:beam_width]
            current_D = min(s['D_values'][-1] for s in beam)
            iteration += 1
            
            if verbose and iteration % 5 == 0:
                best_cost = beam[0]['cost'] if beam else float('inf')
                logger.info(f"  Iteration {iteration}: D={current_D}, beam_size={len(beam)}, best_cost={best_cost:.4f}")
            else:
                logger.debug(f"Beam at D={current_D}: {len(beam)} paths")
    
    # Final scoring and ranking
    complete_chains = []
    for state in beam:
        if state['D_values'][-1] <= target_D:
            # Recalculate final metrics
            metrics = evaluate_chain_metrics(
                state['lnD_values'],
                state['heights'],
                state['edge_meta']
            )
            state['cost'] = combine_cost(metrics, weights)
            state['metrics'] = metrics
            
            # Add labels
            state['labels'] = [df.loc[idx, 'Label'] for idx in state['path']]
            state['dims'] = [int(df.loc[idx, 'Dim']) for idx in state['path']]
            
            complete_chains.append(state)
    
    # Sort by cost and return top k
    complete_chains.sort(key=lambda s: s['cost'])
    
    if complete_chains:
        logger.info(f"Found {len(complete_chains)} complete chains, returning top {min(top_k, len(complete_chains))}")
        logger.info(f"Best chain cost: {complete_chains[0]['cost']:.4f}, CV(Δ³ln(D)): {complete_chains[0]['metrics']['cv_d3']:.4f}")
    else:
        logger.warning("No complete chains found")
    
    return complete_chains[:top_k]


def save_chain_results(chains: List[Dict], output_dir: Path) -> None:
    """
    Save chain search results to files.
    
    Creates:
    - top10/rank_N.csv for each chain
    - ranking.json with metrics for all chains
    """
    import json
    
    output_dir = Path(output_dir)
    chains_dir = output_dir / 'chains' / 'top10'
    chains_dir.mkdir(parents=True, exist_ok=True)
    
    ranking = []
    
    for rank, chain in enumerate(chains, 1):
        # Save chain data
        chain_df = pd.DataFrame({
            'n': range(len(chain['labels'])),
            'label': chain['labels'],
            'dim': chain['dims'],
            'D': chain['D_values'],
            'lnD': chain['lnD_values'],
            'height': chain['heights']
        })
        
        # Add step sizes
        chain_df['s_n'] = 0.0
        if len(chain['lnD_values']) > 1:
            s = np.diff(chain['lnD_values'])
            chain_df.loc[1:, 's_n'] = s
        
        # Save to CSV
        csv_path = chains_dir / f'rank_{rank}.csv'
        chain_df.to_csv(csv_path, index=False)
        
        # Add to ranking
        ranking.append({
            'rank': rank,
            'cost': chain['cost'],
            'metrics': chain['metrics'],
            'n_steps': len(chain['labels']),
            'file': f'rank_{rank}.csv'
        })
        
        logger.info(f"Saved rank {rank} chain to {csv_path}")
    
    # Save ranking
    ranking_path = output_dir / 'chains' / 'ranking.json'
    with open(ranking_path, 'w') as f:
        json.dump(ranking, f, indent=2)
    
    logger.info(f"Saved ranking to {ranking_path}")


def explain_chain(chain: Dict, df: pd.DataFrame) -> str:
    """
    Generate explanation for why this chain was selected.
    
    Returns formatted explanation string.
    """
    lines = []
    lines.append("CHAIN SELECTION EXPLANATION")
    lines.append("=" * 50)
    
    # Overall metrics
    metrics = chain['metrics']
    lines.append("\nOVERALL METRICS:")
    lines.append(f"  • CV(Δ³ln(D)): {metrics['cv_d3']:.4f} (lower is better)")
    lines.append(f"  • Smoothness RMS: {metrics['smooth_s']:.6f}")
    lines.append(f"  • Total height change: {metrics['sum_delta_height']:.0f}")
    lines.append(f"  • Label distance sum: {metrics['label_penalty']:.4f}")
    lines.append(f"  • Large jumps (ΔD=4): {metrics['n_jumps']}")
    lines.append(f"  • Total cost: {chain['cost']:.4f}")
    
    # Step-by-step breakdown
    lines.append("\nSTEP-BY-STEP BREAKDOWN:")
    for i in range(len(chain['labels']) - 1):
        lines.append(f"\nStep {i}: {chain['labels'][i]} → {chain['labels'][i+1]}")
        lines.append(f"  D: {chain['D_values'][i]} → {chain['D_values'][i+1]} (ΔD = {chain['D_values'][i] - chain['D_values'][i+1]})")
        lines.append(f"  Height: {chain['heights'][i]} → {chain['heights'][i+1]} (Δh = {chain['heights'][i+1] - chain['heights'][i]})")
        
        if i < len(chain['edge_meta']):
            meta = chain['edge_meta'][i]
            lines.append(f"  Label distance: {meta.get('label_dist', 0):.3f}")
    
    # Cubic test details
    lnD = np.array(chain['lnD_values'])
    d3 = calculate_delta3(lnD)
    if len(d3) > 0:
        lines.append("\nCUBIC TEST DETAILS:")
        lines.append(f"  Δ³ln(D) values: {d3}")
        lines.append(f"  Mean: {np.mean(d3):.6f}")
        lines.append(f"  Std: {np.std(d3):.6f}")
        lines.append(f"  CV: {metrics['cv_d3']:.4f}")
    
    return "\n".join(lines)
