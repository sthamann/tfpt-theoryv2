#!/usr/bin/env python3
"""
Command-line interface for E8 Orbit Engine.
"""

import typer
from pathlib import Path
from typing import Optional
from rich.console import Console
from rich.table import Table
from rich.progress import track
import numpy as np
import json
import logging

app = typer.Typer(
    name="e8-orbit",
    help="E8 Orbit Engine - Derive quadratic damping from nilpotent orbits",
    add_completion=False
)
console = Console()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


@app.command()
def solve(
    rebuild_e8: bool = typer.Option(
        False,
        "--rebuild-e8",
        help="Rebuild E8 chain from orbit data"
    ),
    orbits: Optional[Path] = typer.Option(
        None,
        "--orbits",
        help="Path to nilpotent orbits CSV file"
    ),
    out: Optional[Path] = typer.Option(
        None,
        "--out",
        help="Output path for E8 chain CSV"
    ),
    verify_rge: bool = typer.Option(
        False,
        "--verify-rge",
        help="Verify RGE crossings for E6/E7/E8"
    ),
    rge: Optional[Path] = typer.Option(
        None,
        "--rge",
        help="Path to RGE gauge coupling data"
    ),
    all: bool = typer.Option(
        False,
        "--all",
        help="Run all tasks (rebuild E8 and verify RGE)"
    ),
    output_dir: Path = typer.Option(
        Path("results"),
        "--output-dir",
        help="Directory for all output files"
    ),
    search_chains: bool = typer.Option(
        True,
        "--search-chains/--no-search-chains",
        help="Search for optimal chain using beam search (default: True)"
    ),
    baseline: bool = typer.Option(
        False,
        "--baseline",
        help="Use baseline predefined chain instead of searching"
    ),
    beam_width: int = typer.Option(
        32,
        "--beam-width",
        help="Beam width for chain search"
    ),
    top_k: int = typer.Option(
        10,
        "--top-k",
        help="Number of top chains to save"
    ),
    label_threshold: float = typer.Option(
        1.0,
        "--label-threshold",
        help="Maximum label distance for edge creation (1.0 = no filter)"
    ),
    max_height_diff: int = typer.Option(
        999,
        "--max-height-diff",
        help="Maximum height difference for edges (999 = no filter)"
    ),
    allowed_steps: str = typer.Option(
        "2,4",
        "--allowed-steps",
        help="Comma-separated allowed ΔD steps (e.g. '2' for strict ΔD=2, or '2,4' for both)"
    ),
    html: bool = typer.Option(
        True,
        "--html/--no-html",
        help="Generate HTML report (default: on)"
    )
):
    """
    Main solve command with multiple modes:
    
    - --rebuild-e8: Build E8 chain and compute gamma
    - --verify-rge: Find E6/E7/E8 crossing scales
    - --all: Run both tasks
    """
    from . import (
        load_orbits, load_rge, build_chain, 
        compute_gamma, quad_and_d3,
        verify_rge_crossings, generate_report,
        compare_models, fit_gamma_log_model, fit_gamma_hyperbolic
    )
    from .chain import validate_chain
    from .fit import assert_first_step, fit_gamma, analyze_fit_quality
    from .verify import comprehensive_verification, phi_third_diff_from_data, verify_lnD_cubic_rolling
    from .io import save_chain as save_chain_func
    
    console.print("\n[bold cyan]E8 Orbit Engine[/bold cyan]")
    console.print("=" * 50)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine what to run
    run_e8 = rebuild_e8 or all
    run_rge = verify_rge or all
    
    chain_result = None
    fit_result = None
    
    # Task 1: Rebuild E8 chain
    if run_e8:
        console.print("\n[yellow]Task 1:[/yellow] Building E8 chain...")
        
        # Set default paths
        if orbits is None:
            orbits = Path("data/nilpotent_orbits.csv")
        if out is None:
            out = output_dir / "data" / "e8_chain_clean.csv"
        
        try:
            # Load orbit data
            df_orbits = load_orbits(orbits)
            console.print(f"  ✓ Loaded {len(df_orbits)} orbits")
            
            # Check for A4+A1 correction
            if "A4+A1" in df_orbits["Label"].values:
                idx = df_orbits[df_orbits["Label"] == "A4+A1"].index[0]
                if df_orbits.loc[idx, "Dim"] == 188:
                    console.print("  [yellow]ℹ A4+A1 dimension was corrected: 206 → 188 (source: E8 table)[/yellow]")
            
            # Build chain
            if baseline:
                console.print("  Using baseline chain...")
                chain = build_chain(df_orbits, use_baseline=True)
            elif search_chains:
                console.print(f"  Searching for optimal chain (beam_width={beam_width})...")
                from .chain_search import beam_search_chains, save_chain_results
                import pandas as pd
                
                # Parse allowed steps
                steps = tuple(int(x.strip()) for x in allowed_steps.split(",") if x.strip())
                
                # Search for chains with relaxed filters
                console.print(f"    Label threshold: {label_threshold} (1.0 = no filter)")
                console.print(f"    Max height diff: {max_height_diff} (999 = no filter)")
                console.print(f"    Allowed ΔD steps: {steps}")
                if steps == (2,):
                    console.print(f"    [green]Strict ΔD=2 mode: expecting 27-step chain[/green]")
                
                chains = beam_search_chains(
                    df_orbits, 
                    beam_width=beam_width, 
                    top_k=top_k,
                    label_threshold=label_threshold,
                    max_height_diff=max_height_diff,
                    allowed_steps=steps,
                    verbose=True
                )
                
                if chains:
                    # Save all top chains
                    save_chain_results(chains, output_dir)
                    
                    # Use best chain
                    best = chains[0]
                    chain = pd.DataFrame({
                        'n': range(len(best['labels'])),
                        'label': best['labels'],
                        'dim_orbit': best['dims'],
                        'D': best['D_values'],
                        'lnD': best['lnD_values']
                    })
                    
                    console.print(f"\n  [green]✓ Found {len(chains)} valid chains![/green]")
                    console.print(f"\n  [bold]Best chain (rank 1):[/bold]")
                    console.print(f"    Total cost: {best['cost']:.4f}")
                    console.print(f"    CV(Δ³ln(D)): {best['metrics']['cv_d3']:.4f} (lower is better)")
                    console.print(f"    Smoothness: {best['metrics']['smooth_s']:.6f}")
                    console.print(f"    Height changes: {best['metrics']['sum_delta_height']:.0f}")
                    console.print(f"    Label penalty: {best['metrics']['label_penalty']:.3f}")
                    console.print(f"    Large jumps (ΔD=4): {best['metrics']['n_jumps']}")
                    
                    # Show brief summary of other chains
                    if len(chains) > 1:
                        console.print(f"\n  Other chains found (saved to results/chains/top10/):")
                        for i, ch in enumerate(chains[1:min(5, len(chains))], 2):
                            console.print(f"    Rank {i}: cost={ch['cost']:.4f}, CV={ch['metrics']['cv_d3']:.4f}")
                else:
                    console.print("  [red]✗ No valid chains found![/red]")
                    console.print("  [yellow]Falling back to baseline chain...[/yellow]")
                    chain = build_chain(df_orbits, use_baseline=True)
            else:
                # Default: use build_chain which will search
                chain = build_chain(df_orbits)
            
            console.print(f"  ✓ Built chain with {len(chain)} orbits")
            console.print(f"  D range: {chain['D'].max()} → {chain['D'].min()}")
            
            # Validate chain
            valid, messages = validate_chain(chain)
            for msg in messages:
                if "✓" in msg:
                    console.print(f"    [green]{msg}[/green]")
                elif "✗" in msg:
                    console.print(f"    [red]{msg}[/red]")
                else:
                    console.print(f"    {msg}")
            
            if not valid:
                if not typer.confirm("Chain validation failed. Continue anyway?"):
                    raise typer.Exit(1)
            
            # Check first step normalization
            try:
                assert_first_step(chain)
                console.print("  ✓ First step normalization correct")
            except AssertionError as e:
                console.print(f"  [yellow]⚠ {e}[/yellow]")
            
            # Compute gamma (without fitting)
            console.print("\n[yellow]Computing γ(n) from data...[/yellow]")
            s, gamma, lam, s_star = compute_gamma(chain)
            console.print(f"  Normalization: s* = ln(248) - ln(60) = {s_star:.6f}")
            console.print(f"  Scaling: λ = 0.834 / s* = {lam:.6f}")
            console.print(f"  γ(0) = {gamma[0]:.6f}")
            
            # Quadratic diagnostics
            console.print("\n[yellow]Quadratic diagnostics (n≥0)...[/yellow]")
            diag = quad_and_d3(chain)
            console.print(f"  Quadratic fit R² = {diag['r2']:.6f} [red](poor fit)[/red]")
            console.print(f"  Third differences: mean = {diag['d3_mean']:.8f}, var = {diag['d3_var']:.8e}")
            console.print(f"  [yellow]⚠ Note: n=0 is an outlier (248→60 step)[/yellow]")
            
            # Full fit for comprehensive results
            fit_result = fit_gamma(chain, gamma0_target=0.834)
            quality = analyze_fit_quality(fit_result)
            
            # Model comparison
            console.print("\n[yellow]Model comparison (n≥1)...[/yellow]")
            try:
                comparison = compare_models(chain, fit_result)
                
                # Create a nice table
                from rich.table import Table
                table = Table(title="Model Performance")
                table.add_column("Model", style="cyan")
                table.add_column("Params", justify="right")
                table.add_column("R²", justify="right")
                table.add_column("AIC", justify="right")
                table.add_column("BIC", justify="right")
                table.add_column("RMSE", justify="right")
                
                for _, row in comparison.iterrows():
                    r2_style = "green" if row['R²'] > 0.9 else "yellow" if row['R²'] > 0.5 else "red"
                    table.add_row(
                        row['Model'],
                        str(int(row['Parameters'])),
                        f"[{r2_style}]{row['R²']:.4f}[/{r2_style}]",
                        f"{row['AIC']:.1f}",
                        f"{row['BIC']:.1f}",
                        f"{row['RMSE']:.6f}"
                    )
                
                console.print(table)
                
                # Best model
                best_model = comparison.iloc[0]['Model']
                console.print(f"\n  [green]✓ Best model (by AIC): {best_model}[/green]")
            except Exception as e:
                console.print(f"  [yellow]Could not compare models: {e}[/yellow]")
                comparison = None
                # Fallback: display fitted parameters
                console.print("\n[yellow]Quadratic fit parameters (diagnostic):[/yellow]")
                console.print(f"  γ₀ = {fit_result['gamma0']:.6f} (target: 0.834)")
                console.print(f"  γ₁ = {fit_result['gamma1']:.6f} (target: 0.108)")
                console.print(f"  γ₂ = {fit_result['gamma2']:.8f} (target: 0.0105627)")
            
            # Topological constraint check
            console.print("\n[yellow]Topological constraint check:[/yellow]")
            # Use anchored γ₀=0.834 (consistent with normalization)
            gamma2_expected = 0.834 / (8 * np.pi**2)
            gamma2_actual = fit_result['gamma2']
            deviation = abs(gamma2_actual - gamma2_expected) / gamma2_expected * 100
            console.print(f"  γ₂ expected (γ₀/8π²): {gamma2_expected:.8f}")
            console.print(f"  γ₂ from quadratic fit: {gamma2_actual:.8f}")
            if deviation > 50:
                console.print(f"  [red]✗ Deviation: {deviation:.1f}% (constraint failed)[/red]")
            else:
                console.print(f"  [yellow]⚠ Deviation: {deviation:.1f}%[/yellow]")
            
            # Calibration-free test from data
            try:
                console.print("\n[yellow]Calibration-free φ test (from data):[/yellow]")
                from .verify import phi_third_diff_from_data
                phi_test = phi_third_diff_from_data(gamma)
                console.print(f"  {phi_test['message']}")
            except Exception as e:
                console.print(f"  [yellow]Phi test skipped: {e}[/yellow]")
            
            # Rolling window analysis
            rolling = verify_lnD_cubic_rolling(chain, window=8)
            console.print(f"  Rolling window: {rolling['message']}")
            
            # Save chain
            out.parent.mkdir(parents=True, exist_ok=True)
            save_chain_func(chain, out)
            console.print(f"\n  ✓ Chain saved to: [cyan]{out}[/cyan]")
            
            # Save summary JSON
            summary = {
                "chain": {
                    "n_orbits": len(chain),
                    "D_range": [int(chain['D'].max()), int(chain['D'].min())],
                    "labels": chain["label"].tolist()
                },
                "normalization": {
                    "s_star": s_star,
                    "lambda": lam,
                    "gamma0": gamma[0]
                },
                "quadratic_diagnostics": {
                    "r2": diag['r2'],
                    "coefficients": diag['coef'].tolist(),
                    "d3_mean": diag['d3_mean'],
                    "d3_var": diag['d3_var']
                },
                "fit_quality": {
                    "gamma0": fit_result['gamma0'],
                    "gamma1": fit_result['gamma1'],
                    "gamma2": fit_result['gamma2'],
                    "r_squared": fit_result['r_squared'],
                    "quality": quality['quality']
                },
                "model_comparison": comparison.to_dict('records') if 'comparison' in locals() else None
            }
            
            summary_path = output_dir / "e8_chain_summary.json"
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2)
            console.print(f"  ✓ Summary saved to: [cyan]{summary_path}[/cyan]")
            
            chain_result = chain
            
        except Exception as e:
            console.print(f"[red]Error in E8 chain task: {e}[/red]")
            if not run_rge:
                raise typer.Exit(1)
    
    # Task 2: Verify RGE crossings
    if run_rge:
        console.print("\n[yellow]Task 2:[/yellow] Verifying RGE crossings...")
        
        # Set default path
        if rge is None:
            rge = Path("data/gauge_couplings.csv")
        
        try:
            # Load RGE data
            df_rge = load_rge(rge)
            console.print(f"  ✓ Loaded {len(df_rge)} RGE points")
            console.print(f"  μ range: {df_rge['mu_GeV'].min():.1e} - {df_rge['mu_GeV'].max():.1e} GeV")
            
            # Find crossings
            crossings = verify_rge_crossings(df_rge)
            
            # Display results
            console.print("\n[green]Crossing scales:[/green]")
            for group in ["E6", "E7", "E8"]:
                if crossings[group]["found"]:
                    mu = crossings[group]["mu_GeV"]
                    target = crossings[group]["target"]
                    console.print(f"  {group}: μ = {mu:.3e} GeV (α₃ = 1/{group[1]}π = {target:.8f})")
                else:
                    console.print(f"  {group}: [red]Not found[/red]")
            
            # Save RGE checks
            rge_path = output_dir / "rge_checks.json"
            with open(rge_path, 'w') as f:
                # Convert to serializable format
                rge_summary = {}
                for group, data in crossings.items():
                    rge_summary[group] = {
                        "target_alpha": data["target"],
                        "mu_GeV": data["mu_GeV"] if data["found"] else None,
                        "log10_mu": data.get("log10_mu", None) if data["found"] else None,
                        "found": data["found"]
                    }
                json.dump(rge_summary, f, indent=2)
            console.print(f"\n  ✓ RGE checks saved to: [cyan]{rge_path}[/cyan]")
            
        except Exception as e:
            console.print(f"[red]Error in RGE verification: {e}[/red]")
            if not run_e8:
                raise typer.Exit(1)
    
    # Generate comprehensive report if we have results
    if html and chain_result is not None and fit_result is not None:
        console.print("\n[yellow]Generating comprehensive report...[/yellow]")
        try:
            # Run verification
            verification = comprehensive_verification(chain_result, fit_result)
            
            # Generate report with plots
            report_path = generate_report(
                chain_result, 
                fit_result, 
                verification, 
                output_dir
            )
            console.print(f"  ✓ Report saved to: [cyan]{report_path}[/cyan]")
            
        except Exception as e:
            console.print(f"[yellow]Warning: Could not generate full report: {e}[/yellow]")
    
    console.print("\n" + "=" * 50)
    console.print("[bold green]✓ Analysis complete![/bold green]")


@app.command()
def fit(
    data_file: Path = typer.Argument(
        ...,
        help="Path to orbit data file (CSV)"
    ),
    max_len: Optional[int] = typer.Option(
        None,
        "--max-len",
        help="Maximum chain length"
    ),
    gamma0: float = typer.Option(
        0.834,
        "--gamma0",
        help="Target value for γ₀"
    ),
    output_dir: Optional[Path] = typer.Option(
        Path("results"),
        "--output",
        "-o",
        help="Output directory for results"
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Verbose output"
    )
):
    """
    Legacy fit command - redirects to solve --rebuild-e8.
    """
    console.print("[yellow]Note: 'fit' command is legacy. Using 'solve --rebuild-e8' instead.[/yellow]")
    
    # Call solve with appropriate options
    solve(
        rebuild_e8=True,
        orbits=data_file,
        output_dir=output_dir,
        all=False
    )


@app.command()
def verify(
    chain_file: Path = typer.Argument(
        ...,
        help="Path to chain CSV file"
    ),
    gamma0: float = typer.Option(
        0.834,
        help="γ₀ value"
    ),
    gamma1: float = typer.Option(
        0.108,
        help="γ₁ value"
    ),
    gamma2: float = typer.Option(
        0.0105627,
        help="γ₂ value"
    )
):
    """
    Verify calibration-free consistency tests with given parameters.
    """
    import pandas as pd
    from .verify import verify_lnD_cubic, phi_third_diff_from_fit
    
    console.print("\n[bold cyan]Verification Mode[/bold cyan]")
    console.print("=" * 50)
    
    # Load chain
    chain = pd.read_csv(chain_file)
    console.print(f"Loaded chain with {len(chain)} entries")
    
    # Test 1: ln(D) cubic
    console.print("\n[yellow]Test 1:[/yellow] ln(D) cubic pattern")
    lnD_result = verify_lnD_cubic(chain)
    console.print(f"  {lnD_result['message']}")
    
    # Test 2: phi third difference
    console.print("\n[yellow]Test 2:[/yellow] φ third difference")
    phi_result = phi_third_diff_from_fit(gamma0, gamma1, gamma2)
    console.print(f"  {phi_result['message']}")
    
    # Test 3: Normalization
    console.print("\n[yellow]Test 3:[/yellow] Normalization checks")
    console.print(f"  γ₂ theoretical = γ₀/(8π²) = {gamma0/(8*np.pi**2):.8f}")
    console.print(f"  γ₂ provided = {gamma2:.8f}")
    console.print(f"  Deviation: {abs(gamma2 - gamma0/(8*np.pi**2))/(gamma0/(8*np.pi**2))*100:.2f}%")
    
    if lnD_result['valid'] and phi_result['valid']:
        console.print("\n[green]✓ All tests passed![/green]")
    else:
        console.print("\n[yellow]⚠ Some tests failed[/yellow]")


@app.command()
def info():
    """
    Display information about the E8 Orbit Engine.
    """
    console.print("\n[bold cyan]E8 Orbit Engine v0.2.0[/bold cyan]")
    console.print("=" * 50)
    console.print("""
This engine derives the quadratic damping function γ(n) from
nilpotent E8 orbits and validates theoretical predictions.

Key features:
- Data-driven derivation of γ(n) ≈ 0.834 + 0.108n + 0.0105627n²
- Calibration-free consistency test: Δ³ln(φ) = -2γ₂
- Topological constraint verification: γ₂ = γ₀/(8π²)
- E6/E7/E8 crossing scale verification from RGE data
- Comprehensive visualization and reporting

Theory reference:
"E8 Lattice Structure and Fine Structure Constants"
Stefan Hamann, 2025

Commands:
  solve   - Main command with --rebuild-e8, --verify-rge, --all options
  fit     - Legacy: Fit quadratic damping from orbit data
  verify  - Run verification tests with given parameters
  info    - Display this information
    """)


if __name__ == "__main__":
    app()