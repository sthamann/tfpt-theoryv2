import sys
import importlib
import numpy as np
from pathlib import Path
from scipy.integrate import ode

from .constants import (
    M_Pl, M_GUT, M_Z, M_SIGMA, M_NR, M_PHI,
    BETA_SM, BETA_E8, BETA_ABOVE_SIGMA, BETA_ABOVE_NR,
    gravity_coefficient, get_beta_coefficients
)

class E8CascadeSolver:
    """Clean interface to the E8 Cascade model RGE solver"""
    
    def __init__(self, model_package, model_module, results_dir, *,
                 enable_gravity_portal=False):
        """
        Initialize the E8 Cascade solver
        
        Parameters:
        -----------
        model_package : str
            Directory containing the PyR@TE-generated model
        model_module : str  
            Module name (without .py extension)
        results_dir : str/Path
            Output directory for results
        enable_gravity_portal : bool
            Enable α³ gravity corrections (experimental)
        """
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        
        # Add model path to Python path if needed
        model_path = Path(model_package).resolve()
        if str(model_path) not in sys.path:
            sys.path.insert(0, str(model_path))
        
        # Import the generated model  
        self.mod = importlib.import_module(model_module)
        self.model = self.mod.RGEsolver('rge', tmin=2, tmax=19, initialScale=2)
        
        # Set physical initial conditions at M_Z
        self._set_initial_conditions()
        
        # Set gravity portal coefficient (keeping sign!)
        self.model.cR3.initialValue = gravity_coefficient(BETA_E8['g3'])
        print(f"c_R3 = {self.model.cR3.initialValue:.6e} (b₃ = {BETA_E8['g3']}, M_GUT = {M_GUT:.1e} GeV)")
        
        self.enable_grav = enable_gravity_portal
        self.cR_factor = 1.0  # Can be modified for parameter studies
    
    def _set_initial_conditions(self):
        """Set realistic initial values at M_Z"""
        # Gauge couplings
        # g1 in SM normalization (g_Y), will convert to GUT norm internally
        self.model.g1.initialValue = 0.357  # g_Y at M_Z
        self.model.g2.initialValue = 0.652
        self.model.g3.initialValue = 1.221
        
        # Yukawa couplings (only 3rd generation)
        self.model.Yu.initialValue = [[0., 0., 0.],
                                      [0., 0., 0.],
                                      [0., 0., 0.94]]  # yt at M_Z
        
        self.model.Yd.initialValue = [[0., 0., 0.],
                                      [0., 0., 0.],
                                      [0., 0., 0.017]]  # yb at M_Z
        
        self.model.Ye.initialValue = [[0., 0., 0.],
                                      [0., 0., 0.],
                                      [0., 0., 0.010]]  # ytau at M_Z
        
        self.model.yN.initialValue = [[0., 0., 0.],
                                      [0., 0., 0.],
                                      [0., 0., 0.]]  # No RH neutrinos yet
        
        # Scalar couplings
        self.model.lambda_.initialValue = 0.13  # Higgs self-coupling
        self.model.lPhi.initialValue = 0
        self.model.lHphi.initialValue = 0
        
        # Mass parameters
        self.model.mu2.initialValue = 0
        self.model.MPhi.initialValue = 0
        
        # VEVs
        self.model.vSM.initialValue = 246.0  # Electroweak VEV
        self.model.vPQ.initialValue = 0
    
    def _beta_alpha_grav(self, g, c_R):
        """
        Gravity-induced correction to gauge coupling beta function
        
        The correction is: Δβ(g_i) = c_Ri * g_i³
        where c_Ri = b_i/(16π²) * (M_GUT/M_Pl)²
        """
        return c_R * g**3
    
    def solve(self):
        """Solve the RGE system"""
        if self.enable_grav:
            # Patch beta functions to include gravity corrections
            original_beta = self.model.betaFunction
            
            # Get robust indices for gauge couplings
            # This works even if parameter ordering changes
            all_params = [c.name for cList in self.model.couplings.values() 
                         for c in cList if not c.is_matrix]
            all_params += [c.name+'_{'+str(i)+str(j)+'}' 
                          for cList in self.model.couplings.values() 
                          for c in cList if c.is_matrix
                          for i in range(1, c.shape[0]+1)
                          for j in range(1, c.shape[1]+1)]
            
            try:
                idx_g1 = all_params.index('g1')
                idx_g2 = all_params.index('g2')
                idx_g3 = all_params.index('g3')
            except ValueError as e:
                raise ValueError(f"Could not find gauge coupling indices: {e}")
            
            gauge_indices = [idx_g1, idx_g2, idx_g3]
            gauge_names = ['g1', 'g2', 'g3']
            
            def patched_beta(t, y):
                dy = original_beta(t, y)
                mu = 10**t  # t = log10(mu)
                
                # Get appropriate beta coefficients for this scale
                beta_coeffs = get_beta_coefficients(mu)
                
                # Add gravity corrections with threshold-dependent coefficients
                for idx, name in zip(gauge_indices, gauge_names):
                    g = y[idx]
                    b = beta_coeffs[name]
                    c_R = gravity_coefficient(b) * self.cR_factor
                    dy[idx] += self._beta_alpha_grav(g, c_R)
                
                return dy
            
            self.model.betaFunction = patched_beta
            
            # Print info at M_GUT scale
            beta_gut = get_beta_coefficients(M_GUT)
            print("✓ Gravity portal enabled (g³ corrections active)")
            print(f"  At M_GUT = {M_GUT:.1e} GeV:")
            for name in gauge_names:
                c_R = gravity_coefficient(beta_gut[name]) * self.cR_factor
                print(f"  c_R({name}) = {c_R:+.6e} (b = {beta_gut[name]:.3f})")
        
        print(f"Solving RGEs from log₁₀(M_Z) = 2 to log₁₀(M_Pl) = 19...")
        print(f"  Thresholds: Σ_F at {M_SIGMA:.0e} GeV, N_R at {M_NR:.0e} GeV, Φ at {M_PHI:.0e} GeV")
        # Use adaptive step size for better efficiency
        # PyR@TE's solve() uses scipy.ode which handles adaptive stepping
        self.model.solve(step=0.1)  # Initial step, will be adapted
        print("✓ RGE solution complete")
        
        # Save results
        self._save_results()
    
    def _save_results(self):
        """Save numerical results to CSV"""
        import pandas as pd
        
        # Extract results
        t_values = self.model.tList
        mu_values = 10**t_values
        
        # Gauge couplings
        g1_SM = self.model.solutions['g1']  # SM normalization (g_Y)
        g2_values = self.model.solutions['g2'] 
        g3_values = self.model.solutions['g3']
        
        # Convert g1 to GUT normalization
        GUT_NORM = np.sqrt(5/3)
        g1_GUT = GUT_NORM * g1_SM
        
        # Compute alpha_i = g_i²/(4π)
        alpha1 = g1_GUT**2 / (4*np.pi)  # Using GUT normalization!
        alpha2 = g2_values**2 / (4*np.pi)
        alpha3 = g3_values**2 / (4*np.pi)
        
        # Create DataFrame
        df = pd.DataFrame({
            'mu_GeV': mu_values,
            'log10_mu': t_values,
            'g1_SM': g1_SM,
            'g1_GUT': g1_GUT,
            'g2': g2_values,
            'g3': g3_values,
            'alpha1_GUT': alpha1,
            'alpha2': alpha2,
            'alpha3': alpha3,
            'alpha1_inv_GUT': 1/alpha1,
            'alpha2_inv': 1/alpha2,
            'alpha3_inv': 1/alpha3
        })
        
        csv_path = self.results_dir / "gauge_couplings.csv"
        df.to_csv(csv_path, index=False)
        print(f"✓ Results saved to {csv_path}")
        
        # Check for unification
        idx_min = np.argmin(np.abs(1/alpha1 - 1/alpha2))
        mu_unif = mu_values[idx_min]
        print(f"\nUnification check:")
        print(f"  Closest approach at μ = {mu_unif:.2e} GeV")
        print(f"  α₁⁻¹ (GUT) = {1/alpha1[idx_min]:.2f}")
        print(f"  α₂⁻¹ = {1/alpha2[idx_min]:.2f}")
        print(f"  |Δ| = {abs(1/alpha1[idx_min] - 1/alpha2[idx_min]):.2f}")
    
    def make_plots(self):
        """Generate standard plots"""
        import matplotlib.pyplot as plt
        
        # Create our own gauge running plot instead of using the built-in
        # which has LaTeX rendering issues
        self._plot_gauge_running()
        
        # Also make a custom unification plot
        self._plot_unification()
    
    def _plot_gauge_running(self):
        """Create gauge coupling running plot"""
        import matplotlib.pyplot as plt
        
        mu = 10**self.model.tList
        g1_SM = self.model.solutions['g1']
        g2 = self.model.solutions['g2']
        g3 = self.model.solutions['g3']
        
        # Convert to GUT normalization
        GUT_NORM = np.sqrt(5/3)
        g1_GUT = GUT_NORM * g1_SM
        
        # Compute alpha values
        alpha1 = g1_GUT**2 / (4*np.pi)
        alpha2 = g2**2 / (4*np.pi)
        alpha3 = g3**2 / (4*np.pi)
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
        
        # Top panel: Gauge couplings
        ax1.plot(mu, g1_GUT, 'b-', linewidth=2, label=r'$g_1$ (GUT)')
        ax1.plot(mu, g2, 'r-', linewidth=2, label=r'$g_2$')
        ax1.plot(mu, g3, 'g-', linewidth=2, label=r'$g_3$')
        ax1.set_xscale('log')
        ax1.set_ylabel(r'$g_i$', fontsize=14)
        ax1.set_title('E8 Cascade: Gauge Coupling Evolution', fontsize=16)
        ax1.legend(fontsize=12)
        ax1.grid(True, alpha=0.3)
        
        # Bottom panel: Inverse alphas
        ax2.plot(mu, 1/alpha1, 'b-', linewidth=2, label=r'$\alpha_1^{-1}$ (GUT)')
        ax2.plot(mu, 1/alpha2, 'r-', linewidth=2, label=r'$\alpha_2^{-1}$')
        ax2.plot(mu, 1/alpha3, 'g-', linewidth=2, label=r'$\alpha_3^{-1}$')
        ax2.set_xscale('log')
        ax2.set_xlabel(r'$\mu$ [GeV]', fontsize=14)
        ax2.set_ylabel(r'$\alpha_i^{-1}$', fontsize=14)
        ax2.set_ylim(0, 70)
        ax2.legend(fontsize=12)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        plot_path = self.results_dir / "gauge_running.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"✓ Plot saved to {plot_path}")
        plt.close()
    
    def _plot_unification(self):
        """Create a dedicated gauge unification plot"""
        import matplotlib.pyplot as plt
        
        mu = 10**self.model.tList
        g1_SM = self.model.solutions['g1']  # SM normalization
        g2 = self.model.solutions['g2']
        g3 = self.model.solutions['g3']
        
        # Convert to GUT normalization
        GUT_NORM = np.sqrt(5/3)
        g1_GUT = GUT_NORM * g1_SM
        
        # Compute inverse fine structure constants
        alpha1_inv = (4*np.pi) / g1_GUT**2  # GUT normalized!
        alpha2_inv = (4*np.pi) / g2**2
        alpha3_inv = (4*np.pi) / g3**2
        
        plt.figure(figsize=(10, 8))
        
        plt.plot(mu, alpha1_inv, 'b-', linewidth=2, label=r'$\alpha_1^{-1}$ (GUT)')
        plt.plot(mu, alpha2_inv, 'r-', linewidth=2, label=r'$\alpha_2^{-1}$')
        plt.plot(mu, alpha3_inv, 'g-', linewidth=2, label=r'$\alpha_3^{-1}$')
        
        plt.xscale('log')
        plt.xlabel(r'$\mu$ [GeV]', fontsize=14)
        plt.ylabel(r'$\alpha_i^{-1}$', fontsize=14)
        plt.title('E8 Cascade: Gauge Coupling Evolution', fontsize=16)
        plt.ylim(0, 70)
        plt.grid(True, alpha=0.3)
        
        # Mark important scales
        plt.axvline(91.2, color='k', linestyle=':', alpha=0.5, label='M_Z')
        plt.axvline(1e16, color='orange', linestyle='--', alpha=0.5, label='GUT')
        plt.axvline(1.22e19, color='purple', linestyle='--', alpha=0.5, label='M_Pl')
        
        plt.legend(loc='best', fontsize=12)
        plt.tight_layout()
        
        plot_path = self.results_dir / "unification.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"✓ Unification plot saved to {plot_path}")
        plt.close() 