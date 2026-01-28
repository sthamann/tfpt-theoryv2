#!/usr/bin/env python3
"""
E8 Cascade Fixpoint Analysis
============================
Analyse der Fixpunkte im E8-Cascade Modell mit 2-Loop RGEs
und gravitativem α³-Effekt durch cR3 Portal
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# Import der generierten RGEs
from RGEs import *

# Physikalische Konstanten
M_Z = 91.1876  # GeV
M_Pl = 1.22e19  # GeV
v_EW = 246.0   # GeV

# Schwellenwerte (GeV)
M_SIGMA = 1e3    # TeV-Skala Triplet
M_NR = 1e15      # Seesaw-Skala
M_PHI = 1e16     # PQ/Axion-Skala

# GUT-Normalisierung für Hypercharge
GUT_NORM = np.sqrt(3/5)

class E8CascadeModel:
    """E8 Cascade Modell mit Fixpunkt-Analyse"""
    
    def __init__(self, loops=2):
        self.loops = loops
        self.n_couplings = 23  # Anzahl der Kopplungen
        
        # Kopplungs-Indices
        self.idx = {
            'g1': 0, 'g2': 1, 'g3': 2,
            'Yu': slice(3, 12), 'Yd': slice(12, 21), 
            'Ye': slice(21, 30), 'yN': slice(30, 39),
            'lambda': 39, 'lPhi': 40, 'lHphi': 41,
            'cR3': 42, 'mu2': 43, 'MPhi': 44,
            'vSM': 45, 'vPQ': 46
        }
        
    def beta_functions(self, t, y, apply_thresholds=True):
        """Berechne alle Beta-Funktionen"""
        mu = np.exp(t)  # t = log(mu)
        
        # Extrahiere Kopplungen
        g1, g2, g3 = y[0:3]
        Yu = y[3:12].reshape(3,3)
        Yd = y[12:21].reshape(3,3)
        Ye = y[21:30].reshape(3,3)
        yN = y[30:39].reshape(3,3)
        lambda_, lPhi, lHphi = y[39:42]
        cR3 = y[42]
        mu2, MPhi = y[43:45]
        vSM, vPQ = y[45:47]
        
        # Beta-Funktionen
        betas = np.zeros_like(y)
        
        # Eichkopplungen (mit GUT-Normalisierung für g1)
        betas[0] = (3/5) * beta_g1(self.loops, g1/GUT_NORM, g2, g3, Yu, Yd, Ye, yN)
        betas[1] = beta_g2(self.loops, g2, g1/GUT_NORM, g3, Yu, Yd, Ye, yN)
        betas[2] = beta_g3(self.loops, g3, g1/GUT_NORM, g2, Yu, Yd)
        
        # Yukawa-Kopplungen
        betas[3:12] = beta_Yu(self.loops, g1/GUT_NORM, g2, g3, Yu, Yd, Ye, yN, lambda_, lHphi).flatten()
        betas[12:21] = beta_Yd(self.loops, g1/GUT_NORM, g2, g3, Yu, Yd, Ye, yN, lambda_, lHphi).flatten()
        betas[21:30] = beta_Ye(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, g3, lambda_, lHphi).flatten()
        betas[30:39] = beta_yN(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, g3, lambda_, lHphi).flatten()
        
        # Skalare Kopplungen
        betas[39] = beta_lambda_(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, lambda_, lHphi, g3)
        betas[40] = beta_lPhi(self.loops, lPhi, lHphi, g1/GUT_NORM, g2, Yu, Yd, Ye, yN)
        betas[41] = beta_lHphi(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, lambda_, lPhi, lHphi, g3)
        
        # Gravity Portal
        betas[42] = beta_cR3(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, lambda_, cR3, g3, lHphi)
        
        # Massenterme
        betas[43] = beta_mu2(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, lambda_, lHphi, cR3, mu2, MPhi, g3)
        betas[44] = beta_MPhi(self.loops, lPhi, lHphi, mu2, MPhi, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, cR3)
        
        # VEVs
        xiGauge = 1.0  # Landau gauge
        betas[45] = beta_vSM(self.loops, g1/GUT_NORM, g2, Yu, Yd, Ye, yN, vSM, xiGauge, g3, lambda_, lHphi)
        betas[46] = beta_vPQ(self.loops, lPhi, lHphi, vPQ)
        
        # Schwellen-Implementierung
        if apply_thresholds:
            if mu < M_SIGMA:
                # SigmaF ausschalten - keine Beiträge zu RGEs
                pass  # Implementierung hängt von Details ab
            if mu < M_NR:
                betas[30:39] = 0  # yN ausschalten
            if mu < M_PHI:
                betas[40] = 0  # lPhi ausschalten
                betas[41] = 0  # lHphi ausschalten
                betas[44] = 0  # MPhi ausschalten
                betas[46] = 0  # vPQ ausschalten
        
        return betas
    
    def find_fixed_points(self, initial_guess=None):
        """Finde Fixpunkte der RGEs"""
        if initial_guess is None:
            # Standard-Startpunkt
            initial_guess = np.zeros(47)
            initial_guess[0:3] = [0.5, 0.5, 0.5]  # Eichkopplungen
            initial_guess[42] = 0.01  # cR3
        
        def fixed_point_eq(y):
            """Fixpunkt-Gleichung: beta(g) = 0"""
            return self.beta_functions(np.log(M_Pl), y, apply_thresholds=False)
        
        # Löse Fixpunkt-Gleichungen
        result = fsolve(fixed_point_eq, initial_guess, full_output=True)
        
        if result[2] == 1:  # Konvergiert
            return result[0]
        else:
            print("Warnung: Fixpunkt-Suche nicht konvergiert!")
            return None
    
    def analyze_stability(self, fixed_point):
        """Analysiere Stabilität des Fixpunkts"""
        # Berechne Jacobi-Matrix numerisch
        eps = 1e-8
        n = len(fixed_point)
        jacobian = np.zeros((n, n))
        
        for i in range(n):
            y_plus = fixed_point.copy()
            y_minus = fixed_point.copy()
            y_plus[i] += eps
            y_minus[i] -= eps
            
            beta_plus = self.beta_functions(np.log(M_Pl), y_plus, apply_thresholds=False)
            beta_minus = self.beta_functions(np.log(M_Pl), y_minus, apply_thresholds=False)
            
            jacobian[:, i] = (beta_plus - beta_minus) / (2 * eps)
        
        # Eigenwerte bestimmen Stabilität
        eigenvalues = np.linalg.eigvals(jacobian)
        
        # Kritische Exponenten
        critical_exponents = -eigenvalues
        
        # Anzahl relevanter Richtungen (positive Realteile)
        n_relevant = np.sum(np.real(critical_exponents) > 0)
        
        return {
            'eigenvalues': eigenvalues,
            'critical_exponents': critical_exponents,
            'n_relevant': n_relevant,
            'stable': n_relevant == 0
        }
    
    def run_to_planck(self, initial_conditions):
        """RG-Evolution von M_Z zu M_Planck"""
        t_span = [np.log(M_Z), np.log(M_Pl)]
        
        sol = solve_ivp(
            lambda t, y: self.beta_functions(t, y),
            t_span, 
            initial_conditions,
            method='DOP853',
            rtol=1e-10,
            atol=1e-12,
            dense_output=True
        )
        
        return sol
    
    def plot_running(self, sol, title="RG Evolution"):
        """Plotte RG-Lauf der Kopplungen"""
        mu = np.logspace(2, 19, 1000)
        t = np.log(mu)
        y = sol.sol(t)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Eichkopplungen
        ax = axes[0, 0]
        ax.plot(mu, y[0], 'b-', label=r'$g_1$ (GUT norm.)')
        ax.plot(mu, y[1], 'r-', label=r'$g_2$')
        ax.plot(mu, y[2], 'g-', label=r'$g_3$')
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu$ [GeV]')
        ax.set_ylabel('Gauge couplings')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Yukawa-Kopplungen (3. Generation)
        ax = axes[0, 1]
        ax.plot(mu, y[11], 'b-', label=r'$y_t$')  # Yu[2,2]
        ax.plot(mu, y[20], 'r-', label=r'$y_b$')  # Yd[2,2]
        ax.plot(mu, y[29], 'g-', label=r'$y_\tau$')  # Ye[2,2]
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu$ [GeV]')
        ax.set_ylabel('Yukawa couplings')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Skalare Kopplungen
        ax = axes[1, 0]
        ax.plot(mu, y[39], 'b-', label=r'$\lambda$')
        ax.plot(mu, y[40], 'r-', label=r'$\lambda_\Phi$')
        ax.plot(mu, y[41], 'g-', label=r'$\lambda_{H\Phi}$')
        ax.plot(mu, y[42], 'm-', label=r'$c_{R^3}$', linewidth=2)
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu$ [GeV]')
        ax.set_ylabel('Scalar couplings')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Einheitliche Kopplung
        ax = axes[1, 1]
        # Berechne α_i = g_i²/(4π)
        alpha1 = y[0]**2 / (4*np.pi)
        alpha2 = y[1]**2 / (4*np.pi)
        alpha3 = y[2]**2 / (4*np.pi)
        ax.plot(mu, 1/alpha1, 'b-', label=r'$\alpha_1^{-1}$')
        ax.plot(mu, 1/alpha2, 'r-', label=r'$\alpha_2^{-1}$')
        ax.plot(mu, 1/alpha3, 'g-', label=r'$\alpha_3^{-1}$')
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu$ [GeV]')
        ax.set_ylabel(r'$\alpha_i^{-1}$')
        ax.set_ylim(0, 60)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Markiere Schwellen
        for ax_row in axes:
            for ax in ax_row:
                ax.axvline(M_SIGMA, color='k', linestyle='--', alpha=0.3, label='Σ_F')
                ax.axvline(M_NR, color='k', linestyle='-.', alpha=0.3, label='N_R')
                ax.axvline(M_PHI, color='k', linestyle=':', alpha=0.3, label='Φ')
        
        plt.suptitle(title)
        plt.tight_layout()
        return fig


def main():
    """Hauptanalyse"""
    print("E8 Cascade Fixpunkt-Analyse")
    print("=" * 40)
    
    # Initialisiere Modell
    model = E8CascadeModel(loops=2)
    
    # 1. Finde Fixpunkte
    print("\n1. Suche Fixpunkte...")
    fp = model.find_fixed_points()
    
    if fp is not None:
        print(f"   Fixpunkt gefunden!")
        print(f"   g1* = {fp[0]:.4f} (GUT norm.)")
        print(f"   g2* = {fp[1]:.4f}")
        print(f"   g3* = {fp[2]:.4f}")
        print(f"   cR3* = {fp[42]:.4f}")
        
        # 2. Stabilitätsanalyse
        print("\n2. Stabilitätsanalyse...")
        stability = model.analyze_stability(fp)
        print(f"   Anzahl relevanter Richtungen: {stability['n_relevant']}")
        print(f"   Fixpunkt ist {'stabil' if stability['stable'] else 'instabil'}")
        
        # Zeige größte kritische Exponenten
        ce_real = np.real(stability['critical_exponents'])
        ce_sorted = np.sort(ce_real)[::-1]
        print(f"   Größte kritische Exponenten: {ce_sorted[:5]}")
    
    # 3. RG-Evolution von M_Z
    print("\n3. RG-Evolution von M_Z...")
    
    # SM-Werte bei M_Z
    initial = np.zeros(47)
    initial[0] = 0.357 * GUT_NORM  # g1 (GUT normalisiert)
    initial[1] = 0.652              # g2
    initial[2] = 1.221              # g3
    initial[11] = 0.95              # yt
    initial[20] = 0.024             # yb
    initial[29] = 0.010             # ytau
    initial[39] = 0.130             # lambda
    initial[42] = 0.01              # cR3
    initial[45] = v_EW              # vSM
    
    # Löse RGEs
    sol = model.run_to_planck(initial)
    
    if sol.success:
        print("   RG-Evolution erfolgreich!")
        
        # Werte bei M_Planck
        y_planck = sol.y[:, -1]
        print(f"\n   Werte bei M_Planck:")
        print(f"   g1 = {y_planck[0]:.4f}")
        print(f"   g2 = {y_planck[1]:.4f}")
        print(f"   g3 = {y_planck[2]:.4f}")
        print(f"   cR3 = {y_planck[42]:.6f}")
        print(f"   λ = {y_planck[39]:.4f}")
        
        # Plotte Ergebnisse
        fig = model.plot_running(sol, "E8 Cascade RG Evolution (2-Loop)")
        plt.savefig('e8cascade_rg_evolution.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # 4. Scan über cR3
    print("\n4. cR3 Scan für Fixpunkt-Tuning...")
    cR3_values = np.linspace(0, 0.02, 21)
    unification_scales = []
    
    for cR3_init in cR3_values:
        initial[42] = cR3_init
        sol = model.run_to_planck(initial)
        
        if sol.success:
            # Finde Vereinigungsskala (wo g1 = g2)
            def diff_g1_g2(t):
                y = sol.sol(t)
                return y[0] - y[1]
            
            # Suche Nullstelle
            try:
                from scipy.optimize import brentq
                t_unif = brentq(diff_g1_g2, np.log(1e15), np.log(M_Pl))
                mu_unif = np.exp(t_unif)
                unification_scales.append(mu_unif)
            except:
                unification_scales.append(np.nan)
    
    # Plotte cR3-Abhängigkeit
    if len(unification_scales) > 0:
        unification_scales = np.array(unification_scales)
        plt.figure(figsize=(8, 6))
        valid_mask = ~np.isnan(unification_scales)
        if np.any(valid_mask):
            plt.plot(np.array(cR3_values)[valid_mask], 
                    unification_scales[valid_mask]/1e16, 
                    'b-', linewidth=2)
        plt.xlabel(r'$c_{R^3}$')
        plt.ylabel(r'$\mu_{\rm GUT}$ [$10^{16}$ GeV]')
        plt.title('Gravitativer Effekt auf GUT-Skala')
        plt.grid(True, alpha=0.3)
        plt.savefig('cR3_scan.png', dpi=300, bbox_inches='tight')
        plt.show()
    else:
        print("   Warnung: Keine Vereinigungsskalen gefunden!")
    
    print("\nAnalyse abgeschlossen!")


if __name__ == "__main__":
    main() 