#!/usr/bin/env python3
"""
Vereinfachter Test der E8 Cascade RGEs
=======================================
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Import der generierten RGEs
from RGEs import *

# Physikalische Konstanten
M_Z = 91.1876  # GeV
M_Pl = 1.22e19  # GeV
v_EW = 246.0   # GeV
GUT_NORM = np.sqrt(3/5)

def run_simple_test():
    """Einfacher Test der RGE-Evolution"""
    
    print("E8 Cascade RGE Test")
    print("===================\n")
    
    # Anfangswerte bei M_Z (nur Eichkopplungen)
    g1_MZ = 0.357  # bereits SM normalisiert (keine weitere GUT-Norm nötig!)
    g2_MZ = 0.652
    g3_MZ = 1.221
    
    print(f"Startwerte bei M_Z = {M_Z} GeV:")
    print(f"  g1 = {g1_MZ:.4f} (SM norm.)")
    print(f"  g2 = {g2_MZ:.4f}")
    print(f"  g3 = {g3_MZ:.4f}")
    
    # Beta-Funktionen mit generierten Funktionen
    def beta_gauge_only(t, y):
        mu = np.exp(t)
        g1, g2, g3 = y
        
        # Dummy Yukawas (alle null für reine Eichevolution)
        Yu = np.zeros((3,3))
        Yd = np.zeros((3,3))
        Ye = np.zeros((3,3))
        yN = np.zeros((3,3))
        
        # Verwende die generierten 1-Loop Beta-Funktionen (g1 bereits SM-normalisiert!)
        beta_g1_val = beta_g1(1, g1, g2, g3, Yu, Yd, Ye, yN)
        beta_g2_val = beta_g2(1, g2, g1, g3, Yu, Yd, Ye, yN)
        beta_g3_val = beta_g3(1, g3, g1, g2, Yu, Yd)
        
        # Konvertiere zu Standard-Python-Floats
        return [float(beta_g1_val), float(beta_g2_val), float(beta_g3_val)]
    
    # Löse RGEs
    t_span = [np.log(M_Z), np.log(M_Pl)]
    y0 = [g1_MZ, g2_MZ, g3_MZ]
    
    print("\nLöse RGEs von M_Z bis M_Planck...")
    try:
        # Verwende RK45 statt DOP853 für bessere Stabilität
        sol = solve_ivp(beta_gauge_only, t_span, y0, 
                        method='RK45', rtol=1e-6, atol=1e-8,
                        dense_output=True, max_step=0.5)
    except Exception as e:
        print(f"Fehler bei solve_ivp: {e}")
        import traceback
        traceback.print_exc()
        sol = None
    
    if sol is not None and sol.success:
        print("✓ RGE-Evolution erfolgreich!")
        
        # Werte bei M_Planck
        g1_Pl, g2_Pl, g3_Pl = sol.y[:, -1]
        print(f"\nWerte bei M_Planck = {M_Pl:.2e} GeV:")
        print(f"  g1 = {g1_Pl:.4f}")
        print(f"  g2 = {g2_Pl:.4f}")
        print(f"  g3 = {g3_Pl:.4f}")
        
        # Plotte Evolution
        mu = np.logspace(2, 19, 1000)
        t = np.log(mu)
        y = sol.sol(t)
        
        plt.figure(figsize=(10, 6))
        
        # Berechne α_i = g_i²/(4π)
        alpha1 = y[0]**2 / (4*np.pi)
        alpha2 = y[1]**2 / (4*np.pi)
        alpha3 = y[2]**2 / (4*np.pi)
        
        plt.plot(mu, 1/alpha1, 'b-', label=r'$\alpha_1^{-1}$ (SM)', linewidth=2)
        plt.plot(mu, 1/alpha2, 'r-', label=r'$\alpha_2^{-1}$', linewidth=2)
        plt.plot(mu, 1/alpha3, 'g-', label=r'$\alpha_3^{-1}$', linewidth=2)
        
        plt.xscale('log')
        plt.xlabel(r'$\mu$ [GeV]', fontsize=14)
        plt.ylabel(r'$\alpha_i^{-1}$', fontsize=14)
        plt.title('E8 Cascade: Eichkopplungs-Evolution (1-Loop)', fontsize=16)
        plt.ylim(0, 70)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=12)
        
        # Markiere wichtige Skalen
        plt.axvline(M_Z, color='k', linestyle=':', alpha=0.5, label='M_Z')
        plt.axvline(1e16, color='orange', linestyle='--', alpha=0.5, label='GUT')
        plt.axvline(M_Pl, color='purple', linestyle='--', alpha=0.5, label='M_Pl')
        
        plt.tight_layout()
        plt.savefig('simple_gauge_evolution.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Teste Vereinigung
        print("\nVereinigungs-Test:")
        # Finde wo g1 ≈ g2
        diff_12 = np.abs(1/alpha1 - 1/alpha2)
        idx_min = np.argmin(diff_12)
        mu_unif = mu[idx_min]
        print(f"  Minimaler |α₁⁻¹ - α₂⁻¹| bei μ = {mu_unif:.2e} GeV")
        print(f"  α₁⁻¹ = {1/alpha1[idx_min]:.2f}, α₂⁻¹ = {1/alpha2[idx_min]:.2f}")
        
    else:
        print("✗ Fehler bei RGE-Evolution!")
        
    # Teste cR3-Effekt
    print("\n\nTeste Gravitationseffekt (cR3):")
    print("================================")
    
    # Dummy-Werte für Yukawas (nur 3. Generation)
    Yu = np.zeros((3,3))
    Yu[2,2] = 0.95  # top
    Yd = np.zeros((3,3))
    Yd[2,2] = 0.024  # bottom
    Ye = np.zeros((3,3))
    Ye[2,2] = 0.010  # tau
    yN = np.zeros((3,3))
    
    # Teste beta_cR3
    cR3 = 0.01
    lambda_ = 0.13
    lHphi = 0.01
    
    beta_cR3_1loop = beta_cR3(1, g1_MZ, g2_MZ, Yu, Yd, Ye, yN, 
                              lambda_, cR3, g3_MZ, lHphi)
    beta_cR3_2loop = beta_cR3(2, g1_MZ, g2_MZ, Yu, Yd, Ye, yN, 
                              lambda_, cR3, g3_MZ, lHphi)
    
    print(f"  cR3 = {cR3}")
    print(f"  β(cR3) 1-Loop = {beta_cR3_1loop:.6f}")
    print(f"  β(cR3) 2-Loop = {beta_cR3_2loop:.6f}")
    print(f"  → cR3 {'wächst' if beta_cR3_1loop > 0 else 'fällt'} mit der Skala")


if __name__ == "__main__":
    run_simple_test() 