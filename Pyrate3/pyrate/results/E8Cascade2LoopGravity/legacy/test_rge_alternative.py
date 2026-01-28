#!/usr/bin/env python3
"""Alternative RGE Test mit manueller Integration"""

import numpy as np
import matplotlib.pyplot as plt
from RGEs import *

# Physikalische Konstanten
M_Z = 91.1876  # GeV
M_Pl = 1.22e19  # GeV
GUT_NORM = np.sqrt(3/5)

print("Alternative E8 Cascade RGE Test")
print("================================\n")

# Anfangswerte bei M_Z
g1_MZ = 0.357  # SM normalisiert
g2_MZ = 0.652
g3_MZ = 1.221

print(f"Startwerte bei M_Z = {M_Z} GeV:")
print(f"  g1 = {g1_MZ:.4f} (SM norm.)")
print(f"  g2 = {g2_MZ:.4f}")
print(f"  g3 = {g3_MZ:.4f}")

# Dummy Yukawas (alle null für reine Eichevolution)
Yu = np.zeros((3,3))
Yd = np.zeros((3,3))
Ye = np.zeros((3,3))
yN = np.zeros((3,3))

# Manuelle Euler-Integration
def manual_rge_evolution(g1_0, g2_0, g3_0, t_start, t_end, n_steps=1000):
    """Manuelle RGE-Evolution mit Euler-Methode"""
    
    t_values = np.linspace(t_start, t_end, n_steps)
    dt = (t_end - t_start) / (n_steps - 1)
    
    g1_vals = [g1_0]
    g2_vals = [g2_0]
    g3_vals = [g3_0]
    
    g1, g2, g3 = g1_0, g2_0, g3_0
    
    for i in range(1, n_steps):
        # Berechne Beta-Funktionen (g1 bereits SM-normalisiert)
        beta_g1_val = beta_g1(1, g1, g2, g3, Yu, Yd, Ye, yN)
        beta_g2_val = beta_g2(1, g2, g1, g3, Yu, Yd, Ye, yN)
        beta_g3_val = beta_g3(1, g3, g1, g2, Yu, Yd)
        
        # Euler-Schritt mit korrektem Faktor
        # PyR@TE RGEs gibt nackte Beta-Funktionen - wir müssen den Loop-Faktor selbst hinzufügen
        factor = 1/(16*np.pi**2) * np.log(10)  # log(10) weil t = log10(mu)
        g1 += dt * factor * beta_g1_val
        g2 += dt * factor * beta_g2_val
        g3 += dt * factor * beta_g3_val
        
        g1_vals.append(g1)
        g2_vals.append(g2)
        g3_vals.append(g3)
    
    return t_values, np.array(g1_vals), np.array(g2_vals), np.array(g3_vals)

# Führe manuelle Evolution durch
t_start = np.log10(M_Z)
t_end = np.log10(M_Pl)

print(f"\nFühre manuelle RGE-Evolution durch...")
t_vals, g1_vals, g2_vals, g3_vals = manual_rge_evolution(g1_MZ, g2_MZ, g3_MZ, t_start, t_end)

print("✓ Evolution erfolgreich!")

# Werte bei M_Planck
g1_Pl = g1_vals[-1]
g2_Pl = g2_vals[-1]
g3_Pl = g3_vals[-1]

print(f"\nWerte bei M_Planck = {M_Pl:.2e} GeV:")
print(f"  g1 = {g1_Pl:.4f}")
print(f"  g2 = {g2_Pl:.4f}")
print(f"  g3 = {g3_Pl:.4f}")

# Plotte Evolution
mu_vals = 10**t_vals

plt.figure(figsize=(10, 6))

# Berechne α_i = g_i²/(4π)
alpha1 = g1_vals**2 / (4*np.pi)
alpha2 = g2_vals**2 / (4*np.pi)
alpha3 = g3_vals**2 / (4*np.pi)

plt.plot(mu_vals, 1/alpha1, 'b-', label=r'$\alpha_1^{-1}$ (SM)', linewidth=2)
plt.plot(mu_vals, 1/alpha2, 'r-', label=r'$\alpha_2^{-1}$', linewidth=2)
plt.plot(mu_vals, 1/alpha3, 'g-', label=r'$\alpha_3^{-1}$', linewidth=2)

plt.xscale('log')
plt.xlabel(r'$\mu$ [GeV]', fontsize=14)
plt.ylabel(r'$\alpha_i^{-1}$', fontsize=14)
plt.title('E8 Cascade: Eichkopplungs-Evolution (1-Loop, manuelle Integration)', fontsize=16)
plt.ylim(0, 70)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)

# Markiere wichtige Skalen
plt.axvline(M_Z, color='k', linestyle=':', alpha=0.5, label='M_Z')
plt.axvline(1e16, color='orange', linestyle='--', alpha=0.5, label='GUT')
plt.axvline(M_Pl, color='purple', linestyle='--', alpha=0.5, label='M_Pl')

plt.tight_layout()
plt.savefig('manual_gauge_evolution.png', dpi=300, bbox_inches='tight')
plt.show()

# Teste Vereinigung
print("\nVereinigungs-Test:")
diff_12 = np.abs(1/alpha1 - 1/alpha2)
idx_min = np.argmin(diff_12)
mu_unif = mu_vals[idx_min]
print(f"  Minimaler |α₁⁻¹ - α₂⁻¹| bei μ = {mu_unif:.2e} GeV")
print(f"  α₁⁻¹ = {1/alpha1[idx_min]:.2f}, α₂⁻¹ = {1/alpha2[idx_min]:.2f}")
print(f"  Differenz: {diff_12[idx_min]:.4f}") 