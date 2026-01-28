#!/usr/bin/env python3
"""Debug-Script für RGE Beta-Funktionen"""

import numpy as np
from RGEs import *

# Test-Werte
g1 = 0.357 * np.sqrt(3/5)  # GUT norm.
g2 = 0.652
g3 = 1.221

# Dummy Yukawas
Yu = np.zeros((3,3))
Yd = np.zeros((3,3))
Ye = np.zeros((3,3))
yN = np.zeros((3,3))

print("Debug Beta-Funktionen")
print("=====================\n")

# Teste beta_gauge_only Funktion
def beta_gauge_only(t, y):
    mu = np.exp(t)
    g1, g2, g3 = y
    
    try:
        # Verwende die generierten 1-Loop Beta-Funktionen
        beta_g1_val = (3/5) * beta_g1(1, g1/np.sqrt(3/5), g2, g3, Yu, Yd, Ye, yN)
        beta_g2_val = beta_g2(1, g2, g1/np.sqrt(3/5), g3, Yu, Yd, Ye, yN)
        beta_g3_val = beta_g3(1, g3, g1/np.sqrt(3/5), g2, Yu, Yd)
        
        print(f"Bei mu = {mu:.2e} GeV:")
        print(f"  beta_g1 = {beta_g1_val}")
        print(f"  beta_g2 = {beta_g2_val}")
        print(f"  beta_g3 = {beta_g3_val}")
        
        return [beta_g1_val, beta_g2_val, beta_g3_val]
    except Exception as e:
        print(f"Fehler in beta_gauge_only: {e}")
        import traceback
        traceback.print_exc()
        return [0, 0, 0]

# Teste bei verschiedenen Skalen
for log_mu in [2, 10, 19]:
    print(f"\nTest bei log(mu) = {log_mu}:")
    result = beta_gauge_only(log_mu, [g1, g2, g3])
    print(f"  Ergebnis: {result}") 