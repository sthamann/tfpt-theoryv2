#!/usr/bin/env python3
"""Test der Beta-Funktionen"""

import numpy as np
from RGEs import *

# Test-Werte
g1 = 0.357  # SM norm.
g2 = 0.652
g3 = 1.221

# Dummy Yukawas
Yu = np.zeros((3,3))
Yd = np.zeros((3,3))
Ye = np.zeros((3,3))
yN = np.zeros((3,3))

print("Beta-Funktionen Test")
print("====================\n")

print(f"Input: g1={g1:.4f}, g2={g2:.4f}, g3={g3:.4f}")
print()

# 1-Loop
b1_1loop = beta_g1(1, g1, g2, g3, Yu, Yd, Ye, yN)
b2_1loop = beta_g2(1, g2, g1, g3, Yu, Yd, Ye, yN)
b3_1loop = beta_g3(1, g3, g1, g2, Yu, Yd)

print("1-Loop Beta-Funktionen:")
print(f"  β(g1) = {b1_1loop:.6f}")
print(f"  β(g2) = {b2_1loop:.6f}")
print(f"  β(g3) = {b3_1loop:.6f}")

# Beta-Koeffizienten
print("\nBeta-Koeffizienten bi = β(gi)/gi³:")
print(f"  b1 = {b1_1loop/g1**3:.3f} (erwartet: 41/6 ≈ 6.833)")
print(f"  b2 = {b2_1loop/g2**3:.3f} (erwartet: -11/6 ≈ -1.833)")
print(f"  b3 = {b3_1loop/g3**3:.3f} (erwartet: -7)")

# SM-Normalisierung sollte jetzt korrekte Werte liefern
print(f"\nHinweis: g1 ist bereits SM-normalisiert") 