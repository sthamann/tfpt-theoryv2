#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Erster PyR@TE-Testlauf: Toy-U(1)-Modell,
1-Loop-RGEs berechnen und auf die Konsole ausgeben.
"""

from pathlib import Path
from pyrate import RGEs

# --- 1. Modellpfad einlesen ----------------------------------------------
model_path = Path("models/toy_u1.model").resolve()

# --- 2. RGEs-Objekt anlegen -----------------------------------------------
# output_dir = 'out' legt TeX / Py-Dateien dort ab (optional)
rge = RGEs(str(model_path), loop_order=1, output_dir='out')

# --- 3. Berechnen ---------------------------------------------------------
rge.compute()             # erstellt Symbol-Ausdrücke

# --- 4. Ergebnisse ausgeben ----------------------------------------------
print("\n=== Gauge-Betafunktionen ===")
for label, beta in rge.betas['gauge'].items():
    print(f"β_{label} =", beta)

print("\n=== Quartik-Betafunktion ===")
print("β_lambda =", rge.betas['quartic']['lambda'])
