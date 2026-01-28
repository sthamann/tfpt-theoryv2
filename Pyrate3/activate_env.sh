#!/bin/bash
# Aktiviere das virtuelle Python-Environment
source venv/bin/activate
echo "Virtual environment aktiviert. Python Version:"
python --version
echo "Installierte Pakete:"
pip list | grep -E "(PyYAML|sympy|h5py|numpy|scipy|matplotlib)"
