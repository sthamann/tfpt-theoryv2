#!/bin/bash
# Helper script to run PyR@TE with virtual environment

# Check if a model file is provided
if [ $# -eq 0 ]; then
    echo "Verwendung: $0 <model_file> [additional_options]"
    echo ""
    echo "Verfügbare Model-Dateien:"
    echo "  Standard PyR@TE Models:"
    ls pyrate/models/*.model 2>/dev/null || echo "    Keine gefunden"
    echo ""
    echo "  Projekt-spezifische Models:"
    ls models/*.model 2>/dev/null || echo "    Keine gefunden"
    echo ""
    echo "Beispiel: $0 pyrate/models/SM.model -py -tex"
    exit 1
fi

# Activate virtual environment and run PyR@TE
source venv/bin/activate
cd pyrate
python pyR@TE.py -m "$@"
cd ..
