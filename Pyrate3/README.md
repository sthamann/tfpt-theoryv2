# Pyrate3 - RGE Computations for TFPT

This directory contains a local installation of **PyR@TE 3** (Python Renormalization Group Equations @ Two-loop-order and beyond) configured for TFPT (Topological Field Physics Theory) research.

## What is PyR@TE 3?

PyR@TE is a Python tool that computes **Renormalization Group Equations (RGEs)** for any renormalizable non-supersymmetric model. After defining gauge groups, particle content, and the Lagrangian in a model file, PyR@TE calculates the RGEs for all couplings at:

- **1-loop** level
- **2-loop** level
- **3-loop** level (gauge couplings only)

The core is based on [arXiv:1906.04625](https://arxiv.org/abs/1906.04625) for general gauge coupling RGEs up to 3-loop order.

## Directory Structure

```
Pyrate3/
├── README.md              # This file
├── requirements.txt       # Python dependencies
├── activate_env.sh        # Script to activate the virtual environment
├── run_pyrate.sh          # Helper script to run PyR@TE with a model
├── run_test.py            # Quick test script for Toy U(1) model
├── models/                # TFPT-specific model files
│   ├── README.md          # Documentation of model files
│   ├── SM_TFPT_2loop_v25.yaml
│   ├── E8Cascade_2Loop_Neu_v24fix.yaml
│   ├── E8Cascade_2Loop_Neu_v24fix_PQ.yaml
│   └── ... (more model files)
└── pyrate/                # PyR@TE 3 source code
    ├── pyR@TE.py          # Main entry point
    ├── README.md          # Official PyR@TE documentation
    ├── doc/               # Tutorials and examples
    │   ├── Tutorial.ipynb
    │   └── PyLieDatabase.ipynb
    ├── models/            # Built-in example models (SM, etc.)
    ├── results/           # Output directory for computed RGEs
    └── src/               # Core source code
        ├── Core/          # Beta function calculations
        ├── Definitions/   # Gauge groups, math utilities
        ├── IO/            # Export to LaTeX, Python, Mathematica
        └── PyLie/         # Group theory module
```

## Installation

### Prerequisites

- Python >= 3.6
- pip (Python package manager)

### Setup Virtual Environment

1. **Create the virtual environment** (if not already done):
   ```bash
   cd Pyrate3
   python3 -m venv venv
   ```

2. **Activate and install dependencies**:
   ```bash
   source venv/bin/activate
   pip install -r requirements.txt
   ```

### Required Python Packages

| Package    | Minimum Version |
|------------|-----------------|
| PyYAML     | >= 5.3          |
| SymPy      | >= 1.5          |
| h5py       | >= 2.10         |
| NumPy      | >= 1.18         |
| SciPy      | >= 1.4          |
| Matplotlib | >= 3.1          |

## Usage

### Quick Start

1. **Activate the virtual environment**:
   ```bash
   source activate_env.sh
   ```
   This shows the Python version and installed package versions.

2. **Run a model**:
   ```bash
   ./run_pyrate.sh models/SM_TFPT_2loop_v25.yaml
   ```

### Command Line Options

Run PyR@TE directly with full control:

```bash
source venv/bin/activate
cd pyrate
python pyR@TE.py -m <model_file> [options]
```

Common options:
- `-py` : Generate Python output
- `-tex` : Generate LaTeX output
- `-m` : Generate Mathematica output

### Running a Test

To verify the installation works:

```bash
source venv/bin/activate
python run_test.py
```

This runs a simple Toy U(1) model at 1-loop and prints the beta functions.

## Model Files

### TFPT-Specific Models (in `models/`)

| File | Description |
|------|-------------|
| `SM_TFPT_2loop_v25.yaml` | Clean 2-loop Standard Model baseline for TFPT v2.5 RG fingerprints |
| `E8Cascade_2Loop_Neu_v24fix.yaml` | Bug-fixed E8-cascade mockup for TFPT v2.4 conventions |
| `E8Cascade_2Loop_Neu_v24fix_PQ.yaml` | Same as above with PQ block (n=10) treatment |
| `E8Cascade_TFPT_G8_Enhanced*.yaml` | Exploratory TFPT + G8 bridge variants (3-loop gauge) |

See `models/README.md` for detailed documentation of each model file.

### Built-in Models (in `pyrate/models/`)

Standard PyR@TE example models including the Standard Model (`SM.model`).

## Output

Results are stored in `pyrate/results/` organized by model name:

- **PythonOutput/**: Python code for numerically solving the RGEs
- **LaTeX files**: `.tex` documents with symbolic RGE expressions
- **Mathematica files**: `.m` files for Mathematica integration
- **CSV files**: Numerical coupling evolution data
- **PNG plots**: Visualizations of coupling running

## Key Results Directories

| Directory | Contents |
|-----------|----------|
| `results/E8Cascade2LoopGravity/` | E8 cascade analysis with gravity contributions |
| `results/E8Cascade2LoopGravityV2/` | Updated V2 analysis |
| `results/SM_TFPT_2loop_v25/` | Standard Model baseline results |
| `results/E8_TFPT_G8_theory_compliant/` | G8 theory analysis |

## References

- **PyR@TE 3 Paper**: [arXiv:2007.12700](https://arxiv.org/abs/2007.12700)
- **3-loop Gauge RGEs**: [arXiv:1906.04625](https://arxiv.org/abs/1906.04625)

When using 3-loop gauge results, please cite both papers.

## Troubleshooting

### Virtual Environment Not Found
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Import Errors
Ensure you're in the correct directory and the virtual environment is activated:
```bash
cd Pyrate3
source venv/bin/activate
```

### YAML Parsing Issues
Use `e+` for positive exponents in YAML files (e.g., `1.0e+16` not `1.0e16`) to avoid parsing edge cases.
