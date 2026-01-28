# Recent Fixes Applied

## 1. Fixed LaTeX Rendering Issue in Plots

**Problem**: The gauge_running.png was empty/white (15KB)

**Cause**: PyR@TE generated `lambda_` coupling without proper LaTeX string, causing matplotlib to fail when rendering "lambda_" as LaTeX.

**Solution**:
- Added `latex='\\lambda'` parameter to lambda_ Coupling definition in `PythonOutput/E8Cascade2LoopGravity.py`
- Replaced PyR@TE's built-in plot function with custom plotting in `solver.py`

**Result**: Plots now generate correctly (252KB vs 15KB)

## 2. Added Output Logging to File

**Problem**: Console output was not saved to file

**Solution**: 
- Created TeeOutput class in `run_e8cascade.py` that writes to both console and file
- All output now saved to `results/results.txt`
- Also captures stderr for error messages

**Result**: Complete run log available for analysis

## Files Modified:
1. `PythonOutput/E8Cascade2LoopGravity.py` - Fixed lambda_ LaTeX
2. `src/e8cascade/solver.py` - Replaced plot function
3. `run_e8cascade.py` - Added output logging

## Usage:
```bash
# Standard run - output saved to results/results.txt
python run_e8cascade.py

# With gravity portal
python run_e8cascade.py --gravity --cR-factor 100 --output results_gravity
```

All console output is now automatically saved to `results.txt` in the output directory! 