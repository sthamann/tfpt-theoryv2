# `Pyrate3/models/`

Project-specific PyR@TE3 model/config files used for TFPT-related RGE studies.

## Notes

- **Numeric YAML scalars**: use `e+` for positive exponents (e.g. `1.0e+16`) to avoid YAML loader edge-cases where `1.0e16` is treated as a string.
- **Model name vs. output folders**: PyR@TE output directories are keyed by the YAML `Name:`. Keep `Name:` aligned with `Pyrate3/pyrate/results/<Name>/...` for reproducibility/fingerprinting.

## Files

- `SM_TFPT_2loop_v25.yaml`
  - **Purpose**: clean **2-loop Standard Model** baseline for the TFPT v2.5 “RG fingerprints”.
  - **Normalization**: keeps `g1` in **SM hypercharge normalization** (`g_Y`). For GUT-normalized quantities use:
    - \(g_{1,\mathrm{GUT}} = \sqrt{5/3}\,g_Y\)
    - \(\alpha_{1,\mathrm{GUT}} = (5/3)\,\alpha_Y\)

- `E8Cascade_2Loop_Neu_v24fix.yaml`
  - **Purpose**: bug-fixed “Neu” 2-loop E8-cascade mockup aligned with TFPT v2.4 conventions:
    - U(1) notes use \(g_{1,\mathrm{GUT}}=\sqrt{5/3}\,g_Y\), \(\alpha_{1,\mathrm{GUT}}=(5/3)\alpha_Y\)
    - \(\lambda_{\mathrm{TFPT}}\) uses \(\gamma(0)=5/6\)

- `E8Cascade_2Loop_Neu_v24fix_PQ.yaml`
  - **Purpose**: same as `E8Cascade_2Loop_Neu_v24fix.yaml`, but treats \(\Phi\) as the **PQ block** (n=10) and anchors \(v_{PQ}\) / \(M_\Phi\) to the TFPT v2.4 table scale (\(\sim 8.86\times 10^{10}\,\mathrm{GeV}\)).

- `E8Cascade 2Loop Gravity.yaml`, `E8Cascade 2Loop Neu.yaml`, `E8_2loop.yaml`
  - **Purpose**: legacy / exploratory E8-cascade + threshold mockups (see inline notes in each file). (`E8Cascade 2Loop Gravity.yaml` is kept aligned with the existing `E8Cascade2LoopGravityV2` PythonOutput.)

- `E8Cascade_TFPT_G8_Enhanced.yaml`, `E8Cascade_TFPT_G8_Enhanced_v2.yaml`
  - **Purpose**: exploratory TFPT + “G8 bridge” variants (loop order 3 for gauge-coupling studies).

- `e8_cascade_1loop.model`, `toy_u1.model`
  - **Purpose**: legacy models used for early tests / scaffolding.

