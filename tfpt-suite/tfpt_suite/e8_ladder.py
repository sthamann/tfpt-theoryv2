from __future__ import annotations

from mpmath import mp

from tfpt_suite.constants import TfptConstants


def e8_phi_n(*, n: int, c: TfptConstants) -> mp.mpf:
    """
    Log-exact E8 ladder (paper v2.4/v2.5):

      varphi_n = varphi0 * exp(-gamma(0)) * (D_n/D_1)^lambda   for n>=1
      D_n = 60 - 2n, D_1 = 58

    Note:
    - This returns dimensionless ladder steps. Dimensionful scales require a block
      calibration factor ζ_B and the (unreduced) Planck mass.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return mp.mpf(c.varphi0)
    D1 = mp.mpf(58)
    Dn = mp.mpf(60 - 2 * n)
    if Dn <= 0:
        raise ValueError(f"Invalid D_n for n={n}: Dn={Dn}")
    return mp.mpf(c.varphi0) * mp.e ** (-mp.mpf(c.gamma0)) * (Dn / D1) ** mp.mpf(c.e8_lambda)
