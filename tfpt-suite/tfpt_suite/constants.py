from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from mpmath import mp


@dataclass(frozen=True)
class TfptConstants:
    """
    Canonical TFPT invariants and frequently used derived quantities (paper v2.5).

    Conventions are taken from `latex/tfpt-theory-fullv25.tex`:
    - c3 := 1/(8π)
    - varphi0 := 1/(6π) + 3/(256π^4)
    - g_{aγγ} := -4 c3 = -1/(2π)
    - β_rad := varphi0/(4π), β_deg := (180/π)*β_rad
    - R^2 scale: M/Mpl := sqrt(8π) * c3^4
    - E8 cascade closure: γ(0) := 5/6, λ := γ(0) / ln(248/60)
    """

    pi: Any
    c3: Any
    varphi0_tree: Any
    delta_top: Any
    varphi0: Any
    delta_star: Any
    b1: Any
    g_a_gamma_gamma: Any
    beta_rad: Any
    beta_deg: Any
    M_over_Mpl: Any
    gamma0: Any
    e8_lambda: Any

    @staticmethod
    def compute() -> "TfptConstants":
        pi = mp.pi

        c3 = mp.mpf(1) / (8 * pi)
        varphi0_tree = mp.mpf(1) / (6 * pi)
        delta_top = mp.mpf(3) / (256 * pi**4)
        varphi0 = varphi0_tree + delta_top
        # Möbius/Z3 flavor anchor (paper v1.06, also referenced as δ⋆): δ⋆ = 3/5 + varphi0/6
        delta_star = mp.mpf(3) / 5 + varphi0 / 6

        # paper convention (GUT normalized hypercharge trace coefficient)
        b1 = mp.mpf(41) / 10

        g_a_gamma_gamma = -4 * c3

        beta_rad = varphi0 / (4 * pi)
        beta_deg = (mp.mpf(180) / pi) * beta_rad

        M_over_Mpl = mp.sqrt(8 * pi) * (c3**4)

        # E8 cascade (paper v2.5 discrete closure)
        gamma0 = mp.mpf(5) / 6
        e8_lambda = gamma0 / mp.log(mp.mpf(248) / mp.mpf(60))

        return TfptConstants(
            pi=pi,
            c3=c3,
            varphi0_tree=varphi0_tree,
            delta_top=delta_top,
            varphi0=varphi0,
            delta_star=delta_star,
            b1=b1,
            g_a_gamma_gamma=g_a_gamma_gamma,
            beta_rad=beta_rad,
            beta_deg=beta_deg,
            M_over_Mpl=M_over_Mpl,
            gamma0=gamma0,
            e8_lambda=e8_lambda,
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "pi": self.pi,
            "c3": self.c3,
            "varphi0_tree": self.varphi0_tree,
            "delta_top": self.delta_top,
            "varphi0": self.varphi0,
            "delta_star": self.delta_star,
            "b1": self.b1,
            "g_a_gamma_gamma": self.g_a_gamma_gamma,
            "beta_rad": self.beta_rad,
            "beta_deg": self.beta_deg,
            "M_over_Mpl": self.M_over_Mpl,
            "gamma0": self.gamma0,
            "e8_lambda": self.e8_lambda,
        }

