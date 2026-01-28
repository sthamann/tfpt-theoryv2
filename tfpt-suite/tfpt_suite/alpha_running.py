from __future__ import annotations

from dataclasses import dataclass

from mpmath import mp


@dataclass(frozen=True)
class AlphaRunningInputs:
    """
    Minimal inputs to map alpha(0) to an effective alpha(MZ) reference quantity.

    Notes:
    - This is a *comparison layer* (SM/QED running), not a TFPT prediction.
    - We implement the leading 1-loop leptonic vacuum polarization piece and take
      Δα_had^(5)(MZ) as an external PDG input.
    - Electroweak decoupling (W/Z scheme shift, top decoupling) is represented
      explicitly as finite pieces, not as a hidden fit parameter.
    """

    MZ_GeV: mp.mpf
    m_e_GeV: mp.mpf
    m_mu_GeV: mp.mpf
    m_tau_GeV: mp.mpf
    delta_alpha_had5_MZ: mp.mpf
    delta_alpha_msbar_on_shell_shift_MZ: mp.mpf = mp.mpf(0)
    delta_alpha_top_decoupling_MZ: mp.mpf = mp.mpf(0)
    delta_alpha_extra_msbar_MZ: mp.mpf = mp.mpf(0)


def alpha_running_inputs_from_pdg(*, pdg: dict) -> AlphaRunningInputs:
    """
    Construct AlphaRunningInputs from a parsed alpha_running_pdg.json payload.
    """
    if not isinstance(pdg, dict):
        raise TypeError("pdg must be a dict")
    try:
        lepton_masses = pdg["lepton_masses_GeV"]
        delta_had = pdg["delta_alpha_had5_MZ"]
    except KeyError as exc:
        raise KeyError(f"alpha_running_pdg.json missing required key: {exc}") from exc

    return AlphaRunningInputs(
        MZ_GeV=mp.mpf(str(pdg["MZ_GeV"])),
        m_e_GeV=mp.mpf(str(lepton_masses["m_e"]["value"])),
        m_mu_GeV=mp.mpf(str(lepton_masses["m_mu"]["value"])),
        m_tau_GeV=mp.mpf(str(lepton_masses["m_tau"]["value"])),
        delta_alpha_had5_MZ=mp.mpf(str(delta_had["mean"])),
        delta_alpha_msbar_on_shell_shift_MZ=mp.mpf(str(pdg.get("delta_alpha_msbar_on_shell_shift_MZ", {}).get("mean", 0.0))),
        delta_alpha_top_decoupling_MZ=mp.mpf(str(pdg.get("delta_alpha_top_decoupling_MZ", {}).get("mean", 0.0))),
        delta_alpha_extra_msbar_MZ=mp.mpf(str(pdg.get("delta_alpha_extra_msbar_MZ", {}).get("mean", 0.0))),
    )


def delta_alpha_lept_1loop(*, alpha0: mp.mpf, MZ_GeV: mp.mpf, m_lepton_GeV: mp.mpf) -> mp.mpf:
    """
    1-loop leptonic vacuum polarization contribution:

      Δα_l(MZ^2) = (α(0) / (3π)) * ( ln(MZ^2 / m_l^2) - 5/3 )

    This is sufficient for an assumption-explicit secondary consistency check
    against ᾱ^(5)(MZ). For publication-grade work, upgrade to the full SM
    electroweak scheme conversion as needed.
    """
    if MZ_GeV <= 0 or m_lepton_GeV <= 0:
        raise ValueError("MZ_GeV and m_lepton_GeV must be positive")
    return (alpha0 / (mp.mpf(3) * mp.pi)) * (mp.log((MZ_GeV**2) / (m_lepton_GeV**2)) - mp.mpf(5) / mp.mpf(3))


def delta_alpha_lept_total_1loop(*, alpha0: mp.mpf, inputs: AlphaRunningInputs) -> mp.mpf:
    return mp.fsum(
        [
            delta_alpha_lept_1loop(alpha0=alpha0, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_e_GeV),
            delta_alpha_lept_1loop(alpha0=alpha0, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_mu_GeV),
            delta_alpha_lept_1loop(alpha0=alpha0, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_tau_GeV),
        ]
    )


def delta_alpha_extra_total(*, inputs: AlphaRunningInputs) -> mp.mpf:
    """
    Aggregate explicit finite pieces beyond leptonic + hadronic VP.
    """
    return mp.fsum(
        [
            mp.mpf(inputs.delta_alpha_msbar_on_shell_shift_MZ),
            mp.mpf(inputs.delta_alpha_top_decoupling_MZ),
            mp.mpf(inputs.delta_alpha_extra_msbar_MZ),
        ]
    )


def alpha_bar5_MZ_from_alpha0(*, alpha0: mp.mpf, inputs: AlphaRunningInputs) -> tuple[mp.mpf, dict[str, mp.mpf]]:
    """
    Map alpha(0) to an effective ᾱ^(5)(MZ) via
      α(MZ) = α(0) / (1 - Δα(MZ))
    with
      Δα(MZ) = Δα_lept(MZ) + Δα_had^(5)(MZ)
    """
    if alpha0 <= 0:
        raise ValueError("alpha0 must be positive")

    dal_lept = delta_alpha_lept_total_1loop(alpha0=alpha0, inputs=inputs)
    dal_had5 = inputs.delta_alpha_had5_MZ
    dal_msbar_shift = inputs.delta_alpha_msbar_on_shell_shift_MZ
    dal_top = inputs.delta_alpha_top_decoupling_MZ
    dal_extra = inputs.delta_alpha_extra_msbar_MZ
    dal_extra_total = delta_alpha_extra_total(inputs=inputs)
    dal_tot = dal_lept + dal_had5 + dal_extra_total
    denom = mp.mpf(1) - dal_tot
    if denom == 0:
        raise ZeroDivisionError("1 - Δα(MZ) = 0")
    alpha_mz = alpha0 / denom
    return alpha_mz, {
        "delta_alpha_lept_1loop": dal_lept,
        "delta_alpha_had5": dal_had5,
        "delta_alpha_msbar_on_shell_shift": dal_msbar_shift,
        "delta_alpha_top_decoupling": dal_top,
        "delta_alpha_extra_msbar": dal_extra,
        "delta_alpha_extra_total": dal_extra_total,
        "delta_alpha_total": dal_tot,
    }

