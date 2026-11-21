"""Woschni (1967) in-cylinder heat transfer, Heywood eq. 12.19:

    h = C B^-0.2 p^0.8 T^-0.55 w^0.8,  w = C1 Sp_mean + C2 (Vd T_r / (p_r V_r)) (p - p_mot)

with p in kPa, B in m, T in K, w in m/s. C1 = 6.18 (gas exchange), 2.28 (compression, combustion and
expansion); C2 = 0 (gas exchange, compression), 3.24e-3 m/(s K) (combustion and expansion). The reference
state (p_r, V_r, T_r) is IVC; p_mot is the isentropic motored pressure from IVC. All constants come from
spec/thermo.yaml and are packed into a float array for the numba kernels.
"""
from __future__ import annotations

import numpy as np
from numba import njit

(HP_ON, HP_CW, HP_BORE, HP_SP, HP_OMEGA, HP_C1, HP_C2, HP_C1PRE, HP_C2PRE, HP_PREF, HP_VREF, HP_TREF,
 HP_NMOT, HP_VD, HP_AHP, HP_ATWHP, HP_TWL, HP_PUNIT, HP_EB, HP_EP, HP_ET, HP_EW) = range(22)
N_HP = 22


def heat_params(tspec: dict, geom, rpm: float, phase: str, on: bool = True,
                ref: tuple[float, float, float, float] | None = None) -> np.ndarray:
    """Pack Woschni + wall parameters. phase: 'closed' (C1/C2 switch at start of combustion) or 'gas_exchange'.

    ref = (p_r, V_r, T_r, n_motored) at IVC for the combustion term.
    """
    ht = tspec["heat_transfer"]
    hp = np.zeros(N_HP)
    hp[HP_ON] = 1.0 if on else 0.0
    hp[HP_CW] = ht["woschni_coeff"]
    hp[HP_BORE] = geom.bore
    hp[HP_SP] = 2.0 * geom.stroke * rpm / 60.0
    hp[HP_OMEGA] = rpm * 2.0 * np.pi / 60.0
    if phase == "closed":
        hp[HP_C1PRE], hp[HP_C2PRE] = ht["C1_compression"], ht["C2_compression"]
        hp[HP_C1], hp[HP_C2] = ht["C1_combustion"], ht["C2_combustion"]
    elif phase == "gas_exchange":
        hp[HP_C1PRE] = hp[HP_C1] = ht["C1_gas_exchange"]
        hp[HP_C2PRE] = hp[HP_C2] = 0.0
    else:
        raise ValueError(phase)
    p_r, V_r, T_r, n = ref if ref is not None else (1.0, 1.0, 1.0, 1.0)
    hp[HP_PREF], hp[HP_VREF], hp[HP_TREF], hp[HP_NMOT] = p_r, V_r, T_r, n
    hp[HP_VD] = geom.swept_volume
    a_head = ht["head_area_factor"] * geom.piston_area
    a_pist = ht["piston_area_factor"] * geom.piston_area
    hp[HP_AHP] = a_head + a_pist
    hp[HP_ATWHP] = a_head * ht["T_head"] + a_pist * ht["T_piston"]
    hp[HP_TWL] = ht["T_liner"]
    hp[HP_PUNIT] = ht["p_unit"]
    hp[HP_EB], hp[HP_EP], hp[HP_ET], hp[HP_EW] = ht["bore_exp"], ht["p_exp"], ht["T_exp"], ht["w_exp"]
    return hp


@njit(cache=True)
def coefficient(p, T, V, burning, hp):
    """Woschni h [W/m^2/K]; `burning` selects the combustion/expansion constants."""
    if burning:
        c1, c2 = hp[HP_C1], hp[HP_C2]
    else:
        c1, c2 = hp[HP_C1PRE], hp[HP_C2PRE]
    w = c1 * hp[HP_SP]
    if c2 > 0.0:
        p_mot = hp[HP_PREF] * (hp[HP_VREF] / V) ** hp[HP_NMOT]
        dp = p - p_mot
        if dp > 0.0:
            w += c2 * hp[HP_VD] * hp[HP_TREF] / (hp[HP_PREF] * hp[HP_VREF]) * dp
    return (hp[HP_CW] * hp[HP_BORE] ** hp[HP_EB] * (p / hp[HP_PUNIT]) ** hp[HP_EP]
            * T ** hp[HP_ET] * w ** hp[HP_EW])


@njit(cache=True)
def wall_terms(V, hp):
    """Total wetted area and sum(A_i T_wall_i); exposed liner area = pi B (V / A_p) = 4 V / B."""
    a_l = 4.0 * V / hp[HP_BORE]
    return hp[HP_AHP] + a_l, hp[HP_ATWHP] + a_l * hp[HP_TWL]
