"""0D component elements applied to a gasdyn.network.Network between kernel calls (M10 turbo stretch).

The numba kernel (gasdyn.fv) knows pipes, orifices and three volume kinds. Turbomachinery is added on top
of it, without touching the kernel, as operator-split 0D elements that act on the volume states every
sub-step (gasdyn.engine `substep_fn`):

* MapCompressor: map-based compressor between a fixed inlet state (ambient minus the air-filter loss) and
  a charge-air plenum. Each sub-step its mass flow at the current shaft speed and plenum pressure is
  taken from the map (quasi-steady: the compressor has no volume of its own) and injected into the plenum
  as mass plus stagnation enthalpy at the intercooler outlet temperature.
* Intercooler: effectiveness (T_out = T_in - eps (T_in - T_coolant)) plus a quadratic, density-corrected
  pressure drop; it sits between the compressor outlet and the charge-air plenum.
* MapTurbine: the turbine-housing pipe ends in an orifice of the turbine's nozzle-equivalent area, so the
  kernel's isentropic orifice boundary *is* the turbine flow map (Watson & Janota, "Turbocharging the
  Internal Combustion Engine", Macmillan 1982, ch. 4). The isentropic work times the blade-speed-ratio
  efficiency is withdrawn from the post-turbine volume, so the gas leaving the turbine is cooled by
  exactly the shaft work (energy is conserved to the time-integration error of the split).
* Wastegate: the orifice area at the end of a pipe in parallel with the turbine.

Maps are duck-typed objects (turbo.maps supplies synthetic ones). All numbers are MODELLED (no physical
engine, CLAUDE.md rule 5).
"""
from __future__ import annotations

from dataclasses import dataclass

from gasdyn import fv

RIGHT_END = 1


def gas_cp(gamma: float, R: float) -> float:
    return gamma * R / (gamma - 1.0)


def plenum_inject(net, iv: int, dm: float, T0_in: float) -> None:
    """Add dm [kg] of gas with stagnation temperature T0_in to a fixed-volume plenum (fresh charge, Y = 0).
    Energy convention of gasdyn.fv: plenum E = m cv T, inflow carries cp T0."""
    g, R = net.v_gam[iv], net.v_R[iv]
    net.S[iv, fv.S_M] += dm
    net.S[iv, fv.S_E] += dm * gas_cp(g, R) * T0_in


def plenum_remove_energy(net, iv: int, dE: float) -> None:
    net.S[iv, fv.S_E] -= dE


def plenum_pT(net, iv: int):
    return net.plenum_state(iv)


def pipe_end_state(net, ip: int, side: int = RIGHT_END):
    """Mass flow towards that end [kg/s] and the static / stagnation state of the end cell.
    Returns (mdot, p, T, p0, T0)."""
    g, R, A = net.p_gam[ip], net.p_R[ip], net.p_area[ip]
    i = int(net.p_off[ip + 1] - 1) if side == RIGHT_END else int(net.p_off[ip])
    rho = net.U[i, 0]
    u = net.U[i, 1] / rho
    p = (g - 1.0) * (net.U[i, 2] - 0.5 * rho * u * u)
    T = p / (rho * R)
    cp = gas_cp(g, R)
    T0 = T + 0.5 * u * u / cp
    p0 = p * (T0 / T) ** (g / (g - 1.0))
    mdot = rho * u * A * (1.0 if side == RIGHT_END else -1.0)
    return mdot, p, T, p0, T0


@dataclass
class Intercooler:
    effectiveness: float
    T_coolant: float
    dp_ref: float          # Pa at m_ref and rho_ref
    m_ref: float           # kg/s
    rho_ref: float         # kg/m^3

    def T_out(self, T_in: float) -> float:
        return T_in - self.effectiveness * (T_in - self.T_coolant)

    def dp(self, mdot: float, rho: float) -> float:
        """Quadratic loss at constant loss coefficient: dp = dp_ref (m / m_ref)^2 (rho_ref / rho)."""
        return self.dp_ref * (mdot / self.m_ref) ** 2 * (self.rho_ref / rho)


class MapCompressor:
    """Compressor + air filter + intercooler feeding the charge-air plenum `iv`.

    cmap must provide `operating_point(omega, p01, T01, demand)` (see turbo.maps.CompressorMap), where
    demand(m) is the compressor total-pressure ratio the downstream system requires at flow m."""

    def __init__(self, cmap, iv: int, p_amb: float, T_amb: float, filter_dp, intercooler: Intercooler):
        self.cmap, self.iv = cmap, iv
        self.p_amb, self.T_amb = p_amb, T_amb
        self.filter_dp = filter_dp          # callable(mdot, rho) -> Pa
        self.ic = intercooler

    def inlet(self, mdot: float):
        rho = self.p_amb / (self.cmap.R * self.T_amb)
        return self.p_amb - self.filter_dp(mdot, rho), self.T_amb

    def step(self, net, omega: float, dt: float) -> dict:
        p_cac, T_cac = plenum_pT(net, self.iv)
        R = self.cmap.R
        rho_ic_ref = p_cac / (R * T_cac)

        def demand(m):
            p1, _ = self.inlet(m)
            return (p_cac + self.ic.dp(m, rho_ic_ref)) / p1

        p1, T1 = self.inlet(0.0)
        op = self.cmap.operating_point(omega, p1, T1, demand)
        m = op["mdot"]
        p1, T1 = self.inlet(m)
        op = dict(op)
        op.update(p1=p1, T1=T1, p2=op["pr"] * p1, p_charge=p_cac, T_charge=T_cac)
        T_ic = self.ic.T_out(op["T2"])
        op["T_ic_out"] = T_ic
        op["q_ic"] = m * gas_cp(self.cmap.gamma, R) * (op["T2"] - T_ic)
        if m > 0.0:
            plenum_inject(net, self.iv, m * dt, T_ic)
        return op


class MapTurbine:
    """Turbine at the right end of pipe `ip`, discharging into plenum `iv_out`.

    tmap must provide `power(mdot, p01, T01, p_out, omega)` -> dict with 'power' [W] (turbo.maps)."""

    def __init__(self, tmap, ip: int, iv_out: int):
        self.tmap, self.ip, self.iv_out = tmap, ip, iv_out

    def step(self, net, omega: float, dt: float) -> dict:
        m, p, T, p0, T0 = pipe_end_state(net, self.ip)
        p4, T4 = plenum_pT(net, self.iv_out)
        r = self.tmap.power(max(m, 0.0), p0, T0, p4, omega)
        P = r["power"]
        if P > 0.0:
            plenum_remove_energy(net, self.iv_out, P * dt)
        r = dict(r)
        r.update(mdot=m, p03=p0, T03=T0, p4=p4, T4=T4)
        return r


class Wastegate:
    """Variable orifice at the end of the bypass pipe: area = opening (0..1) x max area."""

    def __init__(self, slot: int, max_area: float):
        self.slot, self.max_area = slot, max_area

    def set(self, net, opening: float) -> None:
        net.sl_const[self.slot] = self.max_area * min(max(opening, 0.0), 1.0)

    def flow(self, net, ip: int) -> float:
        return pipe_end_state(net, ip)[0]


__all__ = ["Intercooler", "MapCompressor", "MapTurbine", "Wastegate", "plenum_inject", "plenum_remove_energy",
           "pipe_end_state", "gas_cp"]
