"""M9 lockstep co-simulation: `python -m hil.cosim` (or `make -C hil cosim`).

The MVEM plant drives the ECU firmware running on the Vortex-V SoC RTL in Verilator, with the host
reference build of the same sources in lockstep. Reports
  * bit-exactness of the emitted events (tick, type, cyl, arg) and of every slow-loop output struct,
  * cycles per crank-tooth capture and per slow-loop iteration, and the budget at the redline speed,
  * the int8 knock sensor fixture on the RTL against the numpy golden model.
Exit status 0 only if everything is bit-exact. All numbers are MODELLED: there is no physical engine.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

from hil import bridge, runner, rv64
from mvem import params
from redline.spec import ROOT

OUT_PATH = ROOT / "out" / "hil" / "cosim.json"
US_PER_S = 1.0e6


def scenario(name: str) -> runner.Scenario:
    w = dict(rv64.cosim_spec()["windows"][name])
    if w.get("dyno_rpm") == "redline":
        w["dyno_rpm"] = [[0.0, float(params.engine_spec()["engine"]["redline_rpm"])]]
    return runner.Scenario(duration=float(w["duration"]), pedal=w["pedal"], dyno_rpm=w["dyno_rpm"], name=name)


def run_window(name: str, ecu: rv64.Rv64Ecu | None = None, with_glibc: bool = True) -> dict:
    """Run one cosim window in lockstep. Returns counts, differences and per-command cycle arrays."""
    ecu = ecu or rv64.Rv64Ecu()
    shadows = {"host_ref": rv64.RefEcu()}
    if with_glibc:
        shadows["host_glibc"] = bridge.HostEcu()      # the M8 host library (glibc expf / powf)
    ls = rv64.LockstepEcu(ecu, shadows)
    t0 = time.time()
    c0 = ecu.total_cycles()
    res = runner.run(scenario(name), ls)
    cyc = {k: np.array(v, dtype=np.int64) for k, v in ecu.cycles.items()}
    return {"window": name, "events": ls.n_events, "slow_steps": ls.n_slow, "calls": len(ls.inputs),
            "diff": ls.diff, "cycles": cyc, "core_cycles": ecu.total_cycles() - c0, "wall_s": time.time() - t0,
            "status": res.status, "result": res, "rpm": float(res.log["rpm"][-1])}


def tooth_period_s(rpm: float) -> float:
    teeth = int(params.mvem_spec()["trigger"]["teeth_total"])
    return params.SEC_PER_MIN / (rpm * teeth)


def budget(runs: list) -> dict:
    """Cycle budget at the redline speed for the assumed core clock (hil/cosim.yaml)."""
    cfg = rv64.cosim_spec()
    f_clk = float(cfg["core_clock_hz"])
    redline = float(params.engine_spec()["engine"]["redline_rpm"])
    slow_hz = float(params.ecu_spec()["hardware"]["slow_loop_hz"])
    crank = np.concatenate([r["cycles"]["crank"] for r in runs])
    slow = np.concatenate([r["cycles"]["slow"] for r in runs])
    red = [r for r in runs if r["window"] == "redline"]
    crank_red = np.concatenate([r["cycles"]["crank"] for r in red]) if red else crank
    t_tooth = tooth_period_s(redline)
    tooth_cycles = t_tooth * f_clk
    real_teeth = int(params.mvem_spec()["trigger"]["teeth_total"]) - int(params.mvem_spec()["trigger"]["teeth_missing"])
    tooth_rate = redline / params.SEC_PER_MIN * real_teeth
    util = (slow[:, 1].mean() * slow_hz + crank_red[:, 1].mean() * tooth_rate) / f_clk
    worst_tooth = int(crank[:, 1].max())
    worst_slow = int(slow[:, 1].max())
    return {
        "core_clock_hz": f_clk, "redline_rpm": redline,
        "tooth_period_us": t_tooth * US_PER_S, "tooth_period_cycles": tooth_cycles,
        "crank_handler_cycles_max": int(crank[:, 0].max()), "crank_handler_cycles_mean": float(crank[:, 0].mean()),
        "crank_post_to_ack_cycles_max": worst_tooth,
        "crank_redline_post_to_ack_cycles_max": int(crank_red[:, 1].max()),
        "crank_worst_us": worst_tooth / f_clk * US_PER_S,
        "crank_fraction_of_tooth_period": worst_tooth / tooth_cycles,
        "crank_within_margin": bool(worst_tooth <= float(cfg["budget_margin"]) * tooth_cycles),
        "slow_handler_cycles_max": int(slow[:, 0].max()), "slow_handler_cycles_mean": float(slow[:, 0].mean()),
        "slow_post_to_ack_cycles_max": worst_slow,
        "slow_worst_us": worst_slow / f_clk * US_PER_S, "slow_period_us": US_PER_S / slow_hz,
        "slow_fraction_of_period": worst_slow * slow_hz / f_clk,
        "cpu_utilisation_at_redline": float(util),
        # no interrupts on Vortex-V yet: a tooth that arrives during a slow-loop iteration waits for it
        "polled_worst_tooth_latency_us": (worst_slow + worst_tooth) / f_clk * US_PER_S,
        "polled_latency_within_tooth_period": bool(worst_slow + worst_tooth <= tooth_cycles),
    }


def knock_on_rtl(ecu: rv64.Rv64Ecu) -> dict:
    """Every committed fixture vector through the Xmlacc kernel on the RTL vs the golden logit."""
    from ml.knock import golden, train
    fix = json.loads(train.FIXTURE_PATH.read_text())
    policy = golden.load_policy()
    bad, cycles = 0, []
    for v in fix["vectors"]:
        logit, cyc = ecu.knock_infer(v["x_q"])
        gold = int(golden.forward_int8(policy, np.array(v["x_q"], dtype=np.int8))[0][0])
        bad += int(logit != v["logit"] or gold != v["logit"])
        cycles.append(cyc)
    return {"vectors": len(fix["vectors"]), "mismatches": bad, "cycles_max": int(max(cycles)),
            "cycles_mean": float(np.mean(cycles)), "metrics": policy["metrics"]}


def report(runs: list, bud: dict, kn: dict | None) -> str:
    o = ["Redline M9 co-simulation: MVEM plant <-> ECU firmware on Vortex-V RTL (Verilator). All values MODELLED.",
         f"core: vortex_soc, {rv64.VORTEX_V_ROOT}, assumed clock {bud['core_clock_hz'] / 1e6:.0f} MHz", ""]
    for r in runs:
        o.append(f"window {r['window']:8s} {r['rpm']:5.0f} rpm: {r['calls']} ECU calls, {r['events']} events, "
                 f"{r['slow_steps']} slow-loop steps, {r['core_cycles']} core cycles, {r['wall_s']:.1f} s wall")
        for name, d in r["diff"].items():
            ok = d["events"] == 0 and d["slow"] == 0
            o.append(f"    vs {name:10s}: event differences {d['events']}, slow-loop output differences {d['slow']}"
                     f"  -> {'BIT-EXACT' if ok else 'DIFFERS'}")
    o += ["", f"cycle budget at {bud['redline_rpm']:.0f} rpm (tooth period {bud['tooth_period_us']:.1f} us = "
              f"{bud['tooth_period_cycles']:.0f} cycles):",
          f"  crank-tooth capture: handler max {bud['crank_handler_cycles_max']} / mean "
          f"{bud['crank_handler_cycles_mean']:.0f} cycles; post-to-acknowledge max "
          f"{bud['crank_post_to_ack_cycles_max']} cycles = {bud['crank_worst_us']:.1f} us = "
          f"{100 * bud['crank_fraction_of_tooth_period']:.0f} % of the tooth period",
          f"  slow loop (soft-float): handler max {bud['slow_handler_cycles_max']} / mean "
          f"{bud['slow_handler_cycles_mean']:.0f} cycles = {bud['slow_worst_us']:.0f} us of the "
          f"{bud['slow_period_us']:.0f} us period ({100 * bud['slow_fraction_of_period']:.0f} %)",
          f"  CPU utilisation at redline (slow loop + teeth): {100 * bud['cpu_utilisation_at_redline']:.0f} %",
          f"  polled main loop, tooth arriving during a slow-loop step: latency up to "
          f"{bud['polled_worst_tooth_latency_us']:.0f} us "
          f"({'within' if bud['polled_latency_within_tooth_period'] else 'EXCEEDS'} the tooth period)"]
    if kn:
        m = kn["metrics"]["int8"]
        o += ["", f"int8 knock sensor on Xmlacc: {kn['vectors']} fixture vectors, {kn['mismatches']} mismatches vs golden"
                  f"  -> {'BIT-EXACT' if kn['mismatches'] == 0 else 'DIFFERS'}; {kn['cycles_max']} cycles per inference "
                  f"(max) = {kn['cycles_max'] / bud['core_clock_hz'] * US_PER_S:.0f} us",
              f"  accuracy (synthetic test split): AUC {m['auc']:.4f}, detection {m['detection']:.3f} at "
              f"{100 * m['false_alarm']:.2f} % false alarms"]
    return "\n".join(o)


def bit_exact(runs: list, shadow: str = "host_ref") -> bool:
    return all(r["diff"][shadow]["events"] == 0 and r["diff"][shadow]["slow"] == 0 for r in runs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--windows", nargs="+", default=["tipin", "redline"],
                    choices=list(rv64.cosim_spec()["windows"]))
    a = ap.parse_args(argv)
    missing = rv64.toolchain_missing()
    if missing:
        print(f"missing tools: {missing} (run $VORTEX_V_ROOT/scripts/setup-toolchain.sh)", file=sys.stderr)
        return 2
    ecu = rv64.Rv64Ecu()
    runs = [run_window(w, ecu) for w in a.windows]
    bud = budget(runs)
    kn = knock_on_rtl(ecu)
    soc_pass = ecu.exit()
    print(report(runs, bud, kn))
    print(f"firmware exit through tohost: {'PASS' if soc_pass else 'FAIL'}")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    slim = [{k: v for k, v in r.items() if k not in ("cycles", "result")} for r in runs]
    OUT_PATH.write_text(json.dumps({"runs": slim, "budget": bud, "knock": kn}, indent=2, default=str) + "\n")
    print(f"wrote {OUT_PATH}")
    ok = bit_exact(runs) and kn["mismatches"] == 0 and soc_pass
    print("RESULT: " + ("PASS (bit-exact vs the host reference build)" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
