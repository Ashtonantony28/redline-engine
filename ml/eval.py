"""`python -m ml.eval`: learned transient fueling (PPO + Simplex supervisor) against the classical
ECU on tip-in / tip-out, and the knock-sensor summary. All numbers are MODELLED (MVEM plant, host
ECU, synthetic accelerometer); there is no physical engine."""
from __future__ import annotations

import json
import sys

import numpy as np

from hil import runner
from ml import ml_spec
from ml.fuel import ppo
from ml.fuel.env import Loop


def run_cycle(cycle: str, rpm: float, policy=None, supervised: bool = True, cfg: dict | None = None, ecu=None) -> dict:
    """One hil/cycles.yaml cycle at `rpm`. policy: None (classical ECU) or obs -> action in [-1, 1]."""
    cfg = cfg or ml_spec()["fuel_ppo"]
    scn = runner.Scenario.from_cycle(cycle, dyno_rpm=[[0.0, rpm]], init_rpm=rpm)
    loop = Loop(scn, cfg, ecu, supervised)
    sup = loop.sup
    trims = []
    while not loop.sess.done:
        loop.advance(0.0 if policy is None else policy(loop.obs()))
        trims.append(loop.trim)
    t0, t1 = scn.windows["transient"]
    # every combustion event counts, fired or not (a fuel-less event has lambda = inf: a violation)
    ev = loop.sess.eng.events
    lam_tr = np.array([e["lambda"] for e in ev if t0 <= e["t"] < t1])
    lam_all = np.array([e["lambda"] for e in ev if e["t"] >= scn.windows["steady"][0]])
    return {"cycle": cycle, "rpm": rpm, "max_dev": float(np.abs(lam_tr - 1.0).max()),
            "rms_dev": float(np.sqrt(np.mean((lam_tr - 1.0) ** 2))),
            "lambda_min": float(lam_all.min()), "lambda_max": float(lam_all.max()),
            "switches": sup.switches, "time_on_classical_s": sup.time_on_classical,
            "max_abs_trim": float(np.max(np.abs(trims)))}


def table(cfg: dict | None = None, actor=None) -> list:
    cfg = cfg or ml_spec()["fuel_ppo"]
    actor = actor or ppo.Actor.load()
    rows = []
    for cycle in cfg["eval_cycles"]:
        for rpm in cfg["eval_rpm"]:
            rows.append({"cycle": cycle, "rpm": float(rpm),
                         "classical": run_cycle(cycle, float(rpm), None, cfg=cfg),
                         "ppo": run_cycle(cycle, float(rpm), actor, True, cfg),
                         "ppo_unsupervised": run_cycle(cycle, float(rpm), actor, False, cfg)})
    return rows


def render(rows: list, cfg: dict) -> str:
    lo, hi = cfg["safe_band"]
    o = ["Transient fueling: lambda excursion per combustion event in the transient window (MODELLED)",
         f"PPO trim authority +/- {100 * float(cfg['trim_limit']):.0f} %, supervisor switch band "
         f"{cfg['switch_band']}, safe band {cfg['safe_band']}", "",
         f"{'cycle':7s} {'rpm':>5s} | {'classical':^17s} | {'PPO + supervisor':^17s} | {'PPO alone':^17s} | "
         f"{'switches':>8s} | {'lambda range (PPO+sup)':>22s} | safe",
         f"{'':7s} {'':>5s} | {'max|dl|':>8s} {'rms':>8s} | {'max|dl|':>8s} {'rms':>8s} | {'max|dl|':>8s} {'rms':>8s} |"]
    for r in rows:
        c, p, u = r["classical"], r["ppo"], r["ppo_unsupervised"]
        safe = lo <= p["lambda_min"] and p["lambda_max"] <= hi
        o.append(f"{r['cycle']:7s} {r['rpm']:5.0f} | {c['max_dev']:8.4f} {c['rms_dev']:8.4f} | {p['max_dev']:8.4f} "
                 f"{p['rms_dev']:8.4f} | {u['max_dev']:8.4f} {u['rms_dev']:8.4f} | {p['switches']:8d} | "
                 f"{p['lambda_min']:10.4f} .. {p['lambda_max']:<8.4f} | {'yes' if safe else 'NO'}")
    mean = lambda who, key: float(np.mean([r[who][key] for r in rows]))
    o += ["", f"mean over the {len(rows)} cycles: max|dl| classical {mean('classical', 'max_dev'):.4f} -> "
              f"PPO + supervisor {mean('ppo', 'max_dev'):.4f};  rms classical {mean('classical', 'rms_dev'):.4f} -> "
              f"PPO + supervisor {mean('ppo', 'rms_dev'):.4f}"]
    return "\n".join(o)


def main(argv=None) -> int:
    cfg = ml_spec()["fuel_ppo"]
    rows = table(cfg)
    print(render(rows, cfg))
    from ml.knock import golden
    m = golden.load_policy()["metrics"]
    print("\nint8 knock sensor (synthetic accelerometer windows, MODELLED), test split:")
    for name, r in m.items():
        print(f"  {name:16s} ROC AUC {r['auc']:.4f}   detection {r['detection']:.3f} "
              f"(spark >= limit + 2 deg: {r['detection_strong']:.3f}) at {100 * r['false_alarm']:.2f} % false alarms")
    print("  bit-exactness on the Vortex-V RTL and cycles per inference: make -C hil cosim")
    return 0


if __name__ == "__main__":
    sys.exit(main())
