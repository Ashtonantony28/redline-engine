"""Minimal PPO (clip) for the fuel-trim environment: torch, CPU, one thread, seeded.

`python -m ml.fuel.ppo` trains and writes ml/policies/fuel_ppo.json: the deterministic actor
(tanh MLP -> tanh mean) as plain lists, so evaluation needs numpy only."""
from __future__ import annotations

import json
import time

import numpy as np

from ml import POLICY_DIR, ml_spec
from ml.fuel.env import OBS_NAMES, FuelEnv

POLICY_PATH = POLICY_DIR / "fuel_ppo.json"


class Actor:
    """numpy forward pass of the saved policy mean: action in [-1, 1] (x trim_limit = fuel trim)."""

    def __init__(self, policy: dict):
        self.layers = [(np.asarray(L["W"]), np.asarray(L["b"])) for L in policy["actor"]]

    @classmethod
    def load(cls, path=POLICY_PATH) -> "Actor":
        return cls(json.loads(path.read_text()))

    def __call__(self, obs: np.ndarray) -> float:
        a = np.asarray(obs, dtype=np.float64)
        for W, b in self.layers:
            a = np.tanh(W @ a + b)
        return float(a[0])


def train(cfg: dict | None = None, log=print) -> dict:
    import torch
    from torch import nn
    cfg = cfg or ml_spec()["fuel_ppo"]
    seed = int(cfg["seed"])
    torch.manual_seed(seed)
    torch.set_num_threads(1)
    env = FuelEnv(cfg)
    n_obs = len(OBS_NAMES)

    def mlp(n_out):
        dims = [n_obs, *[int(h) for h in cfg["hidden"]]]
        mods = []
        for i in range(len(dims) - 1):
            mods += [nn.Linear(dims[i], dims[i + 1]), nn.Tanh()]
        return nn.Sequential(*mods, nn.Linear(dims[-1], n_out))

    actor, critic = mlp(1), mlp(1)
    with torch.no_grad():                     # start at trim ~ 0: the classical controller
        actor[-1].weight.mul_(0.01)
        actor[-1].bias.zero_()
    log_std = nn.Parameter(torch.full((1,), float(cfg["init_log_std"])))
    opt = torch.optim.Adam([*actor.parameters(), *critic.parameters(), log_std], lr=float(cfg["learning_rate"]))
    T, gamma, lam = int(cfg["rollout_steps"]), float(cfg["gamma"]), float(cfg["gae_lambda"])
    obs, _ = env.reset(seed=seed)
    ep_ret, returns, history, t0, steps = 0.0, [], [], time.time(), 0
    while steps < int(cfg["total_steps"]):
        O = np.zeros((T, n_obs), np.float32)
        A, LP, R, V, D = (np.zeros(T, np.float32) for _ in range(5))
        for i in range(T):
            with torch.no_grad():
                o = torch.as_tensor(obs)
                dist = torch.distributions.Normal(torch.tanh(actor(o)), log_std.exp())
                a = dist.sample()
                O[i], A[i], LP[i], V[i] = obs, a.item(), dist.log_prob(a).item(), critic(o).item()
            obs, r, _, trunc, _ = env.step(a.numpy())
            R[i], ep_ret = r, ep_ret + r
            if trunc:                         # time limit: bootstrap from the value of the last state
                with torch.no_grad():
                    R[i] += gamma * critic(torch.as_tensor(obs)).item()
                D[i] = 1.0
                returns.append(ep_ret)
                ep_ret = 0.0
                obs, _ = env.reset()
        steps += T
        with torch.no_grad():
            last_v = critic(torch.as_tensor(obs)).item()
        adv = np.zeros(T, np.float32)
        g = 0.0
        for i in reversed(range(T)):
            nv = last_v if i == T - 1 else V[i + 1]
            nonterm = 1.0 - D[i]
            delta = R[i] + gamma * nv * nonterm - V[i]
            g = delta + gamma * lam * nonterm * g
            adv[i] = g
        ret = adv + V
        tO, tA, tLP = torch.as_tensor(O), torch.as_tensor(A), torch.as_tensor(LP)
        tAdv = torch.as_tensor((adv - adv.mean()) / (adv.std() + 1e-8))
        tRet = torch.as_tensor(ret)
        gen = torch.Generator().manual_seed(seed + steps)
        for _ in range(int(cfg["epochs"])):
            perm = torch.randperm(T, generator=gen)
            for j in range(0, T, int(cfg["minibatch"])):
                idx = perm[j:j + int(cfg["minibatch"])]
                dist = torch.distributions.Normal(torch.tanh(actor(tO[idx]))[:, 0], log_std.exp())
                ratio = (dist.log_prob(tA[idx]) - tLP[idx]).exp()
                clip = float(cfg["clip"])
                pol = -torch.min(ratio * tAdv[idx], ratio.clamp(1.0 - clip, 1.0 + clip) * tAdv[idx]).mean()
                val = (critic(tO[idx])[:, 0] - tRet[idx]).pow(2).mean()
                loss = pol + float(cfg["value_coef"]) * val - float(cfg["entropy_coef"]) * dist.entropy().mean()
                opt.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_([*actor.parameters(), *critic.parameters()], float(cfg["max_grad_norm"]))
                opt.step()
        recent = float(np.mean(returns[-10:])) if returns else float("nan")
        history.append({"steps": steps, "mean_return_last10": recent})
        log(f"  steps {steps:7d}  episodes {len(returns):4d}  mean return (last 10) {recent:8.4f}  "
            f"std {log_std.exp().item():.3f}  {time.time() - t0:5.0f} s")
    lin = [m for m in actor if isinstance(m, nn.Linear)]
    return {"actor": [{"W": m.weight.detach().numpy().astype(float).round(6).tolist(),
                       "b": m.bias.detach().numpy().astype(float).round(6).tolist()} for m in lin],
            "obs": OBS_NAMES, "trim_limit": float(cfg["trim_limit"]),
            "source": {"generator": "python -m ml.fuel.ppo", "spec": "spec/ml.yaml fuel_ppo", "seed": seed,
                       "total_steps": steps, "train_wall_s": round(time.time() - t0, 1),
                       "note": "MODELLED: trained on the MVEM plant + host ECU"},
            "history": history}


def main(argv=None) -> int:
    print("PPO transient fueling (MODELLED plant):")
    policy = train()
    POLICY_PATH.parent.mkdir(parents=True, exist_ok=True)
    POLICY_PATH.write_text(json.dumps(policy, separators=(",", ":")) + "\n")
    print(f"wrote {POLICY_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
