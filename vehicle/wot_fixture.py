"""Write the committed WOT-envelope fixture fixtures/wot/wot_coarse.csv from the M4 model (sim.wot).

    python -m vehicle.wot_fixture

out/ is not committed, so the vehicle model reads its full-load envelope from this small fixture
instead of out/sim/wot.csv. Regenerating it is a reviewed change (CLAUDE.md). The sweep comes from
spec/vehicle.yaml simulation.wot_fixture_sweep_rpm; each point is cold-started, so the result does
not depend on the number of worker processes. All numbers are MODELLED.
"""
from __future__ import annotations

import argparse
import os

from redline.spec import ROOT
from sim.wot import peaks, run_sweep, sweep_speeds, write_csv
from vehicle.params import vehicle_spec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args(argv)
    sim = vehicle_spec()["simulation"]
    res = run_sweep(sweep_speeds(*sim["wot_fixture_sweep_rpm"]), None, a.jobs)
    path = ROOT / sim["wot_fixture"]
    write_csv(res, path)
    pk = peaks(res)
    print(f"wrote {path} ({len(res)} points, MODELLED): peak {pk['torque_nm']:.1f} Nm at {pk['torque_rpm']:.0f} rpm, "
          f"{pk['power_kw']:.1f} kW at {pk['power_rpm']:.0f} rpm, all converged {pk['all_converged']}")


if __name__ == "__main__":
    main()
