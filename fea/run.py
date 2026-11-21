"""M6 runner: `python -m fea.run --all` (PLAN.md sec. 5).

Runs the analytic anchors, then for each part (conrod, crankshaft, piston, main_cap) the design-iteration
loop and the three-level mesh convergence study; writes out/fea/report.md, out/fea/results.json,
out/fea/iterations.json and the stress/temperature PNGs. `--apply` writes parameters changed by the loop
back to spec/cad.yaml and appends them to fea/design_history.json (a reviewed change). Exit code 0 only
if every part meets its SF targets and the convergence criterion. All results are MODELLED.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from fea.anchors import run_all
from fea.config import base_params, fea_spec, spec_hash
from fea.report import write_report, write_results
from fea.study import PARTS, apply_to_spec, dump_iterations, part_module, run_part
from redline.spec import ROOT

HISTORY_PATH = Path(__file__).resolve().parent / "design_history.json"


def volume_check(name, p, fspec) -> str | None:
    from cad.conrod import conrod
    from cad.piston import piston
    from fea.geometry import conrod_solid, piston_solid
    if name == "conrod":
        a, b = conrod_solid(p, fspec["conrod"], quarter=False).volume, conrod(p).volume
        return f"FEA rod {a / 1e3:.2f} cm^3 vs CAD {b / 1e3:.2f} cm^3 ({a / b - 1:+.1%}; blends added)"
    if name == "piston":
        a, b = piston_solid(p, fspec["piston"], quarter=False).volume, piston(p).volume
        return (f"FEA piston {a / 1e3:.2f} cm^3 vs CAD {b / 1e3:.2f} cm^3 ({a / b - 1:+.1%}; ring grooves filled, "
                f"circlip grooves removed, panel / skirt / relief blends added)")
    return None


def run_probes(name, p, fspec) -> list:
    from fea.config import with_overrides
    from fea.study import ITERATE_LEVEL, part_loads, part_mass
    mod = part_module(name)
    out = []
    for i, ov in enumerate(fspec[mod.SECTION].get("probes", [])):
        q = with_overrides(p, mod.SECTION, ov)
        r = mod.analyse(q, fspec, ITERATE_LEVEL, part_loads(name, q, fspec), tag=f"_probe{i}")
        mass = part_mass(q, "piston_assembly" if name == "piston" else "conrod")
        out.append({"params": ov, "sf": r.sf, "peak_Pa": r.peak, "mass_kg": mass})
        print(f"  [{name}] probe {ov}: SF {r.sf}, mass {mass:.4f} kg")
    return out


def load_history() -> list:
    return json.loads(HISTORY_PATH.read_text()) if HISTORY_PATH.exists() else []


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--all", action="store_true", help="run every part")
    ap.add_argument("--parts", default="", help="comma-separated subset of " + ",".join(PARTS))
    ap.add_argument("--apply", action="store_true", help="write loop-changed parameters to spec/cad.yaml")
    ap.add_argument("--no-anchors", action="store_true")
    a = ap.parse_args(argv)
    names = list(PARTS) if a.all or not a.parts else [s.strip() for s in a.parts.split(",") if s.strip()]
    fspec = fea_spec()
    out = ROOT / fspec["solver"]["outdir"]
    out.mkdir(parents=True, exist_ok=True)
    p = base_params()
    print("Redline M6 FEA (MODELLED results, no physical engine)")
    anchors = None
    if not a.no_anchors:
        anchors = run_all(ROOT / fspec["solver"]["workdir"] / "anchors")
        print("anchors:", ", ".join(f"{k} ok" for k in anchors))
    studies = []
    for name in names:
        print(f"[{name}]")
        st = run_part(name, p, fspec)
        mod = part_module(name)
        st.figures = mod.plots(st.finest, out, p) if st.finest is not None and st.finest.mesh is not None else []
        st.volume_check = volume_check(name, p, fspec)
        st.probes = run_probes(name, p, fspec) if not st.passed else []
        print(f"[{name}] SF " + ", ".join(f"{k} {v:.2f}" for k, v in st.sf.items())
              + f"; convergence {st.convergence.get('change_finest_two', float('nan')):.2%}; "
              + ("PASS" if st.passed else "FAIL") + f" ({st.seconds:.0f} s)")
        studies.append(st)
    history = load_history()
    if a.apply:
        for st in studies:
            if st.changed and st.passed:          # a design that still misses its target is not applied
                mod = part_module(st.name)
                before = {k: p.cad[mod.SECTION][k] for k in st.changed}
                apply_to_spec(mod.SECTION, st.changed)
                history.append({"date": datetime.date.today().isoformat(), "part": st.name,
                                "changes": {k: [before[k], v] for k, v in st.changed.items()},
                                "sf_before": {k: round(v, 3) for k, v in st.iterations[0]["sf"].items()},
                                "sf_after": {k: round(v, 3) for k, v in st.iterations[-1]["sf"].items()},
                                "iterations": st.iterations})
                print(f"applied to spec/cad.yaml [{mod.SECTION}]: {st.changed}")
        HISTORY_PATH.write_text(json.dumps(history, indent=2))
    res = write_results(studies, anchors, spec_hash(), out)
    dump_iterations(studies, out / "iterations.json")
    write_report(studies, anchors, res, out, history)
    print(f"wrote {out / 'report.md'} and results.json; overall {'PASS' if res['all_passed'] else 'FAIL'}")
    return 0 if res["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
