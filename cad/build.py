"""Build every part and the assembly; export STEP/STL/glTF and mass properties.

    python -m cad.build [--crank-deg 0] [--out out/cad]

Writes <out>/<part>.{step,stl,glb} (mm), the assembly at the requested crank angle, and
<out>/mass_properties.json (SI, ADR 0001). All numbers are modelled: there is no physical engine.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from build123d import Compound, Unit, export_gltf, export_step, export_stl

from cad.assembly import assembly, parts
from cad.massprops import mass_props
from cad.measure import (clearance_gas_volume, conrod_centre_distance, counterweight_moments,
                         crank_throws)
from cad.params import params
from cad.piston import piston_parts
from redline.spec import ROOT

DEFAULT_OUT = ROOT / "out" / "cad"
CRANK_AXIS_POINT = (0.0, 0.0, 0.0)  # crank axis is the x axis through the origin (m)
ROTATING_PARTS = ("crankshaft",)


def rel_err(value: float, budget: float) -> float:
    return value / budget - 1.0


def mass_report(p=None) -> dict:
    p = p or params()
    src = parts(p)
    masses = p.engine["masses"]
    report: dict = {"note": "modelled values from CAD; there is no physical engine",
                    "units": "SI (kg, m, kg m^2)", "parts": {}}
    mps = {}
    for name, (shape, mat) in src.items():
        mp = mass_props(shape, p.density(mat), p.k, name)
        mps[name] = mp
        entry = {"material": p.cad["materials"][mat], "density_kgm3": p.density(mat), **mp.to_dict()}
        if name in ROTATING_PARTS:
            entry["inertia_about_crank_axis_kgm2"] = mp.inertia_about(CRANK_AXIS_POINT).tolist()
        report["parts"][name] = entry

    pa = None
    for name in piston_parts(p):
        pa = mps[name] if pa is None else pa + mps[name]
    report["piston_assembly"] = {**pa.to_dict(), "budget_kg": masses["piston_assembly"],
                                 "rel_err": rel_err(pa.mass, masses["piston_assembly"])}

    rod = mps["conrod"]
    L = p.geom.conrod_length
    se_frac_cad = rod.com[2] / L  # big end at z = 0, small end at z = L
    report["conrod"] = {
        "mass_kg": rod.mass, "budget_kg": masses["conrod"], "rel_err": rel_err(rod.mass, masses["conrod"]),
        "centre_distance_m": p.to_si_length(conrod_centre_distance(p, src["conrod"][0])),
        "spec_centre_distance_m": L,
        "com_from_big_end_m": float(rod.com[2]),
        "small_end_fraction_cad": float(se_frac_cad),
        "small_end_fraction_spec": masses["conrod_small_end_fraction"],
        "inertia_about_com_kgm2": rod.inertia_com.tolist(),
    }
    recip_spec = pa.mass + masses["conrod_small_end_fraction"] * rod.mass
    recip_budget = masses["piston_assembly"] + masses["conrod_small_end_fraction"] * masses["conrod"]
    report["reciprocating_mass"] = {"mass_kg": recip_spec, "budget_kg": recip_budget,
                                    "rel_err": rel_err(recip_spec, recip_budget),
                                    "mass_kg_with_cad_small_end_fraction": pa.mass + se_frac_cad * rod.mass}

    crank = mps["crankshaft"]
    ixx = float(crank.inertia_about(CRANK_AXIS_POINT)[0, 0])
    throws = [p.to_si_length(t) for t in crank_throws(p, src["crankshaft"][0])]
    report["crankshaft"] = {
        "mass_kg": crank.mass,
        "inertia_about_crank_axis_kgm2": ixx,
        "budget_kgm2": masses["crank_inertia"],
        "rel_err": rel_err(ixx, masses["crank_inertia"]),
        "throw_m": {"min": min(throws), "max": max(throws), "spec": p.geom.crank_radius},
        "balance": counterweight_moments(p, crank, p.density("crankshaft")),
    }

    v_cad = clearance_gas_volume(p, src["head"][0], src["piston"][0]) / p.k**3
    v_geo = p.geom.clearance_volume
    report["combustion_chamber"] = {
        "roof_angle_deg": p.roof_angle_deg,
        "clearance_volume_cad_m3": v_cad,
        "clearance_volume_geometric_m3": v_geo,
        "rel_err": rel_err(v_cad, v_geo),
        "compression_ratio_cad": (p.geom.swept_volume + v_cad) / v_cad,
        "crown_to_head_face_m": p.to_si_length(p.head_face_z - p.crown_tdc_z),
    }
    return report


def export(shape, stem: Path) -> list[str]:
    written = []
    export_step(shape, str(stem.with_suffix(".step")), unit=Unit.MM)
    written.append(stem.with_suffix(".step").name)
    export_stl(shape, str(stem.with_suffix(".stl")))
    written.append(stem.with_suffix(".stl").name)
    try:
        export_gltf(shape, str(stem.with_suffix(".glb")), unit=Unit.MM, binary=True)
        written.append(stem.with_suffix(".glb").name)
    except Exception as exc:  # glTF is optional: report, do not fail the build
        print(f"glTF export skipped for {stem.name}: {exc}")
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crank-deg", type=float, default=0.0, help="assembly crank angle (ADR 0001)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    p = params()

    files = []
    for name, (shape, _mat) in parts(p).items():
        files += export(shape, out / name)
    piston_set = Compound(children=[s for s, _m in piston_parts(p).values()], label="piston_assembly")
    files += export(piston_set, out / "piston_assembly")
    files += export(assembly(args.crank_deg, p), out / f"assembly_{args.crank_deg:g}deg")

    report = mass_report(p)
    report["assembly_crank_deg"] = args.crank_deg
    (out / "mass_properties.json").write_text(json.dumps(report, indent=2, default=_np))
    files.append("mass_properties.json")

    pa, rod, cr, ch = (report[k] for k in ("piston_assembly", "conrod", "crankshaft", "combustion_chamber"))
    print(f"wrote {len(files)} files to {out}  (modelled values)")
    print(f"piston assembly {pa['mass_kg']:.4f} kg vs {pa['budget_kg']} ({pa['rel_err']:+.2%})")
    print(f"conrod          {rod['mass_kg']:.4f} kg vs {rod['budget_kg']} ({rod['rel_err']:+.2%})")
    print(f"crank inertia   {cr['inertia_about_crank_axis_kgm2']:.5f} kg m^2 vs {cr['budget_kgm2']} "
          f"({cr['rel_err']:+.2%}), static imbalance {cr['balance']['crank_static_imbalance_kgm']:.2e} kg m")
    print(f"clearance vol   {ch['clearance_volume_cad_m3'] * 1e6:.3f} cm^3 vs "
          f"{ch['clearance_volume_geometric_m3'] * 1e6:.3f} ({ch['rel_err']:+.3%}), roof {ch['roof_angle_deg']:.2f} deg")
    return 0


def _np(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


if __name__ == "__main__":
    raise SystemExit(main())
