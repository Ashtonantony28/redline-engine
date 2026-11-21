#!/usr/bin/env python3
"""Check the Redline ECU Zener design against the firmware, ERC rules and datasheet limits.

Verify:  cd ecad && python -m pytest -q check_design.py      (or: python3 check_design.py)

What is checked
  build        `pcb +local build -D warnings` passes with no errors and no warnings
  firmware     every I/O of ecu/include/ecu/*.h (entry points, ecu_inputs_t fields, the throttle
               output, the INJ / DWELL / SPARK events per cylinder) has an io_map.yaml entry whose
               MCU pins, harness pins and conditioning chain exist on the netlist
  ERC          no single-pin nets except explicit no-connects, every pin on a net, IC supply pins on
               the right rails, every MCU pin is either in the pin plan or explicitly open
  limits       datasheet numbers (modules/parts.zen LIMITS, re-read here) against the netlist and the
               spec YAML: TVS clamp vs every part on VBAT, driver currents vs ratings with the loads
               of spec/mvem.yaml and spec/ecu.yaml, ADC full scale vs VDDA, logic levels
  library      every symbol and footprint is a file under ecad/lib or from the pcb stdlib

Environment
  PCB    pcb binary (default: `pcb` on PATH, then ~/.local/bin/pcb); run with the locally built
         toolchain (`pcb +local`), see docs/ecad.md
"""
import ast
import collections
import functools
import json
import os
import re
import shutil
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BOARD = os.path.join(HERE, "boards", "redline_ecu", "redline_ecu.zen")
PCB = os.environ.get("PCB") or shutil.which("pcb") or os.path.expanduser("~/.local/bin/pcb")
TOOLCHAIN = os.environ.get("PCB_TOOLCHAIN", "+local")
ECU_INC = os.path.join(REPO, "ecu", "include", "ecu")

Component = collections.namedtuple("Component", "ref path attrs pins pads")


def available():
    return os.path.isfile(PCB) and os.access(PCB, os.X_OK)


def _run(*args):
    return subprocess.run([PCB, TOOLCHAIN, *args], cwd=HERE, capture_output=True, text=True)


def _attr(v):
    return next(iter(v.values())) if isinstance(v, dict) and len(v) == 1 else v


@functools.lru_cache(maxsize=None)
def build():
    r = _run("build", "-D", "warnings", BOARD)
    return r.returncode, r.stdout + r.stderr


@functools.lru_cache(maxsize=None)
def netlist():
    r = _run("build", "--netlist", BOARD)
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    root = data["root_ref"] + "."
    port_net, nets = {}, {}
    for name, net in data["nets"].items():
        ports = [p[len(root):] for p in net["ports"]]
        nets[name] = {"kind": net["kind"], "ports": ports}
        for p in ports:
            assert p not in port_net, "port %s is on two nets" % p
            port_net[p] = name
    comps = {}
    for key, inst in data["instances"].items():
        if inst.get("kind") != "Component":
            continue
        path = key[len(root):]
        pins, pads = {}, {}
        for pin, child in inst["children"].items():
            net = port_net.get(path + "." + pin)
            pins[pin] = net
            for pad in _attr(data["instances"][child]["attributes"]["pads"]):
                pads[_attr(pad)] = net
        ref = inst["reference_designator"]
        assert ref not in comps, "duplicate designator " + ref
        attrs = {k: _attr(v) for k, v in inst["attributes"].items() if not k.startswith("__")}
        comps[ref] = Component(ref, path, attrs, pins, pads)
    return nets, comps


def io_map():
    with open(os.path.join(HERE, "io_map.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)


def spec(name):
    with open(os.path.join(REPO, "spec", name), encoding="utf-8") as f:
        return yaml.safe_load(f)


def zen_dict(module, name):
    """A dict literal assigned at top level of a .zen file (Starlark literals are Python literals)."""
    tree = ast.parse(open(os.path.join(HERE, "modules", module), encoding="utf-8").read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


LIMITS = zen_dict("parts.zen", "LIMITS")
LOADS = zen_dict("parts.zen", "LOADS")


def value_ohms(text):
    v = text.split()[0].lower().replace("ohm", "")
    mult = {"k": 1e3, "m": 1e6}.get(v[-1:], 1.0)
    return float(v[:-1] if v[-1:] in "km" else v) * mult


def resistors_between(a, b):
    _, comps = netlist()
    return [c for c in comps.values() if c.attrs.get("type") == "resistor" and sorted(c.pins.values()) == sorted([a, b])]


def by_mpn(prefix):
    _, comps = netlist()
    return [c for c in comps.values() if str(c.attrs.get("mpn", "")).startswith(prefix)]


def pad_net(ref, pad):
    _, comps = netlist()
    return comps[ref].pads[pad]


# ---------------------------------------------------------------------------
# firmware I/O inventory, parsed from the headers
# ---------------------------------------------------------------------------


def firmware_io():
    """Names of the firmware-facing I/O the board has to provide."""
    ecu_h = open(os.path.join(ECU_INC, "ecu.h"), encoding="utf-8").read()
    types_h = open(os.path.join(ECU_INC, "ecu_types.h"), encoding="utf-8").read()
    calib = open(os.path.join(ECU_INC, "ecu_calib_gen.h"), encoding="utf-8").read()
    io = {}
    for fn, key in (("ecu_crank_edge", "crank"), ("ecu_cam_edge", "cam"), ("ecu_knock_window", "knock")):
        assert re.search(r"void %s\(" % fn, ecu_h), fn
        io[fn] = key
    inputs = re.search(r"typedef struct \{([^{}]*)\} ecu_inputs_t;", ecu_h).group(1)
    for field in re.findall(r"^\s*(?:float|uint32_t)\s+(\w+);", inputs, re.M):
        if field != "now_tick":  # timer read, not a pin
            io["ecu_inputs_t." + field] = field
    io["ecu_outputs_t.etc_voltage"] = "etc_voltage"
    events = dict(re.findall(r"(ECU_EV_\w+)\s*=\s*(\d+)", types_h))
    assert set(events) == {"ECU_EV_INJ", "ECU_EV_DWELL", "ECU_EV_SPARK"}, events
    ncyl = int(re.search(r"#define ECU_NCYL (\d+)", calib).group(1))
    io["ECU_EV_INJ"] = "inj"
    io["ECU_EV_DWELL"] = "ign"
    io["ECU_EV_SPARK"] = "ign"
    return io, ncyl


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------


def test_build_clean():
    rc, out = build()
    assert rc == 0, out
    low = out.lower()
    assert "warning" not in low and "error" not in low, out


def _check_entry(name, entry):
    nets, comps = netlist()
    mcu = comps["U1"]
    for pin, net in entry.get("mcu", {}).items():
        assert mcu.pins.get(pin) == net, "%s: U1.%s is on %s, io_map says %s" % (name, pin, mcu.pins.get(pin), net)
    for loc, net in entry.get("harness", {}).items():
        ref, pad = loc.split(".")
        assert comps[ref].pads.get(pad) == net, "%s: %s pad %s is on %s, io_map says %s" % (
            name, ref, pad, comps[ref].pads.get(pad), net)
    chains = entry.get("chains") or [entry.get("chain", [])]
    by_path = {c.path: c for c in comps.values()}
    for chain in chains:
        for n in chain:
            assert n in nets, "%s: net %s not in the netlist" % (name, n)
        for a, b in zip(chain, chain[1:]):
            owners_a = {p.rsplit(".", 1)[0] for p in nets[a]["ports"]}
            owners_b = {p.rsplit(".", 1)[0] for p in nets[b]["ports"]}
            shared = owners_a & owners_b
            assert shared, "%s: no component links %s and %s" % (name, a, b)
            assert all(by_path[s].ref for s in shared)


def test_firmware_io_mapped():
    """Every firmware I/O in ecu/include/ecu/*.h is mapped to an MCU pin, a net and a harness pin."""
    io, ncyl = firmware_io()
    m = io_map()
    fw = m["firmware"]
    sources = {e["source"].split("(")[0].split(" ")[0] for e in fw.values()}
    for symbol, key in io.items():
        assert key in fw, "firmware I/O %s (%s) has no io_map entry" % (symbol, key)
        src = fw[key]["source"]
        bare = symbol.split(".")[-1]
        assert bare in src or symbol in src, "%s: io_map source '%s' does not name it" % (symbol, src)
    for key in ("inj", "ign"):
        e = fw[key]
        assert e["per_cylinder"] and len(e["mcu"]) == ncyl and len(e["harness"]) == ncyl and len(e["chains"]) == ncyl
    for name, entry in list(fw.items()) + list(m["board"].items()):
        _check_entry(name, entry)
    _check_entry("timebase", m["timebase"])
    assert sources  # parsed


def test_mcu_pin_plan_complete():
    """Each MCU pin is a mapped function, a power pin, or an explicit no-connect; mapped pins are io_map pins."""
    nets, comps = netlist()
    mcu = comps["U1"]
    m = io_map()
    mapped = {}
    for entry in list(m["firmware"].values()) + list(m["board"].values()) + [m["timebase"]]:
        mapped.update(entry.get("mcu", {}))
    rails = {"3V3", "VDDA", "GND", "MCU.VCAP_1", "MCU.VCAP_2"}
    support = {"PH0": "MCU.HSE_IN", "PH1": "MCU.HSE_OUT", "BOOT0": "BOOT0", "PE3": "MCU.LED_STATUS",
               "PE4": "MCU.LED_FAULT"}
    for pin, net in mcu.pins.items():
        if net in rails:
            continue
        if pin in mapped:
            assert net == mapped[pin]
        elif pin in support:
            assert net == support[pin], (pin, net)
        else:
            assert nets[net]["kind"] == "NotConnected", "U1.%s is on %s but not in io_map.yaml" % (pin, net)
    assert len(mcu.pads) == 100


def test_erc_nets():
    """No floating inputs: single-pin nets only as explicit no-connects; every pin on a net."""
    nets, comps = netlist()
    for name, net in nets.items():
        if net["kind"] == "NotConnected":
            assert len(net["ports"]) == 1, "%s is NotConnected but has %d pins" % (name, len(net["ports"]))
        else:
            assert len(net["ports"]) >= 2, "single-pin net %s (%s)" % (name, net["ports"])
    for c in comps.values():
        for pin, net in c.pins.items():
            assert net is not None, "%s.%s is not on any net" % (c.ref, pin)
    by_path = {c.path: c for c in comps.values()}
    stubs = set()
    for name, net in nets.items():
        owners = [by_path[p.rsplit(".", 1)[0]].ref for p in net["ports"]]
        if len(owners) == 2 and any(o.startswith("TP") for o in owners):
            stubs.add(name)
    assert stubs == set(), stubs


SUPPLY_PINS = {
    "U1": {"VDD": "3V3", "VBAT": "3V3", "VDDA": "VDDA", "VREF+": "VDDA", "VSS": "GND", "VSSA": "GND"},
    "U2": {"ANODE": "VBAT_IN", "CATHODE": "VBAT", "GND": "GND"},
    "U3": {"VIN": "VBAT", "EN": "VBAT", "PGND": "GND", "AGND": "GND", "VCC": "PWR.BUCK_VCC"},
    "U4": {"VIN": "5V", "VOUT": "3V3", "GND": "GND"},
    "U5": {"VIN": "VBAT", "ADJ_EN": "5V", "VOUT": "5V_S1", "GND": "GND"},
    "U6": {"VIN": "VBAT", "ADJ_EN": "5V", "VOUT": "5V_S2", "GND": "GND"},
    "U7": {"VM": "VBAT", "GND": "GND", "SRC": "GND", "DVDD": "ETB.DVDD", "MODE": "ETB.DVDD"},
    "U8": {"VDD": "5V", "GND": "GND"},
    "U9": {"VDD": "5V", "GND": "GND"},
    "U10": {"VDD": "5V", "GND": "GND"},
    "U11": {"VCC": "5V", "GND": "GND"},
    "U12": {"VDD": "5V", "GND": "GND"},
    "U13": {"UB": "VBAT", "VCC": "5V", "VCCS": "5V", "VCCSS": "5V", "GND": "GND", "GNDS": "GND"},
    "U14": {"VCC": "5V", "GND": "GND"},
    "U15": {"VCC": "5V", "VIO": "3V3", "GND": "GND", "STB": "GND"},
    "U16": {"GND": "GND"}, "U17": {"GND": "GND"}, "U18": {"GND": "GND"}, "U19": {"GND": "GND"},
    "U20": {"GND": "GND"}, "U21": {"GND": "GND"},
    "U22": {"VCC": "5V", "GND": "GND"},
}


def test_supply_rails():
    nets, comps = netlist()
    for ref, pins in SUPPLY_PINS.items():
        for pin, rail in pins.items():
            assert comps[ref].pins[pin] == rail, "%s.%s: expected %s, got %s" % (ref, pin, rail, comps[ref].pins[pin])
    ics = {ref for ref in comps if re.fullmatch(r"U\d+", ref)}
    assert ics == set(SUPPLY_PINS), "ICs without a supply check: %s" % sorted(ics ^ set(SUPPLY_PINS))
    for rail in ("VBAT_IN", "VBAT", "5V", "3V3", "VDDA", "5V_S1", "5V_S2", "VBUS"):
        assert nets[rail]["kind"] == "Power", rail
    assert nets["GND"]["kind"] == "Ground"
    # every IC has local decoupling: a capacitor from one of its supply rails to GND in its module
    caps = [c for c in comps.values() if c.attrs.get("type") == "capacitor" or c.ref in ("C1", "C2")]
    for ref in ("U1", "U3", "U4", "U5", "U6", "U7", "U8", "U9", "U10", "U11", "U12", "U13", "U14", "U15", "U22"):
        rails = {n for n in comps[ref].pins.values() if nets[n]["kind"] == "Power"}
        assert any(set(c.pins.values()) == {r, "GND"} for c in caps for r in rails), "%s has no decoupling" % ref


def _parts_on(net):
    nets, comps = netlist()
    by_path = {c.path: c for c in comps.values()}
    refs = {by_path[p.rsplit(".", 1)[0]].ref: by_path[p.rsplit(".", 1)[0]] for p in nets[net]["ports"]}
    return list(refs.values())


ABS_MAX_ON_VBAT = {  # mpn prefix -> LIMITS key with the absolute maximum of the VBAT-connected pin
    "LM74700": ("LM74700-Q1", "anode_max"),
    "LMR36015": ("LMR36015-Q1", "vin_abs_max"),
    "TPS7B4250": ("TPS7B4250-Q1", "vin_abs_max"),
    "DRV8873": ("DRV8873H-Q1", "vm_abs_max"),
    "CJ125": ("CJ125", "ub_abs_max"),
    "DMTH6016": ("DMTH6016LK3Q", "vds_max"),
    "SBR3U60": ("SBR3U60P1", "vr_max"),
}


def _rated_voltage(c):
    v = str(c.attrs.get("voltage", "")).upper().replace("V", "").strip()
    try:
        return float(v.split()[0])
    except (ValueError, IndexError):
        return None


def test_limit_tvs_vs_battery_parts():
    """The TVS clamp (and the jump start) stays inside the absolute maximum of every part on VBAT/VBAT_IN."""
    tvs = LIMITS["SMCJ24CA"]
    assert tvs["vrwm"] >= LOADS["vbat_jump"] and tvs["vbr_min"] > LOADS["vbat_jump"]
    _, comps = netlist()
    assert comps["D1"].attrs["mpn"] == "SMCJ24CA" and set(comps["D1"].pins.values()) == {"VBAT_IN", "GND"}
    checked = set()
    for net in ("VBAT", "VBAT_IN"):
        for c in _parts_on(net):
            mpn = str(c.attrs.get("mpn", ""))
            key = next((v for k, v in ABS_MAX_ON_VBAT.items() if mpn.startswith(k)), None)
            if key:
                assert tvs["vc"] <= LIMITS[key[0]][key[1]], "%s %s: clamp %s V > %s V" % (c.ref, mpn, tvs["vc"], LIMITS[key[0]][key[1]])
                checked.add(key[0])
            elif c.attrs.get("type") == "capacitor" or c.ref in ("C1", "C2"):
                rated = _rated_voltage(c) if c.ref not in ("C1", "C2") else float(c.attrs["value"].split()[1][:-1])
                assert rated and rated >= tvs["vc"], "%s on %s rated %s V < clamp %s V" % (c.ref, net, rated, tvs["vc"])
            elif c.attrs.get("type") == "resistor" or c.ref.startswith(("TP", "J", "D1", "Q1")):
                continue
            else:
                raise AssertionError("%s (%s) on %s has no voltage-limit entry" % (c.ref, mpn, net))
    assert checked >= {"LM74700-Q1", "LMR36015-Q1", "TPS7B4250-Q1", "DRV8873H-Q1", "CJ125", "SBR3U60P1"}, checked
    # the low-side switches see VBAT through their loads: clamp must not trigger on a load dump
    assert LIMITS["ZXMS6006DGQ"]["vclamp_min"] > tvs["vc"]


def test_limit_low_side_drivers():
    """Injector / VVT / relay currents within ZXMS6006DGQ ratings; flyback handled; harness pins reached."""
    zx = LIMITS["ZXMS6006DGQ"]
    mv = spec("mvem.yaml")
    vmax = LOADS["vbat_max_operating"]
    assert mv["injector"]["battery_voltage"] <= vmax
    _, comps = netlist()
    q = {c.ref: c for c in by_mpn("ZXMS6006")}
    assert sorted(q) == ["Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]
    loads = {"Q2": "injector", "Q3": "injector", "Q4": "injector", "Q5": "injector", "Q6": "vvt", "Q7": "relay"}
    outs = {"Q2": "INJ1", "Q3": "INJ2", "Q4": "INJ3", "Q5": "INJ4", "Q6": "VVT", "Q7": "PUMP"}
    for ref, kind in loads.items():
        i = vmax / LOADS[kind + "_r_min"]
        assert i <= zx["id_nom"] and i < zx["id_lim_min_3v"], "%s %s current %.2f A" % (ref, kind, i)
        assert q[ref].pins["D"] == outs[ref] and q[ref].pins["S"] == "GND"
        # 3.3 V MCU drive is above the input threshold and below VIN max
        assert zx["vin_th_max"] < 3.0 < zx["vin_max"]
    assert vmax <= zx["vds_sc"]
    # injectors: active clamp only (fast closing); energy per closing far below EAS
    i_inj = vmax / LOADS["injector_r_min"]
    e = 0.5 * LOADS["injector_l"] * i_inj ** 2
    assert e < 0.05 * zx["eas"], e
    for net in ("INJ1", "INJ2", "INJ3", "INJ4"):
        assert not [c for c in _parts_on(net) if c.attrs.get("type") == "rectifier" or c.ref.startswith("D")], net
    # VVT: freewheel Schottky to VBAT sized for the solenoid current
    d3 = comps["D3"]
    assert d3.pins == {"A": "VVT", "K": "VBAT"}
    assert vmax / LOADS["vvt_r_min"] <= LIMITS["SBR3U60P1"]["if_max"]


def test_limit_ignition_and_heater():
    ucc = LIMITS["UCC27524A-Q1"]
    _, comps = netlist()
    for n in range(1, 5):
        r = resistors_between("OUT.IGN%d_DRV" % n, "IGN%d" % n)
        assert len(r) == 1 and r[0].attrs["package"] == "1206", n
        ohms = value_ohms(r[0].attrs["value"])
        assert 5.25 / ohms <= ucc["iout_dc_max"], "IGN%d short-circuit current" % n
        assert 5.25 / (ohms + LOADS["coil_input_r_min"]) > 0.002, "IGN%d drive too weak for the coil input" % n
    assert ucc["vdd_min"] <= 4.75 and ucc["vin_h_max"] < 3.0
    dm = LIMITS["DMTH6016LK3Q"]
    assert LOADS["heater_i_max"] <= 0.5 * dm["id_ta25"]
    assert LOADS["heater_i_max"] ** 2 * dm["rds_4v5"] < 1.0  # W, conduction loss in the DPAK
    assert comps["Q8"].pins == {"G": "KL.HTR_GATE", "D": "LSU_HTR_N", "S": "GND"}


def test_limit_etb():
    """Throttle motor stall current (spec/mvem.yaml, spec/ecu.yaml) vs DRV8873H-Q1 ITRIP/OCP; IPROPI scaling."""
    d = LIMITS["DRV8873H-Q1"]
    r_motor = spec("mvem.yaml")["etc"]["motor_resistance"]
    v_cmd = spec("ecu.yaml")["throttle"]["max_voltage"]
    i_stall = v_cmd / r_motor
    assert i_stall < d["iocp_min"], i_stall
    r = resistors_between("ETB.IPROPI", "GND")
    assert len(r) == 1
    v_itrip = d["itrip_max"] / d["k_mirror"] * value_ohms(r[0].attrs["value"])
    assert v_itrip <= 3.0, v_itrip
    assert d["iocp_min"] / d["k_mirror"] * value_ohms(r[0].attrs["value"]) < d["ipropi_abs_max"]
    _, comps = netlist()
    u7 = comps["U7"]
    assert u7.pins["nITRIP"] == "GND"  # current regulation on
    assert pad_net("U7", "13") == pad_net("U7", "20") == "VBAT"
    assert {pad_net("U7", p) for p in ("18", "19")} == {"ETB_M_P"} and {pad_net("U7", p) for p in ("14", "15")} == {"ETB_M_N"}
    # two harness contacts per motor lead (Micro-Fit 5 A per contact)
    assert [comps["J1"].pads[p] for p in ("13", "14", "15", "16")] == ["ETB_M_P", "ETB_M_P", "ETB_M_N", "ETB_M_N"]
    assert i_stall / 2 <= 5.0
    # fail-safe: bridge disabled while the MCU is in reset
    assert resistors_between("ETB_DISABLE", "3V3") and resistors_between("ETB_NSLEEP", "GND")


def test_limit_adc_full_scale():
    """Every analog input stays below VDDA min (3.2 V) at its source's maximum."""
    vdda_min = 3.2
    for n in ("TPS1", "TPS2", "PPS1", "PPS2", "MAP"):
        top = resistors_between(n, n + "_ADC")
        bot = resistors_between(n + "_ADC", "GND")
        assert len(top) == 1 and len(bot) == 1, n
        rt, rb = value_ohms(top[0].attrs["value"]), value_ohms(bot[0].attrs["value"])
        assert 5.3 * rb / (rt + rb) <= vdda_min, n
    for mon, rail in (("S5V1_MON", "5V_S1"), ("S5V2_MON", "5V_S2")):
        rt = value_ohms(resistors_between(rail, mon)[0].attrs["value"])
        rb = value_ohms(resistors_between(mon, "GND")[0].attrs["value"])
        assert 5.3 * rb / (rt + rb) <= vdda_min
    rt = value_ohms(resistors_between("VBAT", "VBAT_SENSE")[0].attrs["value"])
    rb = value_ohms(resistors_between("VBAT_SENSE", "GND")[0].attrs["value"])
    assert 18.0 * rb / (rt + rb) <= vdda_min  # VBAT read to 18 V
    rt = value_ohms(resistors_between("KL.KNK_OUT", "KNOCK_ADC")[0].attrs["value"])
    rb = value_ohms(resistors_between("KNOCK_ADC", "GND")[0].attrs["value"])
    assert LIMITS["TPIC8101"]["vout_max"] * rb / (rt + rb) <= vdda_min
    for out in ("LAMBDA_UA", "LAMBDA_UR"):
        bot = value_ohms(resistors_between(out, "GND")[0].attrs["value"])
        top = [c for c in _parts_on(out) if c.attrs.get("type") == "resistor" and "GND" not in c.pins.values()]
        rt = value_ohms(top[0].attrs["value"])
        assert 5.0 * bot / (rt + bot) <= vdda_min, out
    rt = value_ohms(resistors_between("SENS.VR_OUT", "CRANK")[0].attrs["value"])
    rb = value_ohms(resistors_between("CRANK", "GND")[0].attrs["value"])
    assert 5.25 * rb / (rt + rb) <= 3.6  # digital input


def test_limit_crank_and_knock_interfaces():
    lm = LIMITS["LM1815"]
    r = resistors_between("CRANK_P", "SENS.VR_IN")
    assert len(r) == 1 and r[0].attrs["package"] == "1206"
    assert 66.0 / value_ohms(r[0].attrs["value"]) <= lm["iin_max"]
    # 60-2 wheel at redline: tooth period in 1 MHz ticks is far above the 1k/1n input filter
    eng = spec("engine.yaml")
    rpm = next(v for k, v in _walk(eng) if k == "redline_rpm")
    tooth_s = 60.0 / rpm / 60.0
    assert tooth_s > 50 * 1e3 * 1e-9
    # TPIC8101 inputs need 5 V logic: every 5 V-side SPI net is driven by an AHCT125 output
    _, comps = netlist()
    u12 = comps["U12"]
    t = LIMITS["TPIC8101"]
    assert t["vih_frac"] * 4.75 > 3.3
    for pin in ("SCLK", "SDI", "CS", "INT_HOLD"):
        net = u12.pins[pin]
        drivers = [c for c in _parts_on(net) if str(c.attrs.get("mpn", "")).startswith("SN74AHCT125")]
        assert drivers, "TPIC8101 %s (%s) is not level shifted" % (pin, net)
    assert resistors_between("KL.MISO_5V", "SPI2_MISO")  # 5 V SDO into a 5 V-tolerant pin through 1k


def _walk(d, prefix=""):
    if isinstance(d, dict):
        for k, v in d.items():
            yield k, v
            yield from _walk(v)


def test_connectors():
    _, comps = netlist()
    j1 = comps["J1"].pads
    j2 = comps["J2"].pads
    assert [j1["1"], j1["2"]] == ["VBAT_IN", "VBAT_IN"] and {j1[p] for p in ("3", "4", "23", "24")} == {"GND"}
    assert len(j1) == 24 and len(j2) == 24
    assert all(n is not None for n in list(j1.values()) + list(j2.values()))


def test_library_is_local():
    _, comps = netlist()
    for c in comps.values():
        for key in ("symbol_path", "footprint"):
            uri = c.attrs[key]
            assert uri.startswith(("package://workspace/lib/", "package://stdlib/")), "%s %s" % (c.ref, uri)
            if uri.startswith("package://workspace/"):
                assert os.path.isfile(os.path.join(HERE, uri[len("package://workspace/"):])), uri


def summary():
    nets, comps = netlist()
    kinds = collections.Counter(re.match(r"[A-Z]+", r).group(0) for r in comps)
    return {"components": len(comps), "nets": len(nets), "by_prefix": dict(sorted(kinds.items()))}


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS", name)
            except AssertionError as e:
                failed += 1
                print("FAIL", name, "-", e)
    print(summary())
    sys.exit(1 if failed else 0)
