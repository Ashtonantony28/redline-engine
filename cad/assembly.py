"""Engine assembly at a given crank angle (ADR 0001: 0 deg = firing TDC of cylinder 1; ADR 0002 rotation sense).

Every placed part is an *instance* of a source solid from `parts()` (built once, local frame) with a fixed base
location inside its *motion group* and the group's location at the crank angle:

* ``static``      block, head and everything bolted to them;
* ``crank``       crankshaft and what turns with it: rotates by -theta about +x (ADR 0002), so the crankpin of
                  cylinder k sits at (y, z) = r (sin(theta - phase_k), cos(theta - phase_k));
* ``rod_k``       conrod k: at its crankpin, rotated about the big-end axis by the exact slider-crank rod angle;
* ``piston_k``    piston assembly k at the exact slider-crank wrist-pin height (kin/slider_crank.py, M0-verified);
* ``cam_<kind>``  camshaft + pulley: rotates by -(theta + adv)/2 about its axis (adv = intake VVT advance, crank
                  degrees; the exhaust cam has no phaser);
* ``valve_<cyl>_<kind>_<n>``  valve, retainer, collets, shim and bucket: translate down the valve axis by the cam
                  lift of kin/cam.py at the cylinder's phase (firing order 1-3-4-2, ADR 0001/0002).

Springs are placed at their installed (valve-closed) length; their compression with lift is checked separately
(cad/audit.py) and animated by axial scaling in the viewer.

`placed_parts(theta)` returns {instance name: (located shape, material key)}; cad/audit.py uses `instances()` and
`motion()` directly to sweep the moving groups through the cycle without rebuilding any solid.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import math

from build123d import Compound, Location, Pos, Rot

from cad import bearings, front, head as headmod, manifolds, valvetrain as vt
from cad.block import block, head_bolt, head_bolt_xy, head_gasket, liner, oil_pan, pan_bolt, pan_bolt_xy
from cad.conrod import conrod_parts
from cad.crankshaft import crankshaft
from cad.maincap import main_cap_bolt, main_cap_solid
from cad.params import Params, params
from cad.piston import piston_parts
from kin.slider_crank import piston_position

IDENTITY = Location()
CRANK_PERIOD_DEG = 360.0
CYCLE_DEG = 720.0
HALF_TURN = 180.0
LIFT_SPEED = 1.0        # mm of travel per mm of lift (valve groups translate along their axis)


@dataclass(frozen=True)
class Inst:
    name: str                 # unique instance name, e.g. "piston_1"
    kind: str                 # key of the source solid in parts()
    group: str                # motion group (see module docstring)
    base: Location = field(default_factory=Location)   # placement inside the group frame
    cyl: int | None = None    # 1-based cylinder number, if the part belongs to one


@dataclass(frozen=True)
class Motion:
    """Rigid motion of one group: loc(theta_deg, adv_deg) -> Location. period 0 = static.

    Valve groups also carry `lift_loc(u)` (location at valve lift u mm) and `u_max`, so a valve part against a
    static part can be swept over its lift range instead of the crank cycle (cad/audit.py)."""
    group: str
    period: float
    fn: object = None
    lift_loc: object = None
    u_max: float = 0.0
    lift: object = None       # lift(theta_deg, adv_deg) -> mm, valve groups only
    vvt: bool = False         # the group moves with the intake phaser

    def loc(self, theta: float, adv: float = 0.0) -> Location:
        return IDENTITY if self.fn is None else self.fn(theta, adv)

    @property
    def static(self) -> bool:
        return self.period == 0.0


def valve_ids(p: Params):
    """[(cyl idx, kind, n, x_local)] for every valve, n = 1, 2 within the kind."""
    out = []
    for i in range(p.cylinders):
        cnt = {"intake": 0, "exhaust": 0}
        for kind, _s, xl, _y, _d in p.valve_layout():
            cnt[kind] += 1
            out.append((i, kind, cnt[kind], xl))
    return out


@lru_cache(maxsize=None)
def parts(p: Params | None = None) -> dict:
    """Part solids in their local frames: {kind: (shape, material key)}. Cached."""
    p = p or params()
    out = {name: v for name, v in piston_parts(p).items()}
    for k, shape in conrod_parts(p).items():
        out[k] = (shape, {"conrod": "conrod", "conrod_cap": "conrod_cap"}.get(k, "conrod_bolt"))
    out["big_end_shell_upper"] = (bearings.big_end_shell(p, True), "bearing_shell")
    out["big_end_shell_lower"] = (bearings.big_end_shell(p, False), "bearing_shell")
    out["crankshaft"] = (crankshaft(p), "crankshaft")
    out["crank_sprocket"] = (front.crank_sprocket(p), "crank_sprocket")
    out["damper_hub"] = (front.damper_hub(p), "damper_hub")
    out["damper_rubber"] = (front.damper_rubber(p), "damper_rubber")
    out["damper_ring"] = (front.damper_ring(p), "damper_ring")
    out["flywheel"] = (front.flywheel(p), "flywheel")
    out["ring_gear"] = (front.ring_gear(p), "ring_gear")
    for k in range(int(front.fc(p)["flywheel"]["bolts"])):
        out[f"flywheel_bolt_{k + 1}"] = (front.flywheel_bolt(p, k), "flywheel_bolt")
    out["block"] = (block(p), "block")
    out["liner"] = (liner(p), "liner")
    out["main_cap"] = (main_cap_solid(p, half=False), "main_cap")
    out["main_cap_bolt"] = (main_cap_bolt(p), "main_cap_bolt")
    for j in range(p.n_mains):
        out[f"main_shell_upper_{j + 1}"] = (bearings.main_shell(p, j, True), "bearing_shell")
        out[f"main_shell_lower_{j + 1}"] = (bearings.main_shell(p, j, False), "bearing_shell")
    out["head_gasket"] = (head_gasket(p), "head_gasket")
    out["head"] = (headmod.head(p), "head")
    out["head_bolt"] = (head_bolt(p), "head_bolt")
    out["oil_pan"] = (oil_pan(p), "oil_pan")
    out["pan_bolt"] = (pan_bolt(p), "pan_bolt")
    out["spark_plug"] = (headmod.spark_plug(p), "spark_plug")
    for kind in ("intake", "exhaust"):
        out[f"valve_{kind}"] = (vt.valve(p, kind), "valve")
        out[f"seat_insert_{kind}"] = (vt.seat_insert(p, kind), "seat_insert")
        out[f"guide_{kind}"] = (vt.guide(p, kind), "guide")
        out[f"spring_seat_{kind}"] = (vt.spring_seat(p, kind), "spring_seat")
        st = vt.stack(p, kind)
        out[f"spring_{kind}"] = (Pos(0, 0, st["spring_bottom"]) * vt.spring(p, vt.spring_located_length(p, kind, 0.0)),
                                 "valve_spring")
        out[f"retainer_{kind}"] = (vt.retainer(p, kind), "retainer")
        out[f"collet_{kind}"] = (vt.collet(p, kind), "collet")
        out[f"shim_{kind}"] = (vt.shim(p, kind), "shim")
        out[f"bucket_{kind}"] = (vt.bucket(p, kind), "bucket")
        out[f"camshaft_{kind}"] = (vt.camshaft(p, kind), "camshaft")
        out[f"cam_pulley_{kind}"] = (front.cam_pulley(p), "cam_pulley")
    out["cam_cap"] = (vt.cam_cap(p), "cam_cap")
    out["cam_cap_bolt"] = (vt.cam_cap_bolt(p), "cam_cap_bolt")
    for nm, (_c, r) in ((n, (c, r)) for n, c, r in front.belt_circles(p) if n in ("idler", "tensioner")):
        out[nm] = (front.idler_pulley(p, r), "idler")
    out["timing_belt"] = (front.belt(p), "belt")
    out["intake_manifold"] = (manifolds.intake_manifold(p), "intake_manifold")
    out["exhaust_header"] = (manifolds.exhaust_header(p), "exhaust_header")
    out["throttle_body"] = (manifolds.throttle_body(p), "throttle_body")
    out["intake_bracket"] = (manifolds.intake_bracket(p), "intake_bracket")
    out["bracket_bolt_block"] = (manifolds.bracket_bolt(p, "block"), "bracket_bolt")
    out["bracket_bolt_pad"] = (manifolds.bracket_bolt(p, "pad"), "bracket_bolt")
    out["exhaust_stay"] = (manifolds.exhaust_stay(p), "exhaust_stay")
    out["exhaust_stay_bolt_block"] = (manifolds.stay_bolt(p, "block"), "bracket_bolt")
    out["exhaust_stay_bolt_pad"] = (manifolds.stay_bolt(p, "pad"), "bracket_bolt")
    out["manifold_bolt_intake"] = (manifolds.manifold_bolt(p, "intake", 0.0, 0.0), "manifold_bolt")
    out["manifold_bolt_exhaust"] = (manifolds.manifold_bolt(p, "exhaust", 0.0, 0.0), "manifold_bolt")
    return out


def cylinder_pose(p: Params, idx: int, crank_deg: float):
    """(pin_y, pin_z, rod_angle_deg, wrist_pin_z) for cylinder idx at a crank angle."""
    a = math.radians(p.pin_angle_deg(idx, crank_deg))
    r, l = p.crank_r, p.rod_len
    py, pz = r * math.sin(a), r * math.cos(a)
    wz = float(piston_position(a, r, l))
    rod_deg = math.degrees(math.asin(py / l))  # Rot about x by rod_deg maps +z to (-sin, cos)
    return py, pz, rod_deg, wz


VALVE_MOVING = ("valve", "retainer", "collet", "collet_b", "shim", "bucket")
VALVE_STATIC = ("seat_insert", "guide", "spring_seat", "spring")


def valve_group(i: int, kind: str, n: int) -> str:
    return f"valve_{i + 1}_{kind}_{n}"


@lru_cache(maxsize=None)
def instances(p: Params | None = None) -> tuple:
    p = p or params()
    src = parts(p)
    S = "static"
    out = [Inst("block", "block", S), Inst("head_gasket", "head_gasket", S), Inst("head", "head", S),
           Inst("oil_pan", "oil_pan", S), Inst("crankshaft", "crankshaft", "crank")]
    for k in ("crank_sprocket", "damper_hub", "damper_rubber", "damper_ring", "flywheel", "ring_gear"):
        out.append(Inst(k, k, "crank"))
    out += [Inst(k, k, "crank") for k in src if k.startswith("flywheel_bolt_")]
    piston_kinds = list(piston_parts(p))
    for i in range(p.cylinders):
        c = i + 1
        x = p.cyl_x(i)
        out.append(Inst(f"liner_{c}", "liner", S, Pos(x, 0.0, 0.0), c))
        for k in ("conrod", "conrod_cap", "conrod_bolt_a", "conrod_bolt_b", "big_end_shell_upper",
                  "big_end_shell_lower"):
            out.append(Inst(f"{k}_{c}", k, f"rod_{c}", cyl=c))
        for n in piston_kinds:
            out.append(Inst(f"{n}_{c}", n, f"piston_{c}", Pos(x, 0.0, 0.0), c))
        out.append(Inst(f"spark_plug_{c}", "spark_plug", S, Pos(*headmod.plug_xy(p, i), 0.0), c))
    for j in range(p.n_mains):
        m = j + 1
        mx = p.main_x(j)
        out.append(Inst(f"main_cap_{m}", "main_cap", S, Pos(mx, 0.0, 0.0)))
        out.append(Inst(f"main_cap_bolt_{m}a", "main_cap_bolt", S, Pos(mx, 0.0, 0.0)))
        out.append(Inst(f"main_cap_bolt_{m}b", "main_cap_bolt", S, Pos(mx, 0.0, 0.0) * Rot(0, 0, HALF_TURN)))
        out.append(Inst(f"main_shell_upper_{m}", f"main_shell_upper_{m}", S))
        out.append(Inst(f"main_shell_lower_{m}", f"main_shell_lower_{m}", S))
    for k, (x, y) in enumerate(head_bolt_xy(p)):
        out.append(Inst(f"head_bolt_{k + 1}", "head_bolt", S, Pos(x, y, 0.0)))
    for k, (x, y) in enumerate(pan_bolt_xy(p)):
        out.append(Inst(f"pan_bolt_{k + 1}", "pan_bolt", S, Pos(x, y, 0.0)))
    for (i, kind, n, xl) in valve_ids(p):
        tag = f"{i + 1}_{kind}_{n}"
        frame = vt.valve_frame(p, i, kind, xl)
        g = valve_group(i, kind, n)
        for k in VALVE_MOVING:
            kk = "collet" if k == "collet_b" else k
            base = Rot(0, 0, HALF_TURN) if k == "collet_b" else IDENTITY
            out.append(Inst(f"{k}_{tag}", f"{kk}_{kind}", g, base, i + 1))
        for k in VALVE_STATIC:
            out.append(Inst(f"{k}_{tag}", f"{k}_{kind}", S, frame, i + 1))
    cc = vt.vc(p)["cam_cap"]
    for kind in ("intake", "exhaust"):
        out.append(Inst(f"camshaft_{kind}", f"camshaft_{kind}", f"cam_{kind}"))
        out.append(Inst(f"cam_pulley_{kind}", f"cam_pulley_{kind}", f"cam_{kind}"))
        cy, cz = vt.cam_axis_yz(p, kind)
        for j, x in enumerate(vt.journal_xs(p)):
            out.append(Inst(f"cam_cap_{kind}_{j + 1}", "cam_cap", S, Pos(x, cy, cz)))
            for s, ab in ((-1, "a"), (1, "b")):
                out.append(Inst(f"cam_cap_bolt_{kind}_{j + 1}{ab}", "cam_cap_bolt", S,
                                Pos(x, cy + s * p.L(cc["bolt_offset"]), cz)))
        for i in range(p.cylinders):
            for k, (bx, bz) in enumerate(headmod.port_bolt_xz(p, kind, i)):
                out.append(Inst(f"manifold_bolt_{kind}_{i + 1}{'ab'[k]}", f"manifold_bolt_{kind}", S, Pos(bx, 0.0, bz)))
    for nm, c, _r in front.belt_circles(p):
        if nm in ("idler", "tensioner"):
            out.append(Inst(nm, nm, S, Pos(0.0, c[0], c[1])))
    out += [Inst("timing_belt", "timing_belt", S), Inst("intake_manifold", "intake_manifold", S),
            Inst("exhaust_header", "exhaust_header", S), Inst("throttle_body", "throttle_body", S),
            Inst("intake_bracket", "intake_bracket", S)]
    for kind, pts in manifolds.bracket_bolts(p).items():
        for k, (bx, by, bz) in enumerate(pts):
            out.append(Inst(f"bracket_bolt_{kind}_{'ab'[k]}", f"bracket_bolt_{kind}", S, Pos(bx, by, bz)))
    out.append(Inst("exhaust_stay", "exhaust_stay", S))
    for kind, pts in manifolds.stay_bolts(p).items():
        for k, (bx, by, bz) in enumerate(pts):
            out.append(Inst(f"exhaust_stay_bolt_{kind}_{'ab'[k]}", f"exhaust_stay_bolt_{kind}", S, Pos(bx, by, bz)))
    assert all(it.kind in src for it in out)
    assert len({it.name for it in out}) == len(out)
    return tuple(out)


@lru_cache(maxsize=None)
def motions(p: Params | None = None) -> dict:
    p = p or params()
    out = {"static": Motion("static", 0.0), "crank": Motion("crank", CRANK_PERIOD_DEG, lambda th, adv: Rot(-th, 0, 0))}
    for i in range(p.cylinders):
        def rod(th, adv, i=i):
            py, pz, rod_deg, _ = cylinder_pose(p, i, th)
            return Pos(p.cyl_x(i), py, pz) * Rot(rod_deg, 0, 0)

        def pis(th, adv, i=i):
            return Pos(0.0, 0.0, cylinder_pose(p, i, th)[3])
        out[f"rod_{i + 1}"] = Motion(f"rod_{i + 1}", CRANK_PERIOD_DEG, rod)
        out[f"piston_{i + 1}"] = Motion(f"piston_{i + 1}", CRANK_PERIOD_DEG, pis)
    profiles = vt.cam_profiles(p)
    for kind in ("intake", "exhaust"):
        cy, cz = vt.cam_axis_yz(p, kind)
        vvt = kind == "intake"

        def cam(th, adv, cy=cy, cz=cz, vvt=vvt):
            return Pos(0.0, cy, cz) * Rot(-(th + (adv if vvt else 0.0)) / 2.0, 0, 0)
        out[f"cam_{kind}"] = Motion(f"cam_{kind}", CYCLE_DEG, cam, vvt=vvt)
    for (i, kind, n, xl) in valve_ids(p):
        frame = vt.valve_frame(p, i, kind, xl)
        prof = profiles[kind]
        phase = p.firing_phase_deg(i)
        vvt = kind == "intake"

        def lift(th, adv, prof=prof, phase=phase, vvt=vvt):
            a = math.radians(adv) if vvt else 0.0
            return p.L(float(prof.valve_lift_crank(math.radians(th - phase), a)))

        def vloc(th, adv, frame=frame, lift=lift):
            return frame * Pos(0.0, 0.0, -lift(th, adv))

        def lloc(u, frame=frame):
            return frame * Pos(0.0, 0.0, -u)
        g = valve_group(i, kind, n)
        out[g] = Motion(g, CYCLE_DEG, vloc, lift_loc=lloc, u_max=p.max_lift(kind), lift=lift, vvt=vvt)
    return out


def instance_shape(inst: Inst, p: Params, theta: float = 0.0, adv: float = 0.0):
    shape = parts(p)[inst.kind][0]
    return shape.moved(motions(p)[inst.group].loc(theta, adv) * inst.base)


def placed_parts(crank_deg: float = 0.0, p: Params | None = None, adv_deg: float = 0.0) -> dict:
    """{instance name: (located shape, material key)} at the given crank angle (and intake VVT advance)."""
    p = p or params()
    src = parts(p)
    return {it.name: (instance_shape(it, p, crank_deg, adv_deg), src[it.kind][1]) for it in instances(p)}


def assembly(crank_deg: float = 0.0, p: Params | None = None) -> Compound:
    children = []
    for name, (shape, _mat) in placed_parts(crank_deg, p).items():
        shape.label = name
        children.append(shape)
    return Compound(children=children, label=f"redline_engine_{crank_deg:g}deg")
