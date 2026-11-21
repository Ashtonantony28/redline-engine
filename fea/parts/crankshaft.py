"""Crankshaft: single-throw model, fillet stresses under combined bending + torsion from the dyn/
crankpin forces over the whole WOT cycle, Goodman fatigue; FEA fillet Kt vs a Peterson-type estimate.

Model (half throw by symmetry about the pin mid plane x = 0, see fea.geometry.crank_throw_solid):
the throw spans main-bearing centre plane to centre plane; the rear journal end face is a rigid body
(*RIGID BODY) at the bearing centre, i.e. the throw is simply supported between its two mains, the
statically determinate split dyn.bearings uses. Load cases (full-throw quantities, unit magnitude):
  S1  crankpin radial force, outward along the throw (+z)        symmetric BCs (ux = 0 on x = 0)
  S2  crankpin tangential force, rotation direction (+y, ADR 0002) symmetric; torque reacted half at each main
  S3  through torque tau = 1 N m in the driving sense            antisymmetric BCs (uy = uz = 0 on x = 0)
  S4  centrifugal load of the throw itself at omega^2 = 1 rad^2/s^2 (counterweights included)
The rod load enters through a distributing coupling over the central band of the pin.
For throw k with front torque T_in (sum of the torques of the cylinders ahead, dyn.forces) the actual
loading (front main free in torsion, rear main carries T_in + F_t r) is
  rear half  : Fr S1 + Ft S2 + (T_in + Ft r / 2) S3 + omega^2 S4
  front half : the mirror image of Fr S1 + Ft S2 - (T_in + Ft r / 2) S3 + omega^2 S4
(symmetric fields mirror unchanged, the antisymmetric torsion field changes sign).
"""
from __future__ import annotations

import math

import numpy as np

from fea import ccx
from fea.common import LevelResult, cosine_pressure, plane_nodes, solve, workdir
from fea.config import material
from fea.fatigue import goodman_from_history, goodman_sf, peterson_kt_bending, signed_von_mises
from fea.geometry import crank_throw_solid, export, throw_dims, _halfspace
from fea.loads import crank_histories
from fea.mesh import mesh_step

NAME = "crankshaft"
SECTION = "crankshaft"
CANDIDATES = 1500          # surface nodes kept for the cycle-resolved fatigue scan (largest unit stresses)


def _ring(xc, zc, R):
    """gmsh MathEval distance from the circle of radius R about the x-parallel axis through (0, zc) at |x| = xc."""
    return f"Sqrt((Sqrt(y*y+(z-{zc})*(z-{zc}))-{R})^2+(Fabs(x)-{xc})^2)"


def build(p, fspec):
    return crank_throw_solid(p) & _halfspace(x=True)


def mesh_refine(p, fd, level):
    d = throw_dims(p)
    fs = float(fd["mesh"]["fillet_size"][level])          # CAD units (mm), like h_max
    ramp = p.L(fd["fillet_ramp"])
    r = d["crank_r"]
    pin = _ring(d["pw"], r, d["pin_r"])
    main = _ring(d["web_x0"], 0.0, d["main_r"])
    pin_sector = pin + f"+0.5*(Fabs(z-{r})+(z-{r}))"     # penalty grows on the far side of the pin axis
    main_sector = main + "+0.5*(Fabs(z)-z)"                # penalty grows below the crank axis
    out = []
    for expr, rf, size in ((pin_sector, d["r_pin_f"], fs), (main_sector, d["r_main_f"], fs),
                           (pin, d["r_pin_f"], fd["fillet_ring_factor"] * fs),
                           (main, d["r_main_f"], fd["fillet_ring_factor"] * fs)):
        out.append({"expr": expr, "size": size, "dist_min": rf, "dist_max": rf + ramp})
    return out


def analyse(p, fspec, level: int, hist: list | None = None, tag: str = "") -> LevelResult:
    fd = fspec[SECTION]
    hist = hist or crank_histories(p, fspec)
    d = throw_dims(p)
    k = p.k
    S, pw, r, pin_r = d["half_span"] / k, d["pw"] / k, d["crank_r"] / k, d["pin_r"] / k
    step = export(build(p, fspec), workdir(fspec, NAME) / f"{NAME}{tag}.step")
    mm = fd["mesh"]
    h = float(mm["h_max"][level])
    mesh = mesh_step(step, 1.0 / k, h, h * mm["h_min_ratio"], int(mm["curvature_pts"][level]),
                     refine=mesh_refine(p, fd, level))
    mat = material(p.engine, p.cad["materials"]["crankshaft"])
    nu = fspec["materials"]["poisson"][p.cad["materials"]["crankshaft"]]
    deck = ccx.Deck(mesh, ccx.Material("steel_crank", mat["E"], nu, mat["rho"]), solver=fd.get("solver", "SPOOLES"))
    # Supports without multi-point constraints (dense MPC rows make the direct solver very slow): the rear
    # journal end face (bearing centre plane) is held in y and z and free in x, i.e. it can rotate about
    # y / z (simple support) but not about x (torque reaction). The pin load is a full-circumference
    # cosine traction p0 cos(phi - phi_F) over the load band (resultant exact, sign-linear).
    face = plane_nodes(mesh, 0, S)
    mid = plane_nodes(mesh, 0, 0.0)
    band = fd["pin_load_band"] * pw
    _e, _f, _cn, _six, c, _n, _a = mesh.boundary_faces()
    phi = np.arctan2(c[:, 1], c[:, 2] - r)           # angle from +z towards +y about the pin axis
    rr = np.hypot(c[:, 1], c[:, 2] - r)
    band_sel = np.nonzero((np.abs(rr - pin_r) < 0.02 * pin_r) & (c[:, 0] < band))[0]
    half = 0.5
    dl_r, _ = cosine_pressure(mesh, band_sel, np.cos(phi[band_sel]), (0, 0, 1), half)
    dl_t, _ = cosine_pressure(mesh, band_sel, np.sin(phi[band_sel]), (0, 1, 0), half)
    fy = [(int(n), 2, 0.0) for n in face] + [(int(n), 3, 0.0) for n in face]
    sym = [(int(n), 1, 0.0) for n in mid] + fy
    twist = 1.0e-3                                       # rad, prescribed rigid twist of the rear face
    x_anchor = int(face[np.argmin(np.hypot(mesh.nodes[face, 1], mesh.nodes[face, 2]))])
    anti = [(int(n), 2, 0.0) for n in mid] + [(int(n), 3, 0.0) for n in mid] + [(x_anchor, 1, 0.0)]
    anti += [(int(n), 2, -twist * mesh.nodes[n, 2]) for n in face] + [(int(n), 3, twist * mesh.nodes[n, 1]) for n in face]
    deck.steps = [
        ccx.Step(dloads=dl_r, boundaries=sym),                                       # S1
        ccx.Step(dloads=dl_t, boundaries=sym),                                       # S2
        ccx.Step(boundaries=anti),                                                    # S3 (scaled below)
        ccx.Step(centrif=(1.0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), boundaries=sym),   # S4
    ]
    res = solve(deck, fspec, NAME, f"{NAME}{tag}_L{level}")
    N = mesh.n_nodes
    U = [res[i]["STRESS"][:N] for i in range(4)]
    RF = res[2]["FORC"][face]
    Mx = float(np.sum(mesh.nodes[face, 1] * RF[:, 2] - mesh.nodes[face, 2] * RF[:, 1]))
    U[2] = U[2] / Mx          # field of a unit moment +Mx applied at the rear face (reaction of tau = 1)
    # sign check of the through-torque case: torsional shear must be the same magnitude both halves by
    # construction; the rear reaction direction is fixed by statics (front receives the driving torque)
    # ROI: surface nodes away from the rigid face and the load band
    excl = fd["roi_exclusion"]
    bexcl = fd["band_exclusion"]
    surf = np.zeros(N, dtype=bool)
    surf[np.unique(mesh.boundary_faces()[3].ravel())] = True
    x = mesh.nodes[:, 0]
    rpin = np.hypot(mesh.nodes[:, 1], mesh.nodes[:, 2] - r)
    in_band = (x < band + bexcl) & (rpin > pin_r - bexcl)
    roi = surf & (x < S - excl) & ~in_band
    env = np.max([ccx.von_mises(u) for u in U[:3]], axis=0)
    idx = np.nonzero(roi)[0]
    idx = idx[np.argsort(env[idx])[::-1][:CANDIDATES]]
    Uc = np.stack([u[idx] for u in U])                # (4, n, 6)
    se, sut = mat["se"], mat["sut"]
    best = {"sf": math.inf}
    worst_vm = {"vm": 0.0}
    per_rpm = []
    for hp in hist:
        w2 = hp["omega"] ** 2
        for kth in range(hp["Fr"].shape[0]):
            Fr, Ft, Tin = hp["Fr"][kth], hp["Ft"][kth], hp["T_in"][kth]
            tau = Tin + Ft * r / 2
            base = Fr[:, None, None] * Uc[0] + Ft[:, None, None] * Uc[1] + w2 * Uc[3]
            tors = tau[:, None, None] * Uc[2]
            for side, sgn in (("rear", 1.0), ("front", -1.0)):
                Sh = base + sgn * tors                          # (n_theta, n, 6)
                svm = signed_von_mises(Sh)                      # (n_theta, n)
                sa, sm = goodman_from_history(svm.T)
                sf = goodman_sf(sa, sm, se, sut)
                j = int(np.argmin(sf))
                vm = np.abs(svm)
                it, jn = np.unravel_index(np.argmax(vm), vm.shape)
                if vm[it, jn] > worst_vm["vm"]:
                    worst_vm = {"vm": float(vm[it, jn]), "rpm": hp["rpm"], "throw": kth + 1, "side": side,
                                "theta_deg": float(np.degrees(hp["theta"][it])),
                                "at_m": mesh.nodes[idx[jn]].tolist(), "Fr_N": float(Fr[it]), "Ft_N": float(Ft[it]),
                                "T_in_Nm": float(Tin[it])}
                if sf[j] < best["sf"]:
                    best = {"sf": float(sf[j]), "rpm": hp["rpm"], "throw": kth + 1, "side": side,
                            "sigma_a_Pa": float(sa[j]), "sigma_m_Pa": float(sm[j]),
                            "at_m": mesh.nodes[idx[j]].tolist(), "node": int(idx[j])}
        per_rpm.append({"rpm": hp["rpm"], "max_Fr_abs_N": float(np.abs(hp["Fr"]).max()),
                        "max_T_in_Nm": float(np.abs(hp["T_in"]).max())})
    # FEA fillet Kt under the radial (bending) unit case vs a Peterson stepped-bar estimate
    pin_fillet = roi & (np.abs(x - (pw)) < 2 * d["r_pin_f"] / k) & (rpin < pin_r + 2 * d["r_pin_f"] / k)
    vm1 = ccx.von_mises(U[0])
    kt_fea_peak = float(vm1[pin_fillet].max())
    M = 0.5 * (S - pw)                                  # bending moment at the pin/web junction per unit force
    dpin = 2 * pin_r
    sig_nom = M / (math.pi * dpin**3 / 32)
    kt_fea = kt_fea_peak / sig_nom
    kt_pet = peterson_kt_bending(2 * d["web_pin_R"] / k, dpin, d["r_pin_f"] / k)
    # full fields at the worst state for plotting (rear half, its own throw/rpm)
    hp = [hh for hh in hist if hh["rpm"] == worst_vm["rpm"]][0]
    kt_ = worst_vm["throw"] - 1
    it = int(round(math.radians(worst_vm["theta_deg"]) / (hp["theta"][1] - hp["theta"][0])))
    Fr, Ft, Tin = hp["Fr"][kt_, it], hp["Ft"][kt_, it], hp["T_in"][kt_, it]
    sgn = 1.0 if worst_vm["side"] == "rear" else -1.0
    S_worst = Fr * U[0] + Ft * U[1] + sgn * (Tin + Ft * r / 2) * U[2] + hp["omega"] ** 2 * U[3]
    vm_worst = ccx.von_mises(S_worst)
    data = {"fatigue_critical": best, "worst_state": worst_vm, "se_Pa": se, "sut_Pa": sut,
            "kt": {"fea": kt_fea, "peterson": kt_pet, "ratio": kt_fea / kt_pet, "nominal_per_N_Pa": sig_nom,
                   "D_m": 2 * d["web_pin_R"] / k, "d_m": dpin, "r_m": d["r_pin_f"] / k},
            "per_rpm": per_rpm, "yield_like_sf_on_sut": sut / worst_vm["vm"]}
    return LevelResult(level, mesh.n_nodes, mesh.n_elems, worst_vm["vm"], {"fatigue": best["sf"]}, data, mesh,
                       {"vm_worst": vm_worst, "vm_unit_radial": vm1})


def plots(r: LevelResult, outdir, p):
    from fea.common import plot_field
    ws = r.data["worst_state"]
    plot_field(r.mesh, r.fields["vm_worst"] / 1.0e6, outdir / "crankshaft_vm_worst.png",
               f"Crank throw von Mises at the worst state: throw {ws['throw']} ({ws['side']} half), "
               f"{ws['rpm']:.0f} rpm, {ws['theta_deg']:.0f} deg", "von Mises [MPa]", view=(20, 60))
    return ["crankshaft_vm_worst.png"]
