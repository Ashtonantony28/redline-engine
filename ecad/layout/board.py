"""Board geometry and net classes shared by the layout scripts (pure Python, no pcbnew).

Coordinates are millimetres from the board centre, x to the right, y down, looking at the top.
On the KiCad sheet the board centre sits at ORIGIN.
"""
ORIGIN = (100.0, 100.0)
W, H = 120.0, 90.0          # outline (boards/redline_ecu/redline_ecu.zen docstring)
CORNER_R = 3.0
HOLE_XY = (55.5, 40.5)      # M3 holes, one per corner
HOLE_KEEPOUT_R = 3.6        # copper keep-out around the M3 holes (screw head 5.5 mm + washer)
EDGE_INSET = 0.5            # pours stop this far inside the outline (rule: copper to edge 0.4)

# The connector edge: both Micro-Fit headers face out of the bottom edge.
CONN_Y = H / 2              # courtyard front of J1 / J2 on the edge

# Nets and widths (the router and the finish step read these)
CLR = 0.15                  # Default and Power class clearance (redline_ecu.zen)
TRACK = 0.2                 # Default class track
NARROW_POWER = 0.3          # supply nets, and high-current nets where they neck down at a pin
WIDE = 1.0                  # high-current nets away from their pins
# supply and actuator nets routed 0.3 mm (the Q map)
POWER_NETS = ("3V3", "5V", "VDDA", "5V_S1", "5V_S2", "VBUS", "IGN1", "IGN2", "IGN3", "IGN4",
              "OUT.IGN1_DRV", "OUT.IGN2_DRV", "OUT.IGN3_DRV", "OUT.IGN4_DRV", "PUMP", "KL.HTR_GATE",
              "PWR.BUCK_SW", "PWR.BUCK_BOOT", "LSU_IP", "LSU_VM", "LSU_UN", "LSU_IA", "KL.CJ_UB")
# nets that carry amps: 1.0 mm tracks plus a pour on their own layer (route.py finish)
WIDE_NETS = ("VBAT_IN", "VBAT", "ETB_M_P", "ETB_M_N", "INJ1", "INJ2", "INJ3", "INJ4", "VVT",
             "LSU_HTR_N", "PWR.RB_SRC")
PLANE_NETS = ("GND",)       # In1.Cu is one solid ground plane


def corner_holes():
    hx, hy = HOLE_XY
    return [(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)]


def inside_outline(x, y, margin):
    """True if the point is at least `margin` inside the rounded-rectangle outline."""
    hx, hy = W / 2 - margin, H / 2 - margin
    if abs(x) > hx or abs(y) > hy:
        return False
    kx, ky = W / 2 - CORNER_R, H / 2 - CORNER_R
    ax, ay = abs(x), abs(y)
    if ax > kx and ay > ky:
        return ((ax - kx) ** 2 + (ay - ky) ** 2) ** 0.5 <= CORNER_R - margin
    return True
