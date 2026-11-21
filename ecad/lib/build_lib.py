#!/usr/bin/env python3
"""Regenerate ecad/lib/{symbols,footprints} for the Redline ECU board (output is committed).

Sources
  DRONE  /home/user/drone-design/ecad/lib: parts reused from the NF flight controller
         (STM32H743VITx, LMR36015, LDL1117S33R, USB-C, Tag-Connect, BOOT0 switch, ...),
         copied verbatim with their provenance (see drone-design ecad/lib/build_lib.py)
  KSYM   sparse clone of https://gitlab.com/kicad/libraries/kicad-symbols.git
  KFP    sparse clone of https://gitlab.com/kicad/libraries/kicad-footprints.git

Kinds of symbol
  copied   verbatim from KiCad or from drone-design
  authored generated here from the pin table in PINOUT_EVIDENCE / PINS below

Every authored pinout was read from the vendor datasheet PDF (TI and Diodes Inc. sites were
reachable from the build container) except the Bosch CJ125, whose datasheet is not public:
its pin table comes from a secondary source and is marked VERIFY.

Usage: python3 lib/build_lib.py [--ksym DIR] [--kfp DIR] [--drone DIR]
"""
import argparse
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))

PINOUT_EVIDENCE = {
    "LM74700-Q1": "datasheet: TI SNOSD17G, Figure 5-1 / Table 5-1 (DBV SOT-23-6)",
    "TPS7B4250-Q1": "datasheet: TI SLVSCA0C, Pin Functions (DBV SOT-23-5)",
    "DRV8873H-Q1": "datasheet: TI SLVSDY7B, Figure 5-1 / Pin Functions (PWP HTSSOP-24, H variant)",
    "UCC27524A-Q1": "datasheet: TI SLVSCC1C, Figure 5-1 / Table 5-1 (D SOIC-8)",
    "TCAN1044AV-Q1": "datasheet: TI SLLSFJ3D, Figure 5-2 / Table 5-1 (D SOIC-8, V suffix = VIO)",
    "TPIC8101": "datasheet: TI SLIS110C, Pin Functions (DW SOIC-20)",
    "LM1815": "datasheet: TI SNOSBU8F, Figure 1 connection diagram (D SOIC-14)",
    "SN74AHCT125-Q1": "datasheet: TI SN74AHCT125-Q1, Figure 4-1 / Table 4-1 (PW TSSOP-14)",
    "TPD2E2U06-Q1": "datasheet: TI SLLSEJ9E, Figure 5-1 (DBZ SOT-23)",
    "ESD2CAN24-Q1": "datasheet: TI SLVSFW5D, Figure 4-1 (DBZ SOT-23)",
    "ZXMS6006DGQ": "datasheet: Diodes DS40468 rev 3, SOT223 Type DN pinout (1 IN, 2 D, 3 S, tab D)",
    "C_Elec": "generic polarized capacitor: pad 1 = +, pad 2 = - (KiCad CP_Elec footprints)",
    "CJ125": "VERIFY: Bosch CJ125 datasheet is not public; LQFP-32 pin table from the rusEFI "
    "hw_modular cj125_Module library (Cj125_LQPF32.lib); pins 1, 9, 16, 24, 25 are absent there "
    "and are treated as NC",
}

# Authored symbols: name -> (footprint, datasheet, description, left pins, right pins)
# pins: (number, name, KiCad electrical type)
PINS = {
    "LM74700-Q1": (
        "Package_TO_SOT_SMD:SOT-23-6",
        "https://www.ti.com/lit/ds/symlink/lm74700-q1.pdf",
        "Automotive ideal diode controller, 3.2-65 V, drives an external N-FET, SOT-23-6",
        [("6", "ANODE", "power_in"), ("3", "EN", "input"), ("2", "GND", "power_in")],
        [("4", "CATHODE", "input"), ("5", "GATE", "output"), ("1", "VCAP", "passive")],
    ),
    "TPS7B4250-Q1": (
        "Package_TO_SOT_SMD:SOT-23-5",
        "https://www.ti.com/lit/ds/symlink/tps7b4250-q1.pdf",
        "Automotive 50 mA voltage tracker for off-board sensor supplies, 40 V, SOT-23-5",
        [("3", "VIN", "power_in"), ("1", "ADJ_EN", "input"), ("2", "GND", "power_in"), ("5", "GND", "power_in")],
        [("4", "VOUT", "power_out")],
    ),
    "DRV8873H-Q1": (
        "Package_SO:HTSSOP-24-1EP_4.4x7.8mm_P0.65mm_EP3.2x5mm",
        "https://www.ti.com/lit/ds/symlink/drv8873-q1.pdf",
        "Automotive 40 V 10 A H-bridge motor driver, hardware interface, current mirror, HTSSOP-24",
        [("13", "VM", "power_in"), ("20", "VM", "power_in"), ("7", "EN_IN1", "input"), ("8", "PH_IN2", "input"),
         ("11", "nSLEEP", "input"), ("9", "DISABLE", "input"), ("3", "MODE", "input"), ("4", "SR", "input"),
         ("5", "nITRIP", "input"), ("6", "nOL", "input"), ("24", "GND", "power_in"), ("25", "GND", "power_in"),
         ("16", "SRC", "power_in"), ("17", "SRC", "power_in")],
        [("18", "OUT1", "output"), ("19", "OUT1", "output"), ("14", "OUT2", "output"), ("15", "OUT2", "output"),
         ("2", "nFAULT", "open_collector"), ("10", "IPROPI1", "output"), ("12", "IPROPI2", "output"),
         ("1", "DVDD", "power_out"), ("21", "VCP", "passive"), ("22", "CPH", "passive"), ("23", "CPL", "passive")],
    ),
    "UCC27524A-Q1": (
        "Package_SO:SOIC-8_3.9x4.9mm_P1.27mm",
        "https://www.ti.com/lit/ds/symlink/ucc27524a-q1.pdf",
        "Automotive dual 5 A low-side gate driver, TTL inputs, enables, SOIC-8",
        [("6", "VDD", "power_in"), ("2", "INA", "input"), ("4", "INB", "input"), ("1", "ENA", "input"),
         ("8", "ENB", "input"), ("3", "GND", "power_in")],
        [("7", "OUTA", "output"), ("5", "OUTB", "output")],
    ),
    "TCAN1044AV-Q1": (
        "Package_SO:SOIC-8_3.9x4.9mm_P1.27mm",
        "https://www.ti.com/lit/ds/symlink/tcan1044a-q1.pdf",
        "Automotive CAN FD transceiver with VIO level shifting, SOIC-8",
        [("3", "VCC", "power_in"), ("5", "VIO", "power_in"), ("1", "TXD", "input"), ("4", "RXD", "output"),
         ("8", "STB", "input"), ("2", "GND", "power_in")],
        [("7", "CANH", "bidirectional"), ("6", "CANL", "bidirectional")],
    ),
    "TPIC8101": (
        "Package_SO:SOIC-20W_7.5x12.8mm_P1.27mm",
        "https://www.ti.com/lit/ds/symlink/tpic8101.pdf",
        "Dual-channel knock sensor interface: band-pass, rectifier, integrator, SPI, SOIC-20",
        [("1", "VDD", "power_in"), ("20", "CH1P", "input"), ("19", "CH1N", "input"), ("18", "CH1FB", "output"),
         ("15", "CH2P", "input"), ("16", "CH2N", "input"), ("17", "CH2FB", "output"), ("3", "VREF", "output"),
         ("2", "GND", "power_in"), ("5", "NC", "no_connect"), ("6", "NC", "no_connect")],
        [("4", "OUT", "output"), ("7", "INT_HOLD", "input"), ("8", "CS", "input"), ("13", "SCLK", "input"),
         ("12", "SDI", "input"), ("11", "SDO", "tri_state"), ("9", "XIN", "input"), ("10", "XOUT", "output"),
         ("14", "TEST", "input")],
    ),
    "LM1815": (
        "Package_SO:SOIC-14_3.9x8.7mm_P1.27mm",
        "https://www.ti.com/lit/ds/symlink/lm1815.pdf",
        "Adaptive variable-reluctance sensor amplifier, 2-12 V, SOIC-14",
        [("8", "VCC", "power_in"), ("3", "VIN", "input"), ("5", "MODE", "input"), ("7", "PEAK_CAP", "passive"),
         ("14", "RC_TIMING", "passive"), ("2", "GND", "power_in"),
         ("1", "NC", "no_connect"), ("4", "NC", "no_connect"), ("6", "NC", "no_connect"), ("13", "NC", "no_connect")],
        [("10", "GATED_OUT", "output"), ("12", "REF_PULSE_OUT", "open_collector"), ("11", "INPUT_SELECT", "input"),
         ("9", "TIMING_IN", "input")],
    ),
    "SN74AHCT125-Q1": (
        "Package_SO:TSSOP-14_4.4x5mm_P0.65mm",
        "https://www.ti.com/lit/ds/symlink/sn74ahct125-q1.pdf",
        "Automotive quad 3-state buffer, 5 V supply, TTL-level inputs (3.3 V to 5 V shifting), TSSOP-14",
        [("14", "VCC", "power_in"), ("2", "1A", "input"), ("5", "2A", "input"), ("9", "3A", "input"),
         ("12", "4A", "input"), ("1", "1OE", "input"), ("4", "2OE", "input"), ("10", "3OE", "input"),
         ("13", "4OE", "input"), ("7", "GND", "power_in")],
        [("3", "1Y", "tri_state"), ("6", "2Y", "tri_state"), ("8", "3Y", "tri_state"), ("11", "4Y", "tri_state")],
    ),
    "TPD2E2U06-Q1": (
        "Package_TO_SOT_SMD:SOT-23",
        "https://www.ti.com/lit/ds/symlink/tpd2e2u06-q1.pdf",
        "Automotive 2-channel 6 V ESD protection, SOT-23",
        [("1", "IO1", "passive"), ("2", "IO2", "passive")],
        [("3", "GND", "passive")],
    ),
    "ESD2CAN24-Q1": (
        "Package_TO_SOT_SMD:SOT-23",
        "https://www.ti.com/lit/ds/symlink/esd2can24-q1.pdf",
        "Automotive +/-24 V ESD protection for CAN, SOT-23",
        [("1", "IO1", "passive"), ("2", "IO2", "passive")],
        [("3", "GND", "passive")],
    ),
    "ZXMS6006DGQ": (
        "Package_TO_SOT_SMD:SOT-223-3_TabPin2",
        "https://www.diodes.com/assets/Datasheets/ZXMS6006DGQ.pdf",
        "AEC-Q101 60 V self-protected low-side IntelliFET, active clamp, 3.3/5 V logic input, SOT-223",
        [("1", "IN", "input"), ("3", "S", "passive")],
        [("2", "D", "passive")],
    ),
    "C_Elec": (
        "Capacitor_SMD:CP_Elec_8x10",
        "",
        "Aluminium electrolytic capacitor (pins named for Zener, which connects by pin name)",
        [("1", "+", "passive")],
        [("2", "-", "passive")],
    ),
    "CJ125": (
        "Package_QFP:LQFP-32_7x7mm_P0.8mm",
        "https://www.bosch-semiconductors.com/",
        "Bosch wideband lambda (LSU 4.9) interface IC, SPI, LQFP-32",
        [("30", "UB", "power_in"), ("19", "VCC", "power_in"), ("5", "VCCS", "power_in"), ("18", "VCCSS", "power_in"),
         ("31", "UN", "passive"), ("32", "IP", "passive"), ("20", "VM", "passive"), ("2", "IA", "passive"),
         ("22", "UP", "passive"), ("21", "US", "passive"), ("8", "RS", "passive"), ("10", "RM", "passive"),
         ("11", "CM", "passive"), ("29", "GND", "power_in"), ("28", "GNDS", "power_in"),
         ("1", "NC", "no_connect"), ("9", "NC", "no_connect"), ("16", "NC", "no_connect"), ("24", "NC", "no_connect"),
         ("25", "NC", "no_connect")],
        [("23", "UA", "output"), ("12", "UR", "output"), ("27", "RF", "passive"), ("26", "CF", "passive"),
         ("3", "OSZ", "passive"), ("4", "DIAHG", "input"), ("6", "DIAHD", "input"), ("7", "RST", "input"),
         ("17", "SS", "input"), ("15", "SI", "input"), ("14", "SO", "tri_state"), ("13", "SCLK", "input")],
    ),
}

# Parts reused from drone-design ecad/lib (same files, provenance in that repo's build_lib.py)
DRONE_SYMBOLS = ["STM32H743VITx", "LMR36015", "LDL1117S33R", "USB_C_Receptacle_USB2.0_16P",
                 "Conn_ARM_SWD_TagConnect_TC2030-NL", "SW_Push", "D_Schottky_Small", "L_Small", "MountingHole"]
DRONE_FOOTPRINTS = ["LQFP-100_14x14mm_P0.5mm", "Texas_RNX0012C_VQFN-14-11-1EP_2x3mm_P0.5mm_EP0.25x1.825mm",
                    "SOT-223-3_TabPin2", "USB_C_Receptacle_HRO_TYPE-C-31-M-12",
                    "Tag-Connect_TC2030-IDC-NL_2x03_P1.27mm_Vertical", "SW_SPST_TS-1088-xR020",
                    "D_SOD-123", "L_Sunlord_SWPA4030S", "SOT-23-5", "SOT-23-6", "SOT-23"]

# Copied verbatim from KiCad
KICAD_SYMBOLS = {
    "Conn_02x12_Top_Bottom": [("Connector_Generic", "Conn_02x12_Top_Bottom")],
    "Q_NMOS_GDS": [("Transistor_FET", "Q_NMOS_GDS")],
}
KICAD_FOOTPRINTS = [
    "Package_TO_SOT_SMD.pretty/TO-252-2",
    "Connector_Molex.pretty/Molex_Micro-Fit_3.0_43045-2400_2x12_P3.00mm_Horizontal",
    "Package_QFP.pretty/LQFP-32_7x7mm_P0.8mm",
    "Package_SO.pretty/SOIC-20W_7.5x12.8mm_P1.27mm",
    "Package_SO.pretty/SOIC-14_3.9x8.7mm_P1.27mm",
    "Package_SO.pretty/SOIC-8_3.9x4.9mm_P1.27mm",
    "Package_SO.pretty/HTSSOP-24-1EP_4.4x7.8mm_P0.65mm_EP3.2x5mm",
    "Package_SO.pretty/TSSOP-14_4.4x5mm_P0.65mm",
    "Capacitor_SMD.pretty/CP_Elec_8x10",
    "Diode_SMD.pretty/D_PowerDI-123",
    "MountingHole.pretty/MountingHole_3.2mm_M3",
]


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def symbol_body(text):
    start = text.index('\t(symbol "')
    end = text.rstrip().rindex(")")
    return text[start:end].rstrip() + "\n"


HEADER = '(kicad_symbol_lib\n\t(version 20251024)\n\t(generator "kicad_symbol_editor")\n\t(generator_version "10.0")\n'


def lib_of(*bodies):
    return HEADER + "".join(bodies) + ")\n"


def pin_sexpr(ptype, x, y, rot, name, number, hide=False):
    font = "\t\t\t\t\t(effects\n\t\t\t\t\t\t(font\n\t\t\t\t\t\t\t(size 1.27 1.27)\n\t\t\t\t\t\t)\n\t\t\t\t\t)\n"
    return (
        "\t\t\t(pin %s line\n\t\t\t\t(at %.2f %.2f %d)\n\t\t\t\t(length 2.54)\n%s"
        '\t\t\t\t(name "%s"\n%s\t\t\t\t)\n\t\t\t\t(number "%s"\n%s\t\t\t\t)\n\t\t\t)\n'
        % (ptype, x, y, rot, "\t\t\t\t(hide yes)\n" if hide else "", name, font, number, font)
    )


def prop_sexpr(key, value, y, hide=True):
    return (
        '\t\t(property "%s" "%s"\n\t\t\t(at 0 %.2f 0)\n%s\t\t\t(effects\n\t\t\t\t(font\n\t\t\t\t\t(size 1.27 1.27)\n\t\t\t\t)\n\t\t\t)\n\t\t)\n'
        % (key, value.replace('"', '\\"'), y, "\t\t\t(hide yes)\n" if hide else "")
    )


def author(name, footprint, datasheet, descr, left, right, evidence):
    """Box symbol (same generator as drone-design). Pins sharing a name stack (extra ones hidden)."""
    rows = max(len({p[1] for p in left}), len({p[1] for p in right}))
    half_h = 2.54 * (rows + 1) / 2
    half_w = 12.7
    body = '\t(symbol "%s"\n\t\t(exclude_from_sim no)\n\t\t(in_bom yes)\n\t\t(on_board yes)\n' % name
    body += prop_sexpr("Reference", "U", half_h + 2.54, hide=False)
    body += prop_sexpr("Value", name, half_h + 5.08, hide=False)
    body += prop_sexpr("Footprint", footprint, -half_h - 2.54)
    body += prop_sexpr("Datasheet", datasheet, -half_h - 5.08)
    body += prop_sexpr("Description", descr + " [pinout: " + evidence + "]", -half_h - 7.62)
    body += (
        '\t\t(symbol "%s_0_1"\n\t\t\t(rectangle\n\t\t\t\t(start %.2f %.2f)\n\t\t\t\t(end %.2f %.2f)\n'
        "\t\t\t\t(stroke\n\t\t\t\t\t(width 0.254)\n\t\t\t\t\t(type default)\n\t\t\t\t)\n"
        "\t\t\t\t(fill\n\t\t\t\t\t(type background)\n\t\t\t\t)\n\t\t\t)\n\t\t)\n" % (name, -half_w, half_h, half_w, -half_h)
    )
    body += '\t\t(symbol "%s_1_1"\n' % name
    for side, x, rot in ((left, -half_w - 2.54, 0), (right, half_w + 2.54, 180)):
        seen = {}
        for number, pname, ptype in side:
            hide = pname in seen
            if not hide:
                seen[pname] = half_h - 2.54 * (len(seen) + 1)
            body += pin_sexpr("passive" if hide and ptype != "no_connect" else ptype, x, seen[pname], rot, pname, number, hide)
    body += "\t\t)\n\t\t(embedded_fonts no)\n\t)\n"
    return lib_of(body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ksym", default="/home/user/kicadlib/gsym")
    ap.add_argument("--kfp", default="/home/user/kicadlib/gfp")
    ap.add_argument("--drone", default="/home/user/drone-design/ecad/lib")
    a = ap.parse_args()
    sym_out = os.path.join(HERE, "symbols")
    fp_out = os.path.join(HERE, "footprints")
    os.makedirs(sym_out, exist_ok=True)
    os.makedirs(fp_out, exist_ok=True)

    for name in DRONE_SYMBOLS:
        shutil.copyfile(os.path.join(a.drone, "symbols", name + ".kicad_sym"), os.path.join(sym_out, name + ".kicad_sym"))
    for name in DRONE_FOOTPRINTS:
        shutil.copyfile(os.path.join(a.drone, "footprints", name + ".kicad_mod"), os.path.join(fp_out, name + ".kicad_mod"))

    for out, parts in KICAD_SYMBOLS.items():
        bodies = [symbol_body(read(os.path.join(a.ksym, lib + ".kicad_symdir", n + ".kicad_sym"))) for lib, n in parts]
        write(os.path.join(sym_out, out + ".kicad_sym"), lib_of(*bodies))
    for rel in KICAD_FOOTPRINTS:
        src = os.path.join(a.kfp, rel + ".kicad_mod")
        shutil.copyfile(src, os.path.join(fp_out, os.path.basename(src)))

    for name, (fp, ds, descr, left, right) in PINS.items():
        write(os.path.join(sym_out, name + ".kicad_sym"), author(name, fp, ds, descr, left, right, PINOUT_EVIDENCE[name]))

    write(
        os.path.join(HERE, "LICENSE-KiCad.txt"),
        "Symbols and footprints here are copied or derived from the KiCad libraries\n"
        "(https://gitlab.com/kicad/libraries), licensed CC-BY-SA 4.0 with the KiCad library\n"
        "exception, directly or via drone-design/ecad/lib. See build_lib.py for the source of each file.\n",
    )


if __name__ == "__main__":
    main()
