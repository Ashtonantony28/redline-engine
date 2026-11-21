#!/bin/sh
# Layout pipeline for redline_ecu (adapted from drone-design ecad/layout/run.sh).
# Run from anywhere; works in ecad/.  Usage: layout/run.sh <step>...
#   place   footprint dump -> autoplace.py (host) -> place.py (KiCad): outline, parts, In1 plane, keep-outs
#   dump    geometry JSON for plot.py / check_layout.py   -> layout/build/board.json
#   plot    review pictures                               -> layout/build/view_{top,inner}.png
#   pre     scripted copper: a ground via at every SMD ground pad (preroute.py)
#   auto    maze router for everything else (autoroute.py, pure Python + numpy)
#   finish  F/B ground pours, pours on the high-current nets, dangling vias out, silk clean-up, refill
#   drc     kicad-cli DRC                                 -> layout/build/drc.json (+ summary)
#   fab     gerbers, drill, pick-and-place, STEP          -> fab/
#   all     place pre auto finish dump drc plot
# KiCad 10 runs from the kicad/kicad:10.0.8 image; set KICAD_RUN to another wrapper
# (or to "env" when KiCad 10 is installed on this machine).
set -e
cd "$(dirname "$0")/.."
BOARD=boards/redline_ecu/layout/redline_ecu/redline_ecu.kicad_pcb
OUT=layout/build
PY=${PYTHON:-python3}
mkdir -p "$OUT"
kicad() {
    if [ -n "$KICAD_RUN" ]; then $KICAD_RUN "$@"; else
        docker run --rm -u "$(id -u)" -e HOME=/tmp -v "$PWD":"$PWD" -w "$PWD" kicad/kicad:10.0.8 "$@"
    fi
}
steps="$*"
[ "$steps" = all ] && steps="place pre auto finish dump drc plot"
for step in $steps; do
    case "$step" in
    place) rm -f "$OUT/drc.json"
        kicad python3 layout/dump.py "$BOARD" "$OUT/raw.json"
        $PY layout/autoplace.py "$OUT/raw.json" "$OUT/placement.json"
        kicad python3 layout/place.py "$BOARD" "$OUT/placement.json" ;;
    dump) kicad python3 layout/dump.py "$BOARD" "$OUT/board.json" ;;
    plot) $PY layout/plot.py "$OUT/board.json" "$OUT/view" $( [ -f "$OUT/drc.json" ] && echo "$OUT/drc.json" ) ;;
    pre) kicad python3 layout/dump.py "$BOARD" "$OUT/placed.json" &&
        $PY layout/preroute.py "$OUT/placed.json" "$OUT/preroute.json" &&
        kicad python3 layout/route.py apply "$BOARD" "$OUT/preroute.json" locked &&
        kicad python3 layout/dump.py "$BOARD" "$OUT/pre.json" && cp "$BOARD" "$OUT/prerouted.kicad_pcb" ;;
    auto) $PY layout/autoroute.py "$OUT/pre.json" "$OUT/autoroute.json" && cp "$OUT/prerouted.kicad_pcb" "$BOARD" &&
        kicad python3 layout/route.py apply "$BOARD" "$OUT/autoroute.json" unlocked ;;
    finish) kicad python3 layout/route.py finish "$BOARD" ;;
    drc) kicad kicad-cli pcb drc --format json --severity-all --units mm -o "$OUT/drc.json" "$BOARD" >/dev/null &&
        $PY layout/drc_summary.py "$OUT/drc.json" ;;
    fab) FAB=fab; rm -rf "$FAB"; mkdir -p "$FAB/gerbers"
        kicad kicad-cli pcb export gerbers -o "$FAB/gerbers/" --no-protel-ext --subtract-soldermask \
            --layers F.Cu,In1.Cu,In2.Cu,B.Cu,F.Paste,B.Paste,F.SilkS,B.SilkS,F.Mask,B.Mask,Edge.Cuts "$BOARD"
        kicad kicad-cli pcb export drill -o "$FAB/gerbers/" --format excellon --excellon-units mm --excellon-separate-th \
            --generate-map --map-format gerberx2 "$BOARD"
        kicad kicad-cli pcb export pos -o "$FAB/redline_ecu-pos.csv" --format csv --units mm --side both --exclude-dnp "$BOARD"
        kicad kicad-cli pcb export step -o "$FAB/redline_ecu.step" --force --subst-models --no-dnp "$BOARD" ;;
    *) echo "unknown step $step" >&2; exit 2 ;;
    esac
done
