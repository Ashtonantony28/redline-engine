/* Redline ECU on Vortex-V (M9, M9.1): the HIL mailbox, a struct in SoC RAM shared by the firmware
 * (hal_rv64.c) and the Python plant (hil/rv64.py pokes and peeks it through the Verilator model).
 *
 * Frozen-time command mode (M9, kept for the bit-exactness tests): every ECU entry point is a
 * mailbox command that the firmware's main loop polls:
 *   plant: write arg[] / in, then cmd != HIL_CMD_NONE      (capture event or slow-loop tick)
 *   core:  run the handler, fill events[] / out / cycles, then cmd = HIL_CMD_NONE   (acknowledge)
 * The plant clocks the core until the acknowledge; plant time stands still meanwhile (lockstep).
 * `cycles` is the handler's own cost, read from the SoC CYCLE register around the call.
 *
 * Real-time mode (M9.1, ADR 0005; HIL_CMD_RT, the `rt` block): plant time is core cycles. The
 * Vortex-V SoC has no capture/compare timer, so the co-simulation model plays that peripheral
 * through `rt`, cycle-accurately (hil/cosim/vortex_sim.cpp):
 *   rt.timer      the free-running ECU timer [ticks of ECU_TIMER_HZ], advanced by the model
 *   rt.cap[]      input-capture FIFO: the plant stores an edge {kind, captured tick} and bumps
 *                 cap_wr at the cycle the edge occurs; the firmware bumps cap_rd when it has
 *                 serviced an entry. A full FIFO drops the edge (a missed tooth).
 *   events[]      output compares: hal_emit_event stores the event and bumps rt.evt_wr; the model
 *                 stamps the cycle. An event armed after its own tick fires late, at once.
 *   rt.in / out   sensor samples and slow-loop outputs, double-buffered (in_sel / out_sel name the
 *                 valid copy, so neither side reads a half-written struct); slow_seq counts steps.
 * The firmware returns to command mode when the plant sets rt.stop.
 *
 * hil/rv64.py parses the enums below, and reads every field offset from layout[] at run time, so
 * nothing here is retyped in Python. */
#ifndef HIL_MBOX_H
#define HIL_MBOX_H

#include <stdint.h>

#include "ecu/ecu.h"

#define HIL_MBOX_MAGIC 0x52444C39u      /* "RDL9": firmware booted, layout[] valid */
#define HIL_MBOX_EVENTS 64u             /* event slots per command (hil/bridge.py MAX_DRAIN) */
#define HIL_NSTATUS 16u
#define HIL_NCONSTS 11u
#define HIL_NARGS 3u
#define HIL_KNOCK_MAX_IN 256u           /* int8 samples per knock window, upper bound */

enum hil_cmd {
    HIL_CMD_NONE = 0,
    HIL_CMD_RESET = 1,
    HIL_CMD_CRANK = 2,        /* arg0 = captured timer tick */
    HIL_CMD_CAM = 3,          /* arg0 = captured timer tick */
    HIL_CMD_KNOCK = 4,        /* arg0 = cylinder, arg1 = integrated counts */
    HIL_CMD_SLOW = 5,         /* in -> out */
    HIL_CMD_SCHED_SET = 6,    /* arg0 = cylinder, arg1 = advance [angle units], arg2 = pulse [ticks] */
    HIL_CMD_STATUS = 7,       /* -> status[] */
    HIL_CMD_CONSTS = 8,       /* -> consts[] */
    HIL_CMD_KNOCK_INFER = 9,  /* knock_x[] -> knock_logit (int8 MLP on Xmlacc) */
    HIL_CMD_EXIT = 10,        /* store PASS to tohost: ends the simulation */
    HIL_CMD_RT = 11           /* run in real time from rt.cfg[] until rt.stop (acknowledged at the end) */
};

/* ---- real-time mode (M9.1) ---- */
#define HIL_CAP_SLOTS 32u               /* input-capture FIFO depth */

enum hil_cap_kind {
    HIL_CAP_CRANK = 1,        /* a0 = captured timer tick */
    HIL_CAP_CAM = 2,          /* a0 = captured timer tick */
    HIL_CAP_KNOCK = 3         /* a0 = cylinder, a1 = integrated counts */
};

enum hil_rt_cfg {
    HIL_RT_CLK_HZ = 0,        /* core clock [Hz]: converts timer ticks to CYCLE counts */
    HIL_RT_FLAGS,             /* HIL_RT_F_* */
    HIL_RT_POLL_CYCLES,       /* timer period while no edge can be predicted (search, stall, late edge) */
    HIL_RT_LEAD_CYCLES,       /* wake this long before the predicted edge (covers interrupt entry) */
    HIL_RT_SPIN_CYCLES,       /* longest wait inside the handler for a predicted edge */
    HIL_RT_NCFG
};

#define HIL_RT_F_IRQ 1u                 /* captures serviced from the machine timer interrupt; without
                                         * it the background loop polls the FIFO (the M9 structure) */
#define HIL_RT_F_PREDICT 2u             /* arm the timer for the next predicted crank edge */

enum hil_rt_stat {
    HIL_RT_S_IRQS = 0,        /* timer interrupts taken */
    HIL_RT_S_ISR_CYCLES,      /* cycles inside irq_handler, sum / max (entry and exit code not included) */
    HIL_RT_S_ISR_MAX,
    HIL_RT_S_SERVICED,        /* capture entries serviced */
    HIL_RT_S_SPIN_TIMEOUTS,   /* predicted edges that did not arrive within the spin window */
    HIL_RT_S_SLOW_STEPS,
    HIL_RT_S_SLOW_WALL_SUM,   /* slow-loop step, snapshot to commit, including preemption [cycles] */
    HIL_RT_S_SLOW_WALL_MAX,
    HIL_RT_S_SLOW_OVERRUNS,   /* slow-loop periods skipped because the previous step was still running */
    HIL_RT_S_SNAP_RETRIES,    /* snapshots repeated because a capture was serviced meanwhile */
    HIL_RT_S_STALLS,          /* stall verdicts committed to the live decoder */
    HIL_RT_S_IDLE_ITERS,      /* background-loop passes with nothing to do */
    HIL_RT_S_IDLE_CAL_ITERS,  /* calibration: this many idle passes ... */
    HIL_RT_S_IDLE_CAL_CYCLES, /* ... take this many cycles undisturbed */
    HIL_RT_S_RUN_CYCLES,      /* cycles from go to stop */
    HIL_RT_NSTATS
};

typedef struct {
    uint32_t kind;
    uint32_t a0;
    uint32_t a1;
    uint32_t reserved;
} hil_cap_t;

typedef struct {
    uint32_t cfg[HIL_RT_NCFG];
    uint32_t ready;           /* firmware: calibrated, waiting for go */
    uint32_t go;              /* plant: time 0 is now */
    uint32_t stop;            /* plant: leave real-time mode */
    uint32_t timer;           /* model: ECU timer [ticks] */
    uint32_t cap_wr;          /* plant: captures posted (free-running) */
    uint32_t cap_rd;          /* firmware: captures serviced (free-running) */
    uint32_t evt_wr;          /* firmware: events armed (free-running, slot = evt_wr % HIL_MBOX_EVENTS) */
    uint32_t in_sel;          /* plant: valid copy of in[] */
    uint32_t out_sel;         /* firmware: valid copy of out[] */
    uint32_t slow_seq;        /* firmware: slow-loop steps committed */
    hil_cap_t cap[HIL_CAP_SLOTS];
    ecu_inputs_t in[2];
    ecu_outputs_t out[2];
    uint64_t stats[HIL_RT_NSTATS];
} hil_rt_t;

enum hil_layout {
    HIL_L_CMD = 0, HIL_L_ARG, HIL_L_N_EVENTS, HIL_L_OVERFLOW, HIL_L_CYCLES, HIL_L_EVENTS,
    HIL_L_EVENT_SIZE, HIL_L_EVENT_TICK, HIL_L_EVENT_ARG, HIL_L_EVENT_TYPE, HIL_L_EVENT_CYL,
    HIL_L_IN, HIL_L_IN_SIZE, HIL_L_OUT, HIL_L_OUT_SIZE, HIL_L_STATUS, HIL_L_CONSTS,
    HIL_L_KNOCK_X, HIL_L_KNOCK_N_IN, HIL_L_KNOCK_LOGIT,
    HIL_L_RT_CFG, HIL_L_RT_READY, HIL_L_RT_GO, HIL_L_RT_STOP, HIL_L_RT_TIMER, HIL_L_RT_CAP_WR,
    HIL_L_RT_CAP_RD, HIL_L_RT_EVT_WR, HIL_L_RT_IN_SEL, HIL_L_RT_OUT_SEL, HIL_L_RT_SLOW_SEQ,
    HIL_L_RT_CAP, HIL_L_RT_CAP_SIZE, HIL_L_RT_IN, HIL_L_RT_OUT, HIL_L_RT_STATS, HIL_L_N
};

typedef struct {
    uint32_t magic;
    uint32_t layout[HIL_L_N];
    uint32_t cmd;
    uint32_t arg[HIL_NARGS];
    uint32_t n_events;
    uint32_t overflow;
    uint64_t cycles;
    ecu_event_t events[HIL_MBOX_EVENTS];
    ecu_inputs_t in;
    ecu_outputs_t out;
    uint32_t status[HIL_NSTATUS];
    uint32_t consts[HIL_NCONSTS];
    _Alignas(8) int8_t knock_x[HIL_KNOCK_MAX_IN];   /* read as 8-lane words: no misaligned loads */
    int32_t knock_logit;
    hil_rt_t rt;
} hil_mbox_t;

#endif
