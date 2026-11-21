/* Redline ECU on Vortex-V (M9): the HIL mailbox, a struct in SoC RAM shared by the firmware
 * (hal_rv64.c) and the Python plant (hil/rv64.py pokes and peeks it through the Verilator model).
 *
 * Vortex-V has no interrupts or CSRs yet, so capture "interrupts" are delivered as mailbox
 * commands that the firmware's main loop polls:
 *   plant: write arg[] / in, then cmd != HIL_CMD_NONE      (capture event or slow-loop tick)
 *   core:  run the handler, fill events[] / out / cycles, then cmd = HIL_CMD_NONE   (acknowledge)
 * The plant clocks the core until the acknowledge; plant time stands still meanwhile (lockstep).
 * `cycles` is the handler's own cost, read from the SoC CYCLE register around the call.
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
    HIL_CMD_EXIT = 10         /* store PASS to tohost: ends the simulation */
};

enum hil_layout {
    HIL_L_CMD = 0, HIL_L_ARG, HIL_L_N_EVENTS, HIL_L_OVERFLOW, HIL_L_CYCLES, HIL_L_EVENTS,
    HIL_L_EVENT_SIZE, HIL_L_EVENT_TICK, HIL_L_EVENT_ARG, HIL_L_EVENT_TYPE, HIL_L_EVENT_CYL,
    HIL_L_IN, HIL_L_IN_SIZE, HIL_L_OUT, HIL_L_OUT_SIZE, HIL_L_STATUS, HIL_L_CONSTS,
    HIL_L_KNOCK_X, HIL_L_KNOCK_N_IN, HIL_L_KNOCK_LOGIT, HIL_L_N
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
} hil_mbox_t;

#endif
