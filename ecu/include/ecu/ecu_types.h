/* Redline ECU: plain data types shared by the core, the HAL and the host bridge. C11, no heap. */
#ifndef ECU_TYPES_H
#define ECU_TYPES_H

#include <stdint.h>

typedef enum {
    ECU_EV_INJ = 1,     /* injector opens at tick, stays open arg ticks (electrical pulse width) */
    ECU_EV_DWELL = 2,   /* coil starts charging at tick */
    ECU_EV_SPARK = 3    /* coil current cut (spark) at tick */
} ecu_ev_type_t;

typedef struct {
    uint32_t tick;
    uint32_t arg;
    uint8_t type;
    uint8_t cyl;        /* 0-based cylinder number (cylinder 1 = 0) */
} ecu_event_t;

typedef enum {
    ECU_SYNC_NONE = 0,  /* searching for the missing-tooth gap */
    ECU_SYNC_CRANK = 1, /* crank angle known modulo 360 deg */
    ECU_SYNC_FULL = 2   /* cam phase known: cycle angle 0..720 deg */
} ecu_sync_t;

/* fault bits */
#define ECU_FAULT_MAP     (1u << 0)
#define ECU_FAULT_TPS     (1u << 1)
#define ECU_FAULT_LAMBDA  (1u << 2)
#define ECU_FAULT_IAT     (1u << 3)
#define ECU_FAULT_CLT     (1u << 4)
#define ECU_FAULT_PEDAL   (1u << 5)

#endif
