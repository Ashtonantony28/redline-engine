// Verilator model of the Vortex-V SoC (vortex_soc, unmodified) as a shared library for the Python
// HIL bridge (hil/rv64.py, ctypes). The plant talks to the firmware through a mailbox in SoC RAM:
// it pokes a command, clocks the core until the firmware clears it, and peeks the results.
// One clock cycle = one posedge + one negedge; `cycles` counts them since reset release.
//
// M9.1 adds the real-time side (vx_rt_*): plant time is core cycles. While it is enabled the model
// plays the ECU's capture/compare timer peripheral, cycle by cycle:
//   * the free-running timer register (a RAM word) is advanced at clk_hz / timer_hz cycles per tick,
//   * every event the firmware arms (event-ring write index moves) is logged with the cycle at which
//     the store reached RAM, i.e. when a real output compare would have been armed,
//   * every capture the firmware finishes servicing (read index moves), every slow-loop commit and
//     every interrupt the core takes is logged with its cycle.
// The plant posts capture edges by poking the FIFO at the cycle the edge occurs (vx_rt_run_to, then
// vx_write). A retire-PC profiler (vx_prof_*) attributes cycles to the next retiring instruction.
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Vvortex_soc.h"
#include "Vvortex_soc___024root.h"
#include "verilated.h"

namespace {
const int RESET_CYCLES = 4;
const unsigned WORD_BYTES = 8;
const unsigned BITS_PER_BYTE = 8;
const QData BYTE_MASK = 0xFF;

const unsigned RT_EVENT_BYTES_MAX = 32;

struct RtEvent {                // one armed output compare
    uint64_t arm_cycle;         // cycles since the real-time epoch when the event reached RAM
    uint32_t cap_rd;            // captures fully serviced before it (the one in service = cap_rd)
    uint32_t n_raw;
    uint8_t raw[RT_EVENT_BYTES_MAX];   // the ecu_event_t as the firmware stored it
};

struct Rt {
    bool on = false;
    uint64_t epoch = 0;         // `cycles` at plant time 0
    uint64_t clk_hz = 1, timer_hz = 1, acc = 0;
    uint32_t tick = 0;
    uint64_t off_timer = 0, off_cap_rd = 0, off_evt_wr = 0, off_events = 0, off_slow_seq = 0;
    uint32_t ev_size = 0, ev_slots = 1;
    uint32_t cap_rd = 0, evt_wr = 0, slow_seq = 0;
    std::vector<RtEvent> events;
    std::vector<uint64_t> service, slow, irq;   // cycle stamps (since epoch)
};

struct Sim {
    std::unique_ptr<VerilatedContext> ctx;
    std::unique_ptr<Vvortex_soc> top;
    uint64_t cycles = 0;
    Rt rt;
    bool prof_on = false;
    uint64_t prof_pending = 0;
    std::unordered_map<uint64_t, uint64_t> prof;

    Sim() : ctx(new VerilatedContext), top(new Vvortex_soc(ctx.get())) {
        top->clk = 0;
        top->rst_n = 0;
        top->eval();                      // runs the RAM's initial block (zero fill)
    }
    QData &word(uint64_t i) { return top->rootp->vortex_soc__DOT__u_ram__DOT__mem[i]; }
    uint64_t words() const { return sizeof(top->rootp->vortex_soc__DOT__u_ram__DOT__mem) / sizeof(QData); }
    void tick() {
        top->clk = 1;
        top->eval();
        top->clk = 0;
        top->eval();
    }
    uint32_t rd32(uint64_t off) {
        return uint32_t(word(off / WORD_BYTES) >> (BITS_PER_BYTE * (off % WORD_BYTES)));
    }
    void wr32(uint64_t off, uint32_t v) {
        QData &w = word(off / WORD_BYTES);
        const unsigned sh = BITS_PER_BYTE * (off % WORD_BYTES);
        w = (w & ~(QData(0xFFFFFFFFu) << sh)) | (QData(v) << sh);
    }
    // One clock cycle with the profiler and the real-time peripheral model.
    void step() {
        tick();
        cycles++;
        auto *r = top->rootp;
        if (prof_on) {
            prof_pending++;
            if (r->vortex_soc__DOT__g_ooo__DOT__u_core__DOT__dbg_retire_valid) {
                prof[r->vortex_soc__DOT__g_ooo__DOT__u_core__DOT__dbg_retire_pc] += prof_pending;
                prof_pending = 0;
            }
        }
        if (!rt.on) return;
        const uint64_t now = cycles - rt.epoch;
        if (r->vortex_soc__DOT__g_ooo__DOT__u_core__DOT__dbg_irq) rt.irq.push_back(now);
        const uint32_t cap_rd = rd32(rt.off_cap_rd);
        const uint32_t evt_wr = rd32(rt.off_evt_wr);
        while (rt.evt_wr != evt_wr) {            // events are stored before the index (program order)
            RtEvent e;
            e.arm_cycle = now;
            e.cap_rd = rt.cap_rd;
            e.n_raw = rt.ev_size;
            const uint64_t base = rt.off_events + uint64_t(rt.evt_wr % rt.ev_slots) * rt.ev_size;
            for (uint32_t i = 0; i < rt.ev_size; ++i)
                e.raw[i] = uint8_t(word((base + i) / WORD_BYTES) >> (BITS_PER_BYTE * ((base + i) % WORD_BYTES)));
            rt.events.push_back(e);
            rt.evt_wr++;
        }
        while (rt.cap_rd != cap_rd) {
            rt.service.push_back(now);
            rt.cap_rd++;
        }
        const uint32_t slow_seq = rd32(rt.off_slow_seq);
        while (rt.slow_seq != slow_seq) {
            rt.slow.push_back(now);
            rt.slow_seq++;
        }
        rt.acc += rt.timer_hz;                   // the timer the core reads in the next cycle
        if (rt.acc >= rt.clk_hz) {
            do {
                rt.acc -= rt.clk_hz;
                rt.tick++;
            } while (rt.acc >= rt.clk_hz);
            wr32(rt.off_timer, rt.tick);
        }
    }
};

template <typename T>
uint64_t drain(std::vector<T> &v, T *dst, uint64_t max) {
    const uint64_t n = v.size() < max ? v.size() : max;
    if (n) std::memcpy(dst, v.data(), n * sizeof(T));
    v.erase(v.begin(), v.begin() + n);
    return n;
}
}  // namespace

extern "C" {

void *vx_new() { return new Sim; }
void vx_delete(void *h) { delete static_cast<Sim *>(h); }
uint64_t vx_ram_bytes(void *h) { return static_cast<Sim *>(h)->words() * WORD_BYTES; }
uint64_t vx_cycles(void *h) { return static_cast<Sim *>(h)->cycles; }
int vx_done(void *h) { return static_cast<Sim *>(h)->top->done; }
int vx_pass(void *h) { return static_cast<Sim *>(h)->top->pass; }
int vx_leds(void *h) { return static_cast<Sim *>(h)->top->leds; }

// Byte access to the RAM (offset from RAM_BASE). Returns 0 on success.
int vx_write(void *h, uint64_t off, const uint8_t *src, uint64_t n) {
    Sim *s = static_cast<Sim *>(h);
    if (off + n > s->words() * WORD_BYTES) return -1;
    for (uint64_t i = 0; i < n; ++i) {
        QData &w = s->word((off + i) / WORD_BYTES);
        unsigned sh = BITS_PER_BYTE * ((off + i) % WORD_BYTES);
        w = (w & ~(BYTE_MASK << sh)) | (QData(src[i]) << sh);
    }
    return 0;
}

int vx_read(void *h, uint64_t off, uint8_t *dst, uint64_t n) {
    Sim *s = static_cast<Sim *>(h);
    if (off + n > s->words() * WORD_BYTES) return -1;
    for (uint64_t i = 0; i < n; ++i)
        dst[i] = uint8_t(s->word((off + i) / WORD_BYTES) >> (BITS_PER_BYTE * ((off + i) % WORD_BYTES)));
    return 0;
}

// Hold reset for a few cycles and release it (the image must already be in RAM).
void vx_reset(void *h) {
    Sim *s = static_cast<Sim *>(h);
    s->top->rst_n = 0;
    for (int i = 0; i < RESET_CYCLES; ++i) s->tick();
    s->top->rst_n = 1;
    s->cycles = 0;
}

// Clock until the 32-bit RAM word at byte offset `off` equals `value` (or the program ends).
// Returns the cycles spent, or -1 if `max_cycles` passed first.
int64_t vx_run_until(void *h, uint64_t off, uint32_t value, uint64_t max_cycles) {
    Sim *s = static_cast<Sim *>(h);
    const uint64_t wi = off / WORD_BYTES;
    const unsigned sh = BITS_PER_BYTE * (off % WORD_BYTES);
    for (uint64_t n = 0; n < max_cycles; ++n) {
        if (uint32_t(s->word(wi) >> sh) == value || s->top->done) return int64_t(n);
        s->step();
    }
    return -1;
}

// ---- real time (M9.1) ----
// Plant time 0 = now. Offsets are RAM byte offsets of the firmware's real-time mailbox fields.
void vx_rt_begin(void *h, const uint64_t *off, uint32_t ev_size, uint32_t ev_slots, uint64_t clk_hz,
                 uint64_t timer_hz) {
    Sim *s = static_cast<Sim *>(h);
    Rt r;
    r.on = true;
    r.epoch = s->cycles;
    r.clk_hz = clk_hz;
    r.timer_hz = timer_hz;
    r.off_timer = off[0];
    r.off_cap_rd = off[1];
    r.off_evt_wr = off[2];
    r.off_events = off[3];
    r.off_slow_seq = off[4];
    r.ev_size = ev_size > RT_EVENT_BYTES_MAX ? RT_EVENT_BYTES_MAX : ev_size;
    r.ev_slots = ev_slots;
    r.cap_rd = s->rd32(r.off_cap_rd);
    r.evt_wr = s->rd32(r.off_evt_wr);
    r.slow_seq = s->rd32(r.off_slow_seq);
    s->rt = r;
    s->wr32(r.off_timer, 0);
}

void vx_rt_end(void *h) { static_cast<Sim *>(h)->rt.on = false; }
uint64_t vx_rt_now(void *h) {
    Sim *s = static_cast<Sim *>(h);
    return s->cycles - s->rt.epoch;
}

// Clock until `cycle` cycles after the epoch. Returns 0, or 1 if the program ended first.
int vx_rt_run_to(void *h, uint64_t cycle) {
    Sim *s = static_cast<Sim *>(h);
    const uint64_t target = s->rt.epoch + cycle;
    while (s->cycles < target) {
        if (s->top->done) return 1;
        s->step();
    }
    return 0;
}

uint64_t vx_rt_events(void *h, RtEvent *dst, uint64_t max) { return drain(static_cast<Sim *>(h)->rt.events, dst, max); }
uint64_t vx_rt_services(void *h, uint64_t *dst, uint64_t max) { return drain(static_cast<Sim *>(h)->rt.service, dst, max); }
uint64_t vx_rt_slows(void *h, uint64_t *dst, uint64_t max) { return drain(static_cast<Sim *>(h)->rt.slow, dst, max); }
uint64_t vx_rt_irqs(void *h, uint64_t *dst, uint64_t max) { return drain(static_cast<Sim *>(h)->rt.irq, dst, max); }

// ---- retire-PC profiler ----
void vx_prof(void *h, int on) {
    Sim *s = static_cast<Sim *>(h);
    s->prof_on = on != 0;
    if (on) {
        s->prof.clear();
        s->prof_pending = 0;
    }
}

// Copies up to `max` (pc, cycles) pairs; returns the number of distinct PCs.
uint64_t vx_prof_dump(void *h, uint64_t *pcs, uint64_t *counts, uint64_t max) {
    Sim *s = static_cast<Sim *>(h);
    uint64_t i = 0;
    for (const auto &kv : s->prof) {
        if (i < max) {
            pcs[i] = kv.first;
            counts[i] = kv.second;
        }
        ++i;
    }
    return i;
}

}  // extern "C"
