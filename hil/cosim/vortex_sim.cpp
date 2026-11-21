// Verilator model of the Vortex-V SoC (vortex_soc, unmodified) as a shared library for the Python
// HIL bridge (hil/rv64.py, ctypes). The plant talks to the firmware through a mailbox in SoC RAM:
// it pokes a command, clocks the core until the firmware clears it, and peeks the results.
// One clock cycle = one posedge + one negedge; `cycles` counts them since reset release.
#include <cstdint>
#include <cstring>
#include <memory>

#include "Vvortex_soc.h"
#include "Vvortex_soc___024root.h"
#include "verilated.h"

namespace {
const int RESET_CYCLES = 4;
const unsigned WORD_BYTES = 8;
const unsigned BITS_PER_BYTE = 8;
const QData BYTE_MASK = 0xFF;

struct Sim {
    std::unique_ptr<VerilatedContext> ctx;
    std::unique_ptr<Vvortex_soc> top;
    uint64_t cycles = 0;

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
};
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
        s->tick();
        s->cycles++;
    }
    return -1;
}

}  // extern "C"
