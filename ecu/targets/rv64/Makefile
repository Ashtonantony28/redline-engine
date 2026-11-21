# Redline ECU, RV64 target for the Vortex-V SoC (M9). The portable core (ecu/src) is compiled
# unchanged for rv64im / lp64 (no F/D: every float in the slow loop is libgcc soft-float).
#   make -C ecu/targets/rv64            build/ecu_rv64.elf + .hex (64-bit words, vortex_ram format)
#   make -C ecu/targets/rv64 ref        build/libecu_host_ref.so: the host build with the same
#                                       portable libm (ecu_libm.c), the bit-exact reference
#   make -C ecu/targets/rv64 nofloat    prove the sync-path units contain no float (see below)
include soc.mk
PYTHON ?= /home/user/redline-venv/bin/python
CROSS ?= riscv64-unknown-elf-
CC_RV := $(CROSS)gcc
HOSTCC ?= gcc
ECU := ../..
BUILD := build
GEN := gen
STACK_MIN_BYTES ?= 16384

ARCH := -march=rv64im -mabi=lp64 -mcmodel=medany -mstrict-align
CSTD := -std=c11 -Wall -Wextra -Werror -Wshadow -O2 -ffunction-sections
# include/ first: <math.h> resolves to the portable-libm shim for the core's slow loop
INC := -Iinclude -I$(ECU)/include -I$(ECU)/src -I. -I$(GEN) -I$(VORTEX_V_ROOT)/sw/common
# no libc: keep GCC from turning byte loops (including our own memset) into memset/memcpy calls
RVFLAGS := $(ARCH) $(CSTD) -fno-tree-loop-distribute-patterns -nostdlib -nostartfiles -static $(INC)

# Sync path (crank decoder, angle scheduler, knock control) and the int8 sensor kernel: integer only.
INT_SRC := $(ECU)/src/decoder.c $(ECU)/src/sched.c $(ECU)/src/knock.c knock_kernel.c
FLT_SRC := $(ECU)/src/ecu.c $(ECU)/src/tables.c $(ECU)/src/calib.c ecu_libm.c hal_rv64.c
INT_OBJ := $(patsubst %.c,$(BUILD)/rv/%.o,$(notdir $(INT_SRC)))
FLT_OBJ := $(patsubst %.c,$(BUILD)/rv/%.o,$(notdir $(FLT_SRC)))
HDRS := $(wildcard $(ECU)/include/ecu/*.h $(ECU)/src/*.h *.h include/*.h) $(GEN)/knock_weights.h Makefile soc.mk
vpath %.c $(ECU)/src .

.PHONY: all ref nofloat clean
all: $(BUILD)/ecu_rv64.hex nofloat

$(GEN)/knock_weights.h: ../../../ml/policies/knock_int8.json ../../../ml/knock/export.py
	cd ../../.. && $(PYTHON) -m ml.knock.export

$(BUILD)/rv/%.o: %.c $(HDRS)
	@mkdir -p $(dir $@)
	$(CC_RV) $(RVFLAGS) -c $< -o $@

$(BUILD)/rv/crt0.o: $(VORTEX_V_ROOT)/sw/common/crt0.S
	@mkdir -p $(dir $@)
	$(CC_RV) $(ARCH) -c $< -o $@

$(BUILD)/ecu_rv64.elf: $(BUILD)/rv/crt0.o $(INT_OBJ) $(FLT_OBJ) link.ld soc.mk Makefile
	$(CC_RV) $(RVFLAGS) -Wl,--gc-sections -Wl,--no-warn-rwx-segments -Tlink.ld \
	  -Wl,--defsym=RAM_BYTES=$(SOC_RAM_BYTES) -Wl,--defsym=STACK_MIN_BYTES=$(STACK_MIN_BYTES) \
	  $(BUILD)/rv/crt0.o $(INT_OBJ) $(FLT_OBJ) -lgcc -o $@

$(BUILD)/ecu_rv64.bin: $(BUILD)/ecu_rv64.elf
	$(CROSS)objcopy -O binary --gap-fill 0 $< $@

# one 64-bit little-endian word per line from RAM_BASE (vortex_ram INIT_FILE format)
$(BUILD)/ecu_rv64.hex: $(BUILD)/ecu_rv64.bin
	od -An -v -tx8 -w8 $< | tr -d ' ' > $@

# No float in the sync path, RV64 flavour of the host build's -mgeneral-regs-only proof: with
# rv64im/lp64 any float operation becomes a libgcc soft-float call, so the integer-only objects
# must reference no undefined symbol at all beyond the core's own API, and in particular none of
# libgcc's __*sf* / __*df* / __*tf* routines.
nofloat: $(INT_OBJ)
	@for o in $(INT_OBJ); do \
	  if $(CROSS)nm -u $$o | grep -E '__[a-z0-9]*(sf|df|tf)[a-z0-9]*' ; then \
	    echo "FLOAT in sync-path unit $$o"; exit 1; fi; done
	@echo "nofloat: $(notdir $(INT_OBJ)) reference no soft-float routine"

# Host reference: same core sources + host HAL, but libm = ecu_libm.c (via include/math.h).
REF_SRC := $(ECU)/src/decoder.c $(ECU)/src/sched.c $(ECU)/src/knock.c $(ECU)/src/ecu.c \
           $(ECU)/src/tables.c $(ECU)/src/calib.c $(ECU)/host/hal_host.c ecu_libm.c
ref: $(BUILD)/libecu_host_ref.so
$(BUILD)/libecu_host_ref.so: $(REF_SRC) $(HDRS)
	@mkdir -p $(BUILD)
	$(HOSTCC) -std=c11 -Wall -Wextra -Werror -O2 -fPIC -ffp-contract=off -shared \
	  -Iinclude -I$(ECU)/include -I$(ECU)/src $(REF_SRC) -o $@

clean:
	rm -rf $(BUILD)
