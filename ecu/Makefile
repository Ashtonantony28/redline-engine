# Redline ECU (M8): host build of the portable C11 core.
#   make -C ecu test     build and run the C unit tests (no external dependencies)
#   make -C ecu lib      build/libecu_host.so for the Python HIL bridge (hil/bridge.py)
#   make -C ecu gen      regenerate include/ecu/ecu_calib_gen.h from the specs (needs the venv)
#   make -C ecu check-gen  fail if the committed generated header is stale
CC ?= gcc
PYTHON ?= /home/user/redline-venv/bin/python
BUILD := build
CSTD := -std=c11 -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion -Wno-sign-conversion -O2 -fPIC
INC := -Iinclude -Isrc
# Decoder, scheduler and knock control must not use floating point: prove it at compile time.
NOFLOAT := -mgeneral-regs-only

INT_SRC := src/decoder.c src/sched.c src/knock.c
FLT_SRC := src/ecu.c src/tables.c src/calib.c
HOST_SRC := host/hal_host.c
HDRS := $(wildcard include/ecu/*.h src/*.h)

INT_OBJ := $(patsubst %.c,$(BUILD)/%.o,$(INT_SRC))
FLT_OBJ := $(patsubst %.c,$(BUILD)/%.o,$(FLT_SRC) $(HOST_SRC))
CORE_OBJ := $(INT_OBJ) $(FLT_OBJ)

.PHONY: all lib test gen check-gen clean
all: lib

$(INT_OBJ): $(BUILD)/%.o: %.c $(HDRS)
	@mkdir -p $(dir $@)
	$(CC) $(CSTD) $(NOFLOAT) $(INC) -c $< -o $@

$(FLT_OBJ): $(BUILD)/%.o: %.c $(HDRS)
	@mkdir -p $(dir $@)
	$(CC) $(CSTD) $(INC) -c $< -o $@

lib: $(BUILD)/libecu_host.so
$(BUILD)/libecu_host.so: $(CORE_OBJ)
	$(CC) -shared -o $@ $^ -lm

$(BUILD)/ecu_tests: tests/test_ecu.c tests/tharness.h $(CORE_OBJ)
	$(CC) $(CSTD) -Wno-conversion -Wno-float-conversion $(INC) -Itests tests/test_ecu.c $(CORE_OBJ) -lm -o $@

test: $(BUILD)/ecu_tests lib
	./$(BUILD)/ecu_tests

gen:
	cd .. && $(PYTHON) ecu/tools/gen_calib.py

check-gen:
	cd .. && $(PYTHON) ecu/tools/gen_calib.py --check

clean:
	rm -rf $(BUILD)
