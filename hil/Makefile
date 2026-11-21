# Redline M9 / M9.1: co-simulation of the ECU firmware on the Vortex-V RTL (Verilator) against
# the Python MVEM plant. All results are MODELLED.
#   make -C hil cosim     build firmware + model + host reference, then
#                         (frozen time) bit-exact event stream / slow-loop outputs vs the host build,
#                         cycle budget, the int8 knock sensor fixture on the RTL, and
#                         (real time, ADR 0005) latency, angle error, CPU load, missed teeth at
#                         7000 rpm and in the 2000 rpm tip-in, M9 firmware vs M9.1 firmware
#   make -C hil rt-sweep  the real-time runs at 25 / 50 / 100 MHz and the minimum-clock bisection
#   make -C hil model     only build build/libvortex_soc.so
# Vortex-V is consumed read-only from $(VORTEX_V_ROOT) (default in ecu/targets/rv64/soc.mk).
RV64 := ../ecu/targets/rv64
include $(RV64)/soc.mk
PYTHON ?= /home/user/redline-venv/bin/python
VERILATOR ?= verilator
OBJ := build/obj_dir
JOBS ?= 4
CXX_OPT ?= -O2

# rtl/core/files.f + rtl/soc/files.f, in Vortex-V's order (paths are relative to its root)
filelist = $(addprefix $(VORTEX_V_ROOT)/,$(shell grep -v '^\s*//' $(VORTEX_V_ROOT)/$(1) | grep -v '^\s*$$'))
RTL := $(call filelist,rtl/core/files.f) $(call filelist,rtl/soc/files.f)
GPARAMS := $(addprefix -G,$(SOC_PARAMS)) -GRAM_BYTES=$(SOC_RAM_BYTES)

.PHONY: cosim rt rt-sweep model firmware ref
cosim: model firmware ref
	cd .. && $(PYTHON) -m hil.cosim
	cd .. && $(PYTHON) -m hil.realtime

rt: model firmware ref
	cd .. && $(PYTHON) -m hil.realtime

rt-sweep: model firmware ref
	cd .. && $(PYTHON) -m hil.realtime --sweep

model: build/libvortex_soc.so
build/libvortex_soc.so: cosim/vortex_sim.cpp cosim/public.vlt $(RTL) $(RV64)/soc.mk
	@mkdir -p $(OBJ)
	$(VERILATOR) --cc --top-module vortex_soc -Wno-fatal -O3 --x-assign fast --x-initial fast \
	  --noassert $(GPARAMS) -CFLAGS "-fPIC $(CXX_OPT)" --Mdir $(OBJ) cosim/public.vlt $(RTL) --build -j $(JOBS) \
	  > build/verilator.log 2>&1 || (tail -30 build/verilator.log; exit 1)
	VROOT=$$($(VERILATOR) --getenv VERILATOR_ROOT); \
	$(CXX) -shared -fPIC $(CXX_OPT) -std=c++17 -I$(OBJ) -I$$VROOT/include -I$$VROOT/include/vltstd \
	  cosim/vortex_sim.cpp $(OBJ)/Vvortex_soc__ALL.a $(OBJ)/libverilated.a -pthread -o $@

firmware:
	$(MAKE) -s -C $(RV64) all

ref:
	$(MAKE) -s -C $(RV64) ref
