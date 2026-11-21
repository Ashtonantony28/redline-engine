# Shared by ecu/targets/rv64/Makefile (firmware link) and hil/Makefile (Verilator model): the
# Vortex-V checkout (consumed read-only) and the SoC build that the ECU runs on.
VORTEX_V_ROOT ?= /home/user/vortex-v
# RAM size [bytes, power of 2]; vortex_soc's default 16 KiB cannot hold soft-float + the ECU.
SOC_RAM_BYTES ?= 65536
# Core configuration = Vortex-V's board configuration "vortex_ooo_core_sync" (verif/core/cores.py,
# what the Arty A7 top builds): 2-wide speculative OoO, LSQ, pipelined dual Xmlacc, block RAM.
SOC_PARAMS ?= CORE_OOO=1 WIDTH=2 SPEC_BRANCHES=1 LSQ=1 XM_LATENCY=3 SYNC_MEM=1 XM_PORTS=2
