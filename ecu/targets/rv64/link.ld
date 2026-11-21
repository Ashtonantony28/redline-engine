/* Redline ECU on the Vortex-V SoC (M9). One RAM at RAM_BASE, RAM_BYTES long (--defsym from the
 * Makefile, the same value the Verilator model is built with). vortex_soc ends the program on a
 * store to TOHOST_ADDR = RAM_BASE + 0x1000 (M1 contract, used by Vortex-V's crt0.S), so that
 * doubleword is reserved between the start-up code and the rest of .text. */
OUTPUT_ARCH(riscv)
ENTRY(_start)

RAM_BASE = 0x80000000;
TOHOST_OFFSET = 0x1000;
STACK_ALIGN = 16;

SECTIONS
{
  . = RAM_BASE;
  .text : {
    *(.text.init)
    . = TOHOST_OFFSET;
    QUAD(0)
    *(.text .text.*)
  }
  .rodata : ALIGN(8) { *(.rodata .rodata.* .srodata .srodata.*) }
  .data : ALIGN(8) { *(.data .data.* .sdata .sdata.*) }
  .bss (NOLOAD) : ALIGN(8) {
    _bss_start = .;
    *(.bss .bss.* .sbss .sbss.* COMMON)
    . = ALIGN(8);
    _bss_end = .;
  }
  _image_end = .;
  _stack_top = (RAM_BASE + RAM_BYTES) & ~(STACK_ALIGN - 1);
  ASSERT(_image_end + STACK_MIN_BYTES <= _stack_top, "ECU image and stack do not fit SOC_RAM_BYTES")
  /DISCARD/ : { *(.comment) *(.riscv.attributes) *(.note*) *(.eh_frame*) }
}
