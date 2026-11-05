/* x87pc53.c: optional startup hook for card binaries. Sets the x87 precision-control field to
 * 53-bit mantissa (control word 0x027F instead of Linux's 0x037F), so scalar double arithmetic
 * rounds like SSE double on the host. The exponent range stays extended. Linked in only when a
 * run must be compared bit for bit with a host build. */
__attribute__((constructor)) static void x87_pc53(void) {
    unsigned short cw;
    __asm__ __volatile__("fnstcw %0" : "=m"(cw));
    cw = (unsigned short)((cw & ~0x0300u) | 0x0200u);
    __asm__ __volatile__("fldcw %0" : : "m"(cw));
}
