# G1 gate: Fortran runs on the 7120A, 2026-11-05

## Toolchain pieces
- Flang 23.1.1 from the patched tree (results/2026-11-05-g1-flang-host.md).
- flang-rt for the card: `scripts/g1-flang-rt-knc.sh`, standalone runtimes build, `-march=knc`, static,
  against the 7120A stack's musl and libc++ 22.1.8 sysroot: `libflang_rt.runtime.a` 10,143,584 B, 120
  of 120 steps; 15 intrinsic module files (`iso_c_binding.mod` and others).
- Link pieces the Flang driver expects but the stack does not have under those names: compiler-rt
  builtins and `crtbegin`/`crtend` come from the stack's own clang 22.1.8 resource directory
  (`lib/clang/22/lib/linux`), copied into this Flang's resource directory by the same script.
- Link line: `flang --target=x86_64-unknown-linux-musl -march=knc --sysroot=<7120A sysroot> -O2 -static
  -fuse-ld=lld -rtlib=compiler-rt ... -lc++ -lc++abi -lunwind`. Flang has no `-unwindlib` option.

Traps met: without the modules every `use iso_c_binding` fails ("Source file 'iso_c_binding.mod' was
not found"); without the builtins lld asks for `-lgcc`, `-lgcc_eh`, `crtend.o`.

## First program
`print *, "hello from fortran", sqrt(2d0), erfc(0.5d0)`: 1,181,424 B static binary, phi-isa-audit
0 illegal, 0 suspect; on the card `hello from fortran 1.4142135623730951 .4795001221869535`, rc 0.

## Conformance suite: `scripts/g1-conformance.sh`, `tests/fortran/conf.f90` + `conf_c.c`
Card binary 1,371,600 B, audit 0 illegal, 0 suspect. **12 of 12 checks PASS on the card**, also on the
host (gfortran 16.2.1): namelist (sander's mdin read), format-write, format-read, inpcrd-format,
stream-io, allocatable (3 x 49,910 array), derived-type, complex-to-c, complex-from-c, array-to-c,
integer8, strings, math.

Bit comparison of 11 values printed in hex, card against host: **9 of 11 identical**. Identical:
sqrt, exp, log, erfc(0.5), erfc(2.3), sin, cos, atan2, pow. Different: the two sums of 1/j for
j = 1..100,000.

| Sum | Card | Host | Error against a float128 reference |
|---|---|---|---|
| explicit loop | 40282E27A22F3FB0 | 40282E27A22F3F7C | card 0, host -34.4 (units of ref x 2^-52) |
| `sum()` intrinsic | 40282E27A22F3FB0 | 40282E27A22F3F7C | same |

The card keeps the running sum in an 80-bit x87 register, so it is the more accurate one here.

### x87 precision control at 53 bits (`tests/fortran/x87pc53.c`)
A constructor sets the x87 control word's precision field to 53 bits (0x027F for 0x037F). Card
result with it: **10 of 11 identical**, all checks PASS. The loop sum now equals the host's bits.
The `sum()` intrinsic still differs: flang-rt's SUM uses Kahan summation (`flang-rt/lib/runtime/
sum.cpp` line 12, 55), gfortran adds in order. That difference is the algorithm, not the hardware.

Consequence for G2: card results can be compared with the host bit for bit where the code is the same
algorithm, if the pc53 hook is linked. Without it the card is more precise and not bit-identical.
Whether pc53 costs speed is not measured.
