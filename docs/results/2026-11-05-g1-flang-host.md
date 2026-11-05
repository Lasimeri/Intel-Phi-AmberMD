# G1 step 1: Flang from the patched LLVM 23.1.1, 2026-11-05

`scripts/g1-flang-build.sh`. Source: the shared `llvm-project` tree at the SERIES commit
6dfe1677ab8d, patches 0001-0017 applied (48 modified files, `git diff` sha256 4f18bfdf...; plus 14
new untracked files such as `X86InstrKNC.td` and `kncintrin.h`, which `git diff` does not cover).
Projects clang, lld, mlir, flang; X86 only; Release; default triple x86_64-unknown-linux-musl.

| Item | Value |
|---|---|
| Wall time, 8 jobs (Ryzen 5800X) | 22:53.33 |
| Largest process RSS | 3,112,272 kB |
| Lowest host MemAvailable (10 s samples) | 9,255,520 kB |
| Build tree | 3.25 GB at 84 % of the steps |
| Result | `flang version 23.1.1`, target x86_64-unknown-linux-musl, 0 errors |

## Same compiler as the KNC clang in use
12 C test kernels of Intel-Phi-LLVM (`tests/*.c`, `tests/s3/*.c`) compiled with `-march=knc -O2`
by this build's clang and by the reference KNC clang: `.text` identical in 12 of 12, 27,410 bytes,
775 instructions using zmm registers. Control: one changed byte is detected.
(A first attempt compared empty files because `llvm-objcopy` cannot write to `/dev/stdout`; that
false "12 of 12" was discarded.)

## First Fortran compile for KNC
`flang -march=knc -O2 -S` on a print program: 0 lines with `xmm`, `movsd`, `movss`, `cvtsi2sd` or
`mxcsr`; object built. Not linked yet: flang-rt for the card does not exist yet.
