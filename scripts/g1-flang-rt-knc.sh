#!/bin/bash
# G1 step 2: build flang-rt (the Fortran runtime) for the card: x86_64-unknown-linux-musl,
# -march=knc, static, against the 7120A stack's musl + libc++ sysroot. Standalone runtimes build
# with this repository's Flang (g1-flang-build.sh). The LLVM source tree is only read.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
SRC=${PHI_LLVM_SRC:-/media/1TB-NVMe/Xeon-Phi/intel-phi-llvm-build/llvm-project}
OUT=${G1_OUT:-$HOME/.cache/intel-phi-ambermd/llvm-flang}
SYSROOT=${PHI_SYSROOT_7120A:-$HOME/.cache/intel-phi-7120a-build/toolchain/sysroot}
JOBS=${JOBS:-12}
BIN=$OUT/build/bin
BUILD=$OUT/flang-rt-knc
LOG=$R/logs/g1
mkdir -p "$BUILD" "$LOG"
cflags="-march=knc -fno-jump-tables -falign-loops=1"
cmake -S "$SRC/runtimes" -B "$BUILD" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$OUT/knc" \
    -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_SYSTEM_PROCESSOR=x86_64 \
    -DCMAKE_SYSROOT="$SYSROOT" -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY \
    -DCMAKE_C_COMPILER="$BIN/clang" -DCMAKE_CXX_COMPILER="$BIN/clang++" \
    -DCMAKE_Fortran_COMPILER="$BIN/flang" -DCMAKE_Fortran_COMPILER_WORKS=yes \
    -DCMAKE_C_COMPILER_TARGET=x86_64-unknown-linux-musl -DCMAKE_CXX_COMPILER_TARGET=x86_64-unknown-linux-musl \
    -DCMAKE_Fortran_COMPILER_TARGET=x86_64-unknown-linux-musl \
    -DCMAKE_C_FLAGS="$cflags" -DCMAKE_CXX_FLAGS="$cflags -stdlib=libc++" -DCMAKE_Fortran_FLAGS="-march=knc" \
    -DCMAKE_AR="$BIN/llvm-ar" -DCMAKE_RANLIB="$BIN/llvm-ranlib" -DCMAKE_NM="$BIN/llvm-nm" \
    -DLLVM_BINARY_DIR="$OUT/build" -DLLVM_ENABLE_RUNTIMES=flang-rt \
    -DLLVM_DEFAULT_TARGET_TRIPLE=x86_64-unknown-linux-musl \
    -DFLANG_RT_ENABLE_STATIC=ON -DFLANG_RT_ENABLE_SHARED=OFF -DFLANG_RT_INCLUDE_TESTS=OFF \
    > "$LOG/flang-rt-cmake.log" 2>&1
cmake --build "$BUILD" -j "$JOBS" > "$LOG/flang-rt-build.log" 2>&1
find "$BUILD" -name 'libflang_rt*.a' -printf '%s %p\n'

# Make the driver find everything for the card target without extra flags: flang-rt and its
# intrinsic modules from this build, compiler-rt builtins and crtbegin/crtend from the 7120A
# stack's own clang (it never uses libgcc).
STACK_RT=${PHI_STACK_CLANG_RT:-/media/1TB-NVMe/Xeon-Phi/intel-phi-7120a-build/toolchain/llvm/lib/clang/22/lib/linux}
RD=$("$BIN/flang" -print-resource-dir)
T=x86_64-unknown-linux-musl
mkdir -p "$RD/lib/$T" "$RD/finclude/flang/$T"
cp "$BUILD/lib/clang/23/lib/$T/libflang_rt.runtime.a" "$RD/lib/$T/"
cp "$BUILD/lib/clang/23/finclude/flang/$T/"*.mod "$RD/finclude/flang/$T/"
cp "$STACK_RT/libclang_rt.builtins-x86_64.a" "$RD/lib/$T/libclang_rt.builtins.a"
cp "$STACK_RT/clang_rt.crtbegin-x86_64.o" "$RD/lib/$T/clang_rt.crtbegin.o"
cp "$STACK_RT/clang_rt.crtend-x86_64.o" "$RD/lib/$T/clang_rt.crtend.o"
ls "$RD/finclude/flang/$T" | wc -l | awk '{print "modules installed:", $1}'
