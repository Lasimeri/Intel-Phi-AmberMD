#!/bin/bash
# G1 step 1: build Flang (with MLIR, clang, lld) from the same patched LLVM 23.1.1 tree the
# KNC clang comes from. The source tree is only read. Build and install go to the host disk.
# Memory is sampled every 10 s into logs/g1/mem.tsv (MemAvailable, SwapFree in kB).
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
SRC=${PHI_LLVM_SRC:-/media/1TB-NVMe/Xeon-Phi/intel-phi-llvm-build/llvm-project}
SERIES=${PHI_LLVM_SERIES:-"/media/1TB-NVMe/Xeon-Phi/Intel Phi LLVM/patches/SERIES"}
OUT=${G1_OUT:-$HOME/.cache/intel-phi-ambermd/llvm-flang}
JOBS=${JOBS:-8}
LOG=$R/logs/g1
mkdir -p "$OUT" "$LOG"

want=$(awk '/^# LLVM commit:/ {print $4}' "$SERIES")
have=$(git -C "$SRC" rev-parse HEAD)
[ "$want" = "$have" ] || { echo "g1: $SRC is at $have, SERIES pins $want; refusing" >&2; exit 1; }
git -C "$SRC" diff | sha256sum | awk '{print "patched tree diff sha256", $1}' | tee "$LOG/source.txt"
git -C "$SRC" diff --stat | tail -1 | tee -a "$LOG/source.txt"

( while :; do awk -v t="$(date +%s)" '/^MemAvailable:/{a=$2} /^SwapFree:/{s=$2} END{print t "\t" a "\t" s}' /proc/meminfo; sleep 10; done ) > "$LOG/mem.tsv" &
sampler=$!
trap 'kill $sampler 2>/dev/null' EXIT

cmake -S "$SRC/llvm" -B "$OUT/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$OUT/install" \
    -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DLLVM_USE_LINKER=lld \
    -DLLVM_ENABLE_PROJECTS="clang;lld;mlir;flang" -DLLVM_TARGETS_TO_BUILD=X86 \
    -DLLVM_DEFAULT_TARGET_TRIPLE=x86_64-unknown-linux-musl \
    -DLLVM_ENABLE_ASSERTIONS=OFF -DLLVM_INCLUDE_TESTS=OFF -DLLVM_INCLUDE_BENCHMARKS=OFF \
    -DLLVM_INCLUDE_EXAMPLES=OFF -DLLVM_INCLUDE_DOCS=OFF -DLLVM_ENABLE_BINDINGS=OFF \
    -DMLIR_INCLUDE_TESTS=OFF -DFLANG_INCLUDE_TESTS=OFF -DFLANG_INCLUDE_DOCS=OFF \
    -DFLANG_ENABLE_FLANG_RT=OFF \
    -DLLVM_ENABLE_ZLIB=ON -DLLVM_ENABLE_ZSTD=OFF -DLLVM_ENABLE_LIBXML2=OFF -DLLVM_ENABLE_TERMINFO=OFF \
    -DLLVM_PARALLEL_LINK_JOBS=2 -DLLVM_PARALLEL_COMPILE_JOBS="$JOBS" \
    > "$LOG/cmake.log" 2>&1
/usr/bin/time -v cmake --build "$OUT/build" -j "$JOBS" --target flang clang lld llvm-ar llvm-ranlib llvm-nm llvm-objdump \
    > "$LOG/build.log" 2>&1
echo "g1: flang built: $("$OUT/build/bin/flang" --version | head -1)"
awk 'NR==1||$2<m{m=$2} END{print "lowest MemAvailable", m, "kB"}' "$LOG/mem.tsv"
