#!/bin/bash
# md-test.sh [build|host|card] [m steps threads print]: the MD monitoring test.
#   build  the test simulation (tests/md/ljmd.f90 + src/card/team.c + src/card/mdstat.c) for the
#          card (Flang + KNC clang, static, audited) and for the host (gfortran + gcc)
#   host   run it on the host against a scratch file in place of the shared region, stop it
#          through the stop word after 1 s, and dump the block (logic check)
#   card   run it on the card in the background and record phitop frames every 5 s until it
#          ends: logs/md/<stamp>/ljmd.log and phitop.txt. Stop early with scripts/md-stop.sh.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
O=${G1_OUT:-$HOME/.cache/intel-phi-ambermd/llvm-flang}
BIN=$O/build/bin
SYSROOT=${PHI_SYSROOT_7120A:-$HOME/.cache/intel-phi-7120a-build/toolchain/sysroot}
: "${PHI_TREE:=/media/1TB-NVMe/Xeon-Phi/Intel Phi 7120A}"
: "${PHI_CARD:=0}"
: "${PHITOP:=phitop}"
out=$R/out/md
mkdir -p "$out"
mode=${1:-build}; shift || true
M=${1:-23} STEPS=${2:-2000} THREADS=${3:-120} PRINT=${4:-100} FRAME=${5:-1}
cflags=(-std=gnu11 -O2 -Wall -Wextra -Werror)
knc=(--target=x86_64-unknown-linux-musl -march=knc --sysroot="$SYSROOT" -O2)

build() {
    for c in team mdstat mdframe; do
        gcc "${cflags[@]}" -c "$R/src/card/$c.c" -o "$out/$c.host.o"
        "$BIN/clang" "${knc[@]}" "${cflags[@]}" -fno-jump-tables -c "$R/src/card/$c.c" -o "$out/$c.knc.o"
    done
    gfortran -O2 -Wall -Werror -J "$out" "$R/tests/md/ljmd.f90" "$out/team.host.o" "$out/mdstat.host.o" "$out/mdframe.host.o" -o "$out/ljmd.host" -pthread
    "$BIN/flang" "${knc[@]}" -static -fuse-ld=lld -rtlib=compiler-rt -module-dir "$out" "$R/tests/md/ljmd.f90" \
        "$out/team.knc.o" "$out/mdstat.knc.o" "$out/mdframe.knc.o" -o "$out/ljmd.knc" -lc++ -lc++abi -lunwind
    local line
    line=$("$PHI_TREE/host/target/debug/phi-isa-audit" "$out/ljmd.knc" | /usr/bin/grep -a '^result:')
    echo "audit: $line"
    [[ $line == "result: 0 illegal, 0 suspect,"* ]] || { echo "md-test: audit failed" >&2; exit 1; }
}

case "$mode" in
build) build ;;
host)
    f=$out/mdstat-host.bin
    dd if=/dev/zero of="$f" bs=1M count=2 status=none
    AMBER_MDSTAT_FILE=$f AMBER_MDSTAT_OFFSET=0 AMBER_MDFRAME_OFFSET=65536 "$out/ljmd.host" 10 100000 4 1000 &
    pid=$!
    sleep 1
    printf 'STOP' | dd of="$f" bs=1 seek=128 conv=notrunc status=none
    wait "$pid"
    od -A d -t x4 -N 16 "$f"
    od -A d -t d8 -j 16 -N 72 "$f"
    od -A d -t u4 -j 88 -N 8 "$f"
    ;;
card)
    log=$R/logs/md/$(date +%Y%m%d-%H%M%S)
    mkdir -p "$log"
    phi -c "$PHI_CARD" run mkdir -p /tmp/md < /dev/null
    phi -c "$PHI_CARD" put "$out/ljmd.knc" /tmp/md/ljmd < /dev/null
    phi -c "$PHI_CARD" run sh -c "cd /tmp/md && ./ljmd $M $STEPS $THREADS $PRINT $FRAME" < /dev/null > "$log/ljmd.log" 2>&1 &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
        { date +%T; timeout 20 "$PHITOP" -c "$PHI_CARD" --batch 2 -i 2 | tail -32; } >> "$log/phitop.txt" 2>&1 || true
        sleep 3
    done
    wait "$pid" || true
    echo "logs in $log"
    tail -5 "$log/ljmd.log"
    ;;
*) echo "usage: md-test.sh [build|host|card] [m steps threads print]" >&2; exit 2 ;;
esac
