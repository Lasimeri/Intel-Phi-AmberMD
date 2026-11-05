#!/bin/bash
# G1 gate: Fortran conformance on the card. Builds tests/fortran/conf.f90 + conf_c.c twice:
#   card: this repository's Flang + its clang, -march=knc, static musl, flang-rt for KNC
#   host: gfortran + gcc (the G0 compilers)
# Audits the card binary, runs both, and compares every HEX line bit for bit.
# Output: logs/g1/conf-card.txt, conf-host.txt, conf-hexdiff.txt
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
O=${G1_OUT:-$HOME/.cache/intel-phi-ambermd/llvm-flang}
BIN=$O/build/bin
SYSROOT=${PHI_SYSROOT_7120A:-$HOME/.cache/intel-phi-7120a-build/toolchain/sysroot}
: "${PHI_TREE:=/media/1TB-NVMe/Xeon-Phi/Intel Phi 7120A}"
: "${PHI_CARD:=0}"
T=$R/tests/fortran; out=$R/out/g1; log=$R/logs/g1
mkdir -p "$out" "$log"
knc=(--target=x86_64-unknown-linux-musl -march=knc --sysroot="$SYSROOT" -O2)

(cd "$out" && "$BIN/clang" "${knc[@]}" -fno-jump-tables -c "$T/conf_c.c" -o conf_c.knc.o \
  && "$BIN/flang" "${knc[@]}" -static -fuse-ld=lld -rtlib=compiler-rt -module-dir "$out" \
     "$T/conf.f90" conf_c.knc.o -o conf.knc -lc++ -lc++abi -lunwind)
line=$("$PHI_TREE/host/target/debug/phi-isa-audit" "$out/conf.knc" | /usr/bin/grep -a '^result:')
echo "audit: $line"
[[ $line == "result: 0 illegal, 0 suspect,"* ]] || { echo "g1: audit failed" >&2; exit 1; }

(cd "$out" && gcc -O2 -c "$T/conf_c.c" -o conf_c.host.o && gfortran -O2 -J "$out" "$T/conf.f90" conf_c.host.o -o conf.host)
(cd "$out" && ./conf.host) > "$log/conf-host.txt" || true

phi -c "$PHI_CARD" put "$out/conf.knc" /tmp/g1conf/conf < /dev/null >/dev/null 2>&1 || {
  phi -c "$PHI_CARD" run mkdir -p /tmp/g1conf < /dev/null; phi -c "$PHI_CARD" put "$out/conf.knc" /tmp/g1conf/conf < /dev/null; }
phi -c "$PHI_CARD" run sh -c 'cd /tmp/g1conf && ./conf; echo "exit $?"' < /dev/null > "$log/conf-card.txt" || true
phi -c "$PHI_CARD" run rm -rf /tmp/g1conf < /dev/null

echo "== card"; cat "$log/conf-card.txt"
echo "== host"; /usr/bin/grep -E '^(RESULT|FAIL)' "$log/conf-host.txt" || true
join <(awk '$1=="HEX"{print $2, $3}' "$log/conf-card.txt" | sort) \
     <(awk '$1=="HEX"{print $2, $3}' "$log/conf-host.txt" | sort) \
  | awk '{print $1, $2, $3, ($2==$3 ? "same" : "DIFF")}' > "$log/conf-hexdiff.txt"
echo "== hex card vs host"; cat "$log/conf-hexdiff.txt"
echo "bit-identical: $(awk '$4=="same"' "$log/conf-hexdiff.txt" | wc -l) of $(wc -l < "$log/conf-hexdiff.txt")"
