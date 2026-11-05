#!/bin/bash
# G3 probe: build src/probe/g3-team.c for the card with the stack's knc-cc, audit it,
# check the logic on the host, then run it on the card. Output: logs/g3/team-<host|card>.tsv
# Columns: form round threads min p50 p99 max (cycles) ns_per_tick.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
: "${PHI_TREE:=/media/1TB-NVMe/Xeon-Phi/Intel Phi 7120A}"
: "${PHI_CARD:=0}"
THREADS=${THREADS:-120} ITERS=${ITERS:-10000} ROUNDS=${ROUNDS:-3}
out=$R/out/probe; log=$R/logs/g3
mkdir -p "$out" "$log"
flags=(-std=gnu11 -O2 -Wall -Wextra -Werror -pthread)

cc "${flags[@]}" -o "$out/g3-team.host" "$R/src/probe/g3-team.c"
( . "$PHI_TREE/toolchain/env.sh"; knc-cc "${flags[@]}" -static -o "$out/g3-team.knc" "$R/src/probe/g3-team.c" )
line=$("$PHI_TREE/host/target/debug/phi-isa-audit" "$out/g3-team.knc" | /usr/bin/grep -a '^result:')
echo "audit: $line"
[[ $line == "result: 0 illegal, 0 suspect,"* ]] || { echo "g3: audit failed" >&2; exit 1; }

hdr=$'form\tround\tthreads\tmin\tp50\tp99\tmax\tns_per_tick'
{ echo "$hdr"; "$out/g3-team.host" 16 2000 1; } > "$log/team-host.tsv"
echo "host logic check: ok"

phi -c "$PHI_CARD" put "$out/g3-team.knc" /tmp/g3-team < /dev/null
{ echo "$hdr"; phi -c "$PHI_CARD" run /tmp/g3-team "$THREADS" "$ITERS" "$ROUNDS" < /dev/null; } > "$log/team-card.tsv"
phi -c "$PHI_CARD" run rm -f /tmp/g3-team < /dev/null
cat "$log/team-card.tsv"
