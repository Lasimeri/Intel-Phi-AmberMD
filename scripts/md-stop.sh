#!/bin/bash
# md-stop.sh [card]: ask the MD run on a card to stop after its current step, through the shared
# host-memory region (works while `phi run` holds the card's control socket). Writes "STOP" at
# offset 128 of the MD stats page (src/card/mdstat.h); the run sees it at the end of the step.
set -euo pipefail
card=${1:-0}
f=/dev/shm/phi-hostmem; [ "$card" = 0 ] || f=/dev/shm/phi-hostmem-$card
[ -f "$f" ] || { echo "md-stop: $f not found" >&2; exit 1; }
off=$((1320 * 1048576 + 128))
magic=$(od -A n -t x4 -j $((1320 * 1048576)) -N 4 "$f" | tr -d ' ')
[ "$magic" = 53444d41 ] || { echo "md-stop: no MD stats block in $f" >&2; exit 1; }
printf 'STOP' | dd of="$f" bs=1 seek=$off conv=notrunc status=none
echo "md-stop: stop requested on card $card"
