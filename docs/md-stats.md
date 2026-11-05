# The MD stats block: how phitop sees an MD run on the card

An MD program on the card publishes its progress in the shared host-memory region
(`/dev/phihost` on the card, `/dev/shm/phi-hostmem` on the host; `phi-hostmem-N` for card N),
and phitop, from the 7120A stack, draws it as two `md` lines beside the core grid and the
process table. Writer: `src/card/mdstat.c` (`mdstat.h`), callable from Fortran through `bind(C)`.

## Place
Part offset **1320 MiB** of the region, first page. Claimed because [1292, 1344) MiB is unused by
the other parts on this card: the AVX-512 wrapper [0, 768), the APL engine [768, 1024), pwgd from
1024 to 1292 (its 4 MiB of rings plus 2,048 answer slots of 135,168 B), Cortex [1344, 1472), the
vision and Blackstar parts from 1472, r4d [1984, 2048). The card's swap starts at 2 GiB.

## Layout, 128 bytes, little-endian
| Offset | Type | Field |
|---|---|---|
| 0 | u32 | magic `AMDS` (0x53444D41) |
| 4 | u32 | version, 1 |
| 8 | u64 | seq: odd while the writer changes the block, even when whole |
| 16 | u64 | steps done |
| 24 | u64 | steps requested |
| 32 | u64 | time step in 1/1000 fs |
| 40 | u64 | atoms |
| 48 | u64 | team threads |
| 56 | i64 | total energy, 1/1000 kcal/mol |
| 64 | u64 | temperature, 1/1000 K |
| 72 | u64 | sum over team threads of ns inside parallel work |
| 80 | u64 | wall ns inside MD steps |
| 88 | u32 | state: 0 setting up, 1 running, 2 done, 3 failed, 4 stopped by the host |
| 92 | u32 | pid on the card |
| 96 | char[32] | run name, NUL padded |
| 128 | u32 | **stop word, written by the host**: `STOP` (0x504F5453) asks the run to end after its current step; cleared when a run opens the block |

The card maps the page uncached (`O_SYNC`), so its stores reach host memory in program order:
seq odd, the fields, seq even. phitop copies the block between two equal even reads of seq,
up to three tries.

## What phitop shows
```
md    ljmd argon  step 42 of 300   0.0 steps/s  0.03 ns/day  13070.075 ms/step  left 0:56:12  running
      E -62123.823 kcal/mol  T 72.69 K   48668 atoms, 120 threads, team busy 58.6%
```
- steps/s and ns/day: steps over phitop's sample interval (wall clock; ns/day = steps x dt x 8640 /
  interval in us, with dt in 1/1000 fs).
- ms/step: the program's own time per step over the interval.
- left: remaining steps x the program's time per step. `--` while no step finished in the interval.
- team busy: the threads' time inside parallel work over (time inside steps x threads).
- Rates need two samples from the same pid; the energy, temperature and counts are shown as written.

## Stopping a run
`scripts/md-stop.sh [card]` writes the stop word. This works while `phi run` holds the card's
control socket, which blocks every `phi run`, `phi put` and `phictl exec` until the job ends.
