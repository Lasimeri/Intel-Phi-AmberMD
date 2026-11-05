# phitop md lines and the 3090 Ti view, 2026-11-05

## Test simulation
`tests/md/ljmd.f90`: liquid argon, Lennard-Jones, 48,668 atoms, box 132.126 A, cut 8.5 A, 4 fs steps,
velocity Verlet, Berendsen. Built with this repository's Flang for the card, linked with
`src/card/team.c` (the G3 8-ary tree team), `mdstat.c`, `mdframe.c`. Audit 0 illegal, 0 suspect.
It is a load generator with real physics, not Amber.

## A bug found by phitop
First card run: 2 steps a minute, `ljmd` at 198 % CPU, the grid 0.8 % busy. musl answers
`sysconf(_SC_NPROCESSORS_ONLN)` from the calling thread's affinity mask; the team pinned the master
before creating the workers, so 119 workers saw 1 CPU and all pinned themselves to CPU 0. Fixed by
reading the CPU count once before pinning. The run held the control socket and could not be killed
through `phi` or SSH; the card was rebooted (30.7 s). Since then: the stop word (`scripts/md-stop.sh`).

## Monitoring, 600 steps, 120 threads (raw: 2026-11-05-ljmd-card-600.log, -phitop-md-lines.txt)
| | ljmd's own log | phitop |
|---|---|---|
| speed | 3.047 to 3.217 steps/s (since start) | 2.9 to 3.4 steps/s, 311 to 321 ms/step |
| ns/day | 3.2 x 4 fs x 86,400 = 1.106 | 1.02 to 1.19 |
| team busy | 80.5 to 84.8 % | 83.3 to 85.0 % |
| time left at step 95 | 505 x 0.315 s = 159 s | 0:02:39 |
Process line: `ljmd` 11,943.7 % (119.4 hardware threads), 120 threads. Grid: threads 0 and 1 of
cores 0-59 full, threads 2 and 3 empty, core 60 untouched. The grid counts spinning as busy; team
busy is the share of time in parallel work.

## View, 60 s with frames every step (raw: 2026-11-05-amview-60s.txt, image -amview-snap.png)
RTX 3090 Ti, OpenGL 4.6.0 NVIDIA 615.78.08, window 1920 x 1080.
| | |
|---|---|
| frames drawn | 60.00 every second (59.99 in the first) |
| loop interval | p50 16.63 to 16.66 ms, p99 16.67 to 18.31 ms (23.07 ms in the first second) |
| GPU time per frame | 0.52 to 3.62 ms (rose from 0.5 to 3.4 in the first 8 s; cause not checked) |
| MD frames received | 2.95 to 4.00 a second: every step of the run |
| pixels covered (`--snap`) | 37.6 % |
Card side: 700 frames, 1.643 ms per frame to fill and publish (team writes straight into the
uncached region). The run with frames reached 3.311 steps/s, the run without 3.217: frames do not
slow it measurably.
