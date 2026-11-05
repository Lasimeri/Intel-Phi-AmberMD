# Frames: the simulation drawn by the host's GPU

The card publishes atom positions into the shared host-memory region; `amview` on the host reads
the newest whole frame and draws it on the RTX 3090 Ti with OpenGL 4.6 at 1920 x 1080, 60 frames a
second. Writer: `src/card/mdframe.c` (`mdframe.h`). Viewer: `src/host/amview.c`.

## Place and layout
Part offset 1320 MiB + 64 KiB, after the stats page (`docs/md-stats.md`), limited to 52 MiB so it
ends below 1344 MiB. A 4 KiB header page, then 3 slots of natoms x 16 bytes, each rounded up to
4 KiB. An atom is 4 floats: x, y, z in Angstrom and w, a value the viewer colours by (ljmd: the
speed in A/ps). 48,668 atoms: 778,688 B a slot.

| Offset | Type | Field |
|---|---|---|
| 0 | u32 | magic `AMDF` (0x46444D41) |
| 4 | u32 | version, 1 |
| 8 | u64 | latest complete frame number, 0 = none |
| 16 | u64 | atoms |
| 24 | u32 | slots, 3 |
| 28 | u32 | bytes per atom, 16 |
| 32 | u64 | bytes per slot |
| 40 | u64 | offset of slot 0, 4096 |
| 64 + 64 s | 64 B | slot s: u64 seq (frame number, or all ones while written), u64 step, float box[3], float wmax |

## Protocol
Frame k goes to slot k % 3: the writer sets the slot's seq to all ones, the team writes its atoms
(each thread its own range, straight into the region), then the step, box and wmax, then seq = k,
then latest = k. The reader takes latest, copies slot latest % 3 to host memory, and uploads it
only if the slot's seq was k before and after the copy and latest has not moved two frames on.

## amview
- Every atom is one instanced quad; the fragment shader intersects the view ray with the sphere
  and writes the exact depth (`layout(depth_greater)`), so spheres intersect correctly.
- Colour: w / wmax on a blue-cyan-yellow-red ramp. The periodic box is drawn as lines.
- Pacing: swap interval 0 and `clock_nanosleep` to absolute 1/60 s deadlines.
- Mouse drag rotates, the wheel zooms, space stops the rotation, q or Esc quits.
- Each second it prints frames drawn, loop interval p50 and p99, GPU time per frame from GL timer
  queries, MD frames received and the run's step; `--snap FILE` saves its own frame once as PPM.
