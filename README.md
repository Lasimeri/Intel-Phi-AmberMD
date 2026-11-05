# Intel-Phi-AmberMD

Molecular dynamics from the FOSS parts of the Amber suite (AmberTools, chiefly `sander` with its
constant-pH methods) on Intel Xeon Phi Knights Corner coprocessors, with no Intel compilers,
no MPSS and no MKL. Compiler: the patched LLVM from
[Intel-Phi-LLVM](https://github.com/Lasimeri/Intel-Phi-LLVM).

First target: one Xeon Phi 7120A. Later target: a rack of 3120A cards.

**Status (2026-11-05):** G0 passed: host sander passes the constant-pH tests and sets the bar,
4.659 ns/day on a 49,910-atom PrP box (8 ranks, Ryzen 5800X). G3 probe ran on the 7120A: a
120-thread fork-join costs 8.37 us with an 8-ary gather tree, 1.41 ms with pthread barriers.
G1 passed: Flang builds Fortran for the card; 12 of 12 conformance checks pass on the 7120A.
sander does not run on the card yet: that is G2.
Monitoring and view: phitop shows MD runs on the card (two `md` lines, docs/md-stats.md), and
`amview` draws them on the host's RTX 3090 Ti at 1080p60 (docs/md-frames.md), tested with an argon run.
Results: [docs/results](docs/results).
Read [docs/plan.md](docs/plan.md).

`pmemd` is not used and its source must never be added here: the Amber licence forbids
distributing copies. Only AmberTools code, which is GPL or LGPL, is a base for this work.

`scripts/fetch-sources.sh` fetches and verifies the AmberTools26 source.

Licence: GPL-3.0-or-later (see LICENSE). AmberTools components keep their own licences, listed in
their `AmberTools/LICENSE`.
