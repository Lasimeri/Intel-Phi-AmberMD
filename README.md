# Intel-Phi-AmberMD

Molecular dynamics from the FOSS parts of the Amber suite (AmberTools, chiefly `sander` with its
constant-pH methods) on Intel Xeon Phi Knights Corner coprocessors, with no Intel compilers,
no MPSS and no MKL. Compiler: the patched LLVM from
[Intel-Phi-LLVM](https://github.com/Lasimeri/Intel-Phi-LLVM).

First target: one Xeon Phi 7120A. Later target: a rack of 3120A cards.

**Status: plan only. Nothing has been built for the card and nothing has run on a card.**
Read [docs/plan.md](docs/plan.md).

`pmemd` is not used and its source must never be added here: the Amber licence forbids
distributing copies. Only AmberTools code, which is GPL or LGPL, is a base for this work.

`scripts/fetch-sources.sh` fetches and verifies the AmberTools26 source.

Licence: GPL-3.0-or-later (see LICENSE). AmberTools components keep their own licences, listed in
their `AmberTools/LICENSE`.
