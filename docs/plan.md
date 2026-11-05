# Plan: FOSS AmberTools molecular dynamics on one Xeon Phi 7120A

State on 2026-11-05 13:10: G0 passed on the host (results/2026-11-05-g0-host.md). One C probe has
run on the card (G3 team fork-join, results/2026-11-05-g3-team.md). sander has not been built for
the card: G1 (Flang) is building. Figures below are measurements where a results file is cited,
otherwise cited measurements by others or arithmetic, and are labelled. Source tree facts were
read from the checksummed tarball.

Order of targets: one 7120A first (61 cores, 244 threads, 1.238 GHz, 16 GB GDDR5, PCIe Gen2
x8). A rack of 3120A cards comes after the 7120A gates pass (stage G8) and is not started here.

## 1. Which code can be the base

| Engine | Licence | Language, size | Constant pH | KNC code present | Verdict |
|---|---|---|---|---|---|
| pmemd (Amber26) | Own licence: $0 non-commercial, no distribution of copies outside the licensee's organization, patches may be distributed | Fortran | yes | Amber 14 to 18 era `-mic_native`, `-mic_offload` (Intel compilers, Intel MPI, MPSS). Amber26 manual: 0 lines name `mic_native` | **Excluded.** Not FOSS, source cannot enter this repository |
| sander (AmberTools26) | GPL-3 by the LICENSE file (see 1.1) | Fortran + C, 242,211 lines in the directory | yes: `icnstph=1` GB, `icnstph=2` explicit water with GB titration | none (0 files match `__MIC__`, `MIC_offload`, `mic_native`, `k1om`) | **Base** |
| mdgx (AmberTools26) | GPL-3, libmdgx LGPL-3 | C, 87,653 lines; links FFTW, NetCDF and Fortran netlib | none (0 hits for protonation, titration) | none | Rejected as the base; available as a C cross-check |
| gem.pmemd (AmberTools26) | GPL-3 by the same file | Fortran, 95,431 lines | not checked; aimed at AMOEBA and GEM force fields | none | Rejected for ff14SB-class work |
| msander (AmberClassic) | GPL-2+ | Fortran | removed upstream, per its README | none | Rejected |
| mAMBER (paper: sander offload, 4.17x) | GPL-3 repository | Repository tree has 13,159 paths, 0 under `AmberTools/src/sander/`, no `short_ene`, no `offload` hit in code search | - | not in the tree | Citation only, no code to reuse |
| GROMACS 2018.8 | LGPL-2.1 | C++ | external package ST-CpHMD-v4.1.2_GMX2018 from the authors' lab | `src/gromacs/simd/impl_x86_mic` (9 files), present v2016.6 to v2021.7, absent from v2022.6 | Not Amber. See decision D1 |

### 1.1 Discrepancies found

- ambermd.org/AmberTools.php says "The sander program has the LGPL license". `AmberTools/LICENSE`
  names only libsander, libpbsa and libmdgx as LGPL-3, and everything not named as GPL-3. This
  repository follows the stricter reading, GPL-3.0-or-later. Both allow public redistribution.
- The ambermd.org download form is not needed. `ambertools26_rc7.tar.bz2`, the URL pinned by the
  conda-forge recipe, has md5 `370620f5ddb646e552a182ce8625b78e`, equal to the md5 published for
  the AmberTools26 release (measured `md5sum`, 625,923,237 B). `scripts/fetch-sources.sh` pins its
  sha256 `5d46eef3...8c9c`.
- GROMACS is a separate licence condition: "derived work must not be called official GROMACS".

## 2. The cited science does not use Amber

Thesis (Vila-Vicosa 2010, hdl 10451/5451, methods read in full) and the J. Phys. Chem. B 2012
paper. The paper is closed access (OpenAlex `is_oa: false`), so only its abstract was read; its
parameters below come from the thesis.

| Parameter | Thesis | Amber form with sander | Comparable? |
|---|---|---|---|
| Engine | modified GROMACS 3.2.1 | sander | no |
| Force field | GROMOS 53a6 (43a1 for comparison), united atom | ff14SB-class, all atom. Reference energies parametrized on ff99SB (manual 28.3.1) | no |
| Water | SPC, 16,089 molecules, rhombic dodecahedron | TIP3P class, truncated octahedron | partly |
| Long range | GRF, twin range 8/14 A, eps 54, 0.15 M | PME, cut 8 A | no |
| Titration | PB (MEAD 2.2.0) + MC (PETIT 1.5, 1e5 steps) | GB (`igb=2`) + MC | no |
| Cycle | 2 ps with 0.2 ps solvent relaxation | `ntcnstph`, `ntrelax` (test case uses 10 and 10) | yes |
| Titrated at pH 2 | Glu, Asp, C terminus | AS4, GL4. No C-terminal titratable residue in `titratable_residues.py` (AS4 GL4 TYR HIP LYS CYS DAP DCP DG DT AP CP G U HEH PRN TYX) | partly |
| Volume | 47 simulations, 1.25 us | - | - |

Result to keep in view: reversal was seen only from the early part of the misfolding transition.
No full reversal from a misfolded conformation was observed in sub-microsecond sampling. An Amber
run is a new experiment with the same question, not a reproduction of those numbers. Simulating a
conformational return says nothing about clearing infectious prions.

## 3. Hardware and toolchain state

Measured or read on 2026-11-05:

| Item | State |
|---|---|
| Compiler | KNC clang 23.1.1, musl target, patches 0001 to 0017 (Intel-Phi-LLVM e673241). `LLVM_ENABLE_PROJECTS=clang;lld`, so no Flang binary |
| Sources | `llvm-project` tree holds flang, flang-rt, mlir, openmp |
| KNC ISA in C | `kncintrin.h` plus `__builtin_ia32_knc_*` from `BuiltinsX86KNC.td` |
| Scalar FP | x87 only, no SSE, no CMOV; scalar FP returned in st0 (patch 0002) |
| Host packages missing | gcc-fortran, openmpi, tcsh |
| Host packages present | cmake 4.4.5, ninja 1.13.2, flex, bison, patch, bc, perl, fftw 3.3.11 |
| Disk | `/` 189.0 GB free, 1TB NVMe 26.6 GB free (98 %). Flang build tree goes on `/` |

Card limits that bind this project: 12 GB maximum per program, core 60 left to the OS, huge pages
touched first by the team, no polling over SSH, one card shared with other projects.

## 4. Where the resistance is

1. **No Fortran for KNC exists in this stack.** sander is Fortran 90. Everything else waits on this.
2. **x87 assumptions in Flang and flang-rt.** `flang/lib/Optimizer/CodeGen/Target.cpp` classifies
   REAL arguments as SSE for x86-64, `flang-rt/lib/runtime/exceptions.cpp` and `edit-input` reference
   `xmmintrin`, `_mm_` or `mxcsr`. COMPLEX ABI against C `_Complex` under patch 0002 is unverified.
3. **`libomp` assembly uses SSE registers.** `openmp/runtime/src/z_Linux_asm.S` has 8 `xmm` or
   `mxcsr` lines (`__kmp_xchg_real32/64`).
4. **sander's threading is the wrong shape.** The only `!$omp` loop in `short_ene.F90` (line 405)
   sits inside the per-atom routine: one fork-join per atom, about 50,000 per step at this size.
   On a 120-thread team this cannot work. The whole direct-space loop has to be restructured.
5. **Amdahl dominates, not the pair kernel.** See section 5. One serial O(N) loop consumes the
   entire step budget.
6. **Double-precision gaps on KNC.** No F64 gather, no F64 rcp or rsqrt, so F32 seed plus Newton,
   with a seed range limit of about 1.18e-38 to 3.4e38.
7. **Excess precision.** Scalar x87 results differ from SSE host results, so card and host cannot
   be compared bit for bit without control over x87 precision.
8. **Input tools are Python or Perl.** `cpinutil.py`, `cpinutil.pl`, parmed. Ready cpin files exist
   in `test/cnstph`. New systems need a C generator, `mkcpin`.

## 5. Budget (arithmetic, not measurement)

System: thesis water count, 16,089 x 3 = 48,267 atoms, plus an estimated 2,250 protein atoms for
residues 90-231: N about 50,500. Density about 0.1 atoms per cubic A. Step 2 fs.

| Quantity | Value | Formula |
|---|---|---|
| Steps per second for 1 ns/day | 5.787 | 5e5 / 86,400 |
| Step budget at 20 ns/day | 8.64 ms | 1 / (20 x 5.787) |
| Pairs inside 8 A | 5.41 M | 50,500 x (4/3)(pi)(512)(0.1) / 2 |
| Direct space, 60 cores | 0.29, 0.87, 2.19 ms at 4, 12, 30 cycles per pair | pairs x c / (60 x 1.238e9) |
| One thread on one core | 619 M issue slots/s | KNC does not issue from one thread on consecutive cycles |
| **One serial loop at 100 cycles per atom** | **8.16 ms** | 50,500 x 100 / 619e6 |

Published anchors:

| Source | Number | Relevance |
|---|---|---|
| Amber list, 2014-09-02 | 0.05 ns/day, serial `pmemd.mic_native`, 1 of 61 cores, atom count not given | Scalar floor |
| Amber list, 2017-03-14 | 2 ns in 38 h = 1.263 ns/day, 262,572 atoms, 7120P with `pmemd.mic_offload.MPI`, 20 host ranks | Scaled to 50,500 atoms: 6.57 ns/day (linear estimate, host cores included) |
| Needham 2016 | about 35 % gain, 1,067,095 atoms, 2 x E5-2697 v2 + 7120P; best above 400,000 atoms | Offload does not fit 50,000 atoms: native mode only |
| Pennycook 2013 | KNC 5.2x force and 3.8x neighbour list from 512-bit SIMD over scalar; gathers dominate | Kernel shape |
| Ma 2015 | Inter-intra hybrid vectorization 1.25x force, 1.40x list over auto-vectorized code on 7120P | Kernel shape |
| GROMACS forum 2020 | Forum reply: modern CPUs are faster than Xeon Phi | The host baseline G0 will quantify this |

Targets: beat the host, 4.659 ns/day measured at G0 (best of 8 ranks on the 5800X), and the 6.57 ns/day
estimate, then 20 ns/day. Volume check: 1.25 us takes 190
days at 6.57 ns/day, 62.5 days at 20 ns/day, on one card. The 20 ns/day target is an estimate
that no measurement backs yet.

## 5a. Card support: a patch series, never edits in place

AmberTools has no notion of this card's software stack: no `knc` CPU, no `knc-cc` and
`knc.cfg`, a musl sysroot instead of glibc, x87-only scalar FP, no MPSS, and programs reach the
card through `phi put` and `phi run`. Support goes in as numbered patches against the pinned
tarball (`patches/ambertools/SERIES`, applied by a script that refuses a tree whose checksum
differs), plus files that live only in this repository. Survey of the source on 2026-11-05:

| What | Found in AmberTools26 | Patch or file |
|---|---|---|
| Cross build | `CROSSCOMPILE` option already exists for `GNU`, `CLANG`, `MANUAL`, needs `TARGET_TRIPLE` and `CMAKE_FIND_ROOT_PATH` (`cmake/AmberCompilerConfig.cmake`) | Repository file `cmake/knc-toolchain.cmake`: `knc-cc --config knc.cfg`, Flang with `-march=knc`, triple `x86_64-unknown-linux-musl`, root = the stack sysroot |
| Compiler identity | `ALL_COMPILER_VALUES` has no Flang; Fortran flags keyed on `GNU`, `Intel`, `PGI` | Patch: a `KNC` compiler value, LLVMFlang flag branch |
| Host tuning | `-mtune=native` added for GNU C, C++ and Fortran (`cmake/CompilerFlags.cmake` lines 76, 118, 154) | Patch: never for KNC |
| Configure-time runs | `check_fortran_source_runs` uses `try_run`; `TargetArch.cmake` try_run only compiles | Repository file `scripts/phi-emu.sh` as `CMAKE_CROSSCOMPILING_EMULATOR`: puts the probe on the card and runs it, so answers come from the card, not guesses |
| glibc, SSE, inline asm | 0 files in sander, sqm, pbsa, sff, dlfind, lib match `__GLIBC__`, `execinfo`, `mallinfo`, `fpu_control`, `feenableexcept`, `asm`, `cpuid`, `_mm_*`, `xmmintrin`, `__SSE` | none needed so far |
| NetCDF 4.6.1 and netcdf-fortran 4.4.4 | built in-tree (`INTERNAL`) by autotools | Patch or prebuilt: configure with `--host=x86_64-unknown-linux-musl` against the stack sysroot |
| Shell-outs | `system('rm ...')`, `system('mv ...')` only in QM/MM paths | not on the MD path; QM/MM stays disabled |
| Tests | `Run.*` are `/bin/sh`, `dacdif` uses awk | sander runs on the card, `phi get` brings mdout back, dacdif runs on the host |

Card rules that the scripts enforce, not just document: total reservation at most 12 GB, `/tmp` on
the card is tmpfs and counts toward it; 120 threads on cores 0 to 59, core 60 left alone; one
`phi run` at a time (it holds the control socket); no SSH polling of the card; card timestamps are
never used (its clock is not the host's).

## 6. Stages and gates

Each stage ends with a measured result in `docs/results/`, committed and pushed. A failed gate is
recorded as failed. Card stages need the card to be free.

### G0: host reference (no card): PASSED 2026-11-05, see results/2026-11-05-g0-host.md
Build serial `sander` with gfortran on the host from the pinned tarball, minimal CMake (`-DMPI=FALSE
-DCUDA=FALSE -DDOWNLOAD_MINICONDA=FALSE`, no Python). 14 jobs at most.
- Gate: `test/cnstph/implicit` and `test/cnstph/explicit` pass `dacdif` against the shipped
  `mdout.save`.
- Measure: ns/day on the Ryzen 5800X for the explicit case and a 50,500-atom PrP box. This is the
  number the card must be compared with, interleaved and repeated.

### G1: Fortran for KNC
Add `flang;mlir` to a build of the patched LLVM 23.1.1 on `/`. RAM use of this build is unmeasured:
start at 8 jobs, log MemAvailable, raise toward 14 only if it stays clear, since host RAM plus swap
above 90 % kills the session. Build flang-rt for `x86_64-unknown-linux-musl` with `-march=knc`. Patch
what the audit finds. Add `LLVMFlang` to AmberTools' compiler handling (it has no entry today).
- Gate on the card: formatted I/O round trip, NAMELIST read (sander's input), module allocatables,
  COMPLEX arguments and returns against KNC clang `_Complex`, Fortran-to-C struct call.
  `phi-isa-audit` reports 0 illegal instructions.
- Fallback if Flang resists: gfortran from upstream GCC with an x87 return ABI patch. Not `f2c`,
  because sander is Fortran 90.

### G2: unmodified sander on the card
Cross-build serial sander and its in-tree dependencies (`libsqm`, `libpbsa`, `sff_fortran`,
`dlfind`, NetCDF 4.6.1, netlib) with Flang. Disable RISM, MPI, CUDA, Python, plumed, torchani.
- Gate: same two cnstph tests on the card. Differences against the host recorded per field. Test
  x87 precision control (`fldcw` to 53-bit at entry) as the fix for excess precision: keep it only
  if it brings the energy differences down and the speed holds.
- Measure: scalar ns/day, plus a per-routine time profile. The profile decides G3 to G6 order.

### G3: thread team
Single process, 120 threads on cores 0 to 59, no MPI on one card. Two candidates, chosen by measured
fork-join latency on the card: `libomp` with the SSE lines replaced, or a spin team with `delay`,
callable from Fortran through `iso_c_binding`.
- Gate: empty parallel region cost in cycles, at least 1,000 samples, both candidates, reported.
- Spin team measured 2026-11-05 (results/2026-11-05-g3-team.md): 8-ary gather tree p50 10,340 cycles
  (8.37 us), plain pthread barrier 1,744,350. libomp not measured yet.

### G4: direct-space kernel in MVEX
C with KNC builtins, called from `short_ene`'s caller, replacing the per-atom loop. Generated by a
small C generator in the manner of `knc-c-gen`, so all variants come from one source.
- F64 first, because that allows comparison with the reference. SPDP second (F32 pair, F64
  accumulators), with Amber's published tolerance in view: `ambermd.org/intel/midpoint.htm`,
  max deviation 8.6e-4 (DHFR) to 4.0e-3 (STMV) for SPDP.
- Techniques, each kept only if measured at least as fast as the plain form on the 7120A and
  matching forces: water-water tile with 9 pairs fused and no gather; branchless cutoff through
  compare masks; exclusions as a bit mask tested in `k` registers; hand-placed SoA arrays on huge
  pages, touched by the team first, with a prime stride to avoid L2 aliasing; two-stage prefetch
  (`vprefetch1`, then `vprefetch0`); per-thread force buffers merged once per step.
- Gate: per-atom force difference against scalar sander below a stated tolerance on a stored frame;
  NVE drift over 1 ns compared to the host; pair cost in cycles per pair.

### G5: everything else that is O(N)
Pair list, PME reciprocal (B-spline spread, 3D FFT, gather), bonded terms, SHAKE or SETTLE,
integrator, thermostat, I/O path. sander's default FFT is `pubfft` (scalar Fortran, F64). FFTW's
`--enable-kcvi` is single precision only (FFTW 3.3.11 docs; the bundled `fftw-3.3` has no kcvi), so
an F64 radix kernel in MVEX is written here, gated against `pubfft` output.
- Gate: profile on the card shows no single-thread routine above 10 % of a step. ns/day reported
  with spread over interleaved repeats, against G0.

### G6: constant pH at speed
Thread the GB titration step (`egb.F90` 3,238 lines, `constantph.F90` 1,437 lines) and the relaxation
path. Write `mkcpin` in C from the reference data in `titratable_residues.py`.
- Gate: `test/cnstph` on the threaded build. Lysozyme (HEWL) pKa series from Amber tutorial 18 at
  several pH values compared to the tutorial and to experiment.

### G7: PrP pilot on the 7120A
Build the model with tleap: PDB 1QM0 holds residues 125-228 only, so residues 90-124 come from 2PRP
as the thesis did (structures 6 and 9). Solvate, equilibrate at pH 7, run pH 2 replicates, then
pH 7 jumps from early and late frames, the protocol of section 2.
- Gate: stable pH 7 control in RMSD and DSSP. Physical conclusions only with replicate counts stated.

### G8: rack (out of scope until G5 passes)
One independent trajectory per card; no domain decomposition across PCIe. pH-REMD through a host
coordinator instead of MPI between cards. Rebuild with 3120A constants (57 cores, 6 GB). The number
of cards and GPUs in the rack has to be confirmed then: the record on hand and the stated plan differ.

## 7. Decisions

| Id | Question | Default taken here |
|---|---|---|
| D1 | Amber stays the target, or GROMACS 2018.8 with ST-CpHMD, which is the method and code family of the cited papers, has KNC SIMD kernels, and needs no Fortran? | Amber |
| D2 | Vendor the sander subset into this repository, or fetch by checksum? | Fetch for G0 and G1. At G2 the exact file list is known: vendor it under `upstream/`, GPL-3 permits it, so the build survives the URL disappearing |
| D3 | Project licence | GPL-3.0-or-later: changes to GPL-3 sander cannot be distributed under anything else |

## 8. Design

```
╔═ Desktop host (Ryzen 5800X, 31 GB) ══════════════════════════════════════════════════╗
║ ┌──────────────────────────┐   ┌──────────────────────────┐   ┌──────────────────┐   ║
║ │ AmberTools26 tarball     │──▶│ KNC LLVM 23.1.1 + Flang  │──▶│ sander.knc       │   ║
║ │ sha256 5d46eef3...       │   │ flang-rt, musl (G1)      │   │ phi put          │   ║
║ └────────────┬─────────────┘   └──────────────────────────┘   └────────┬─────────┘   ║
║              ▼                                                         │             ║
║ ┌──────────────────────────┐   ┌──────────────────────────┐            │             ║
║ │ gfortran sander (G0)     │──▶│ mdout, forces, ns/day    │            │             ║
║ │ reference build          │   │ compared in results/     │            │             ║
║ └──────────────────────────┘   └────────────▲─────────────┘            │             ║
╚═════════════════════════════════════════════│══════════════════════════│═════════════╝
                                              │ mdout, profile           │ PCIe Gen2 x8 
╔═ Xeon Phi 7120A (2f:00.0, 12 GB cap) ═══════│══════════════════════════▼═════════════╗
║ ┌───────────────────────────────────────────┴────────────────────────────────────┐   ║
║ │ sander Fortran shell: input, cpH MC, integrator control, output (G2)           │   ║
║ └───────┬─────────────────────┬──────────────────────┬───────────────────────────┘   ║
║         ▼                     ▼                      ▼                               ║
║ ┌───────────────┐   ┌───────────────────┐   ┌────────────────────┐                   ║
║ │ direct space  │   │ PME: spline, FFT  │   │ GB titration, SHAKE│                   ║
║ │ MVEX (G4)     │   │ MVEX F64 (G5)     │   │ bonded (G5, G6)    │                   ║
║ └───────┬───────┘   └─────────┬─────────┘   └──────────┬─────────┘                   ║
║         ▼                     ▼                        ▼                             ║
║ ┌────────────────────────────────────────────────────────────────────────────────┐   ║
║ │ thread team, 60 cores x 2 threads (G3), SoA arrays on huge pages               │   ║
║ └────────────────────────────────────────────────────────────────────────────────┘   ║
╚══════════════════════════════════════════════════════════════════════════════════════╝
```

## 9. Sources

- Amber KNC support: https://ambermd.org/intel/ ; Amber 14 manual section: https://ambermd.org/MICv2_Amber14.pdf
- Amber licences: https://ambermd.org/GetAmber.php , https://ambermd.org/AmberTools.php
- Amber 2026 Reference Manual: https://ambermd.org/doc12/Amber26.pdf (chapter 28, constant pH)
- Tarball URL and sha256: https://github.com/conda-forge/ambertools-feedstock/blob/master/recipe/meta.yaml
- AmberClassic: https://github.com/Amber-MD/AmberClassic ; mAMBER: https://github.com/zhangxiaoyu11/mAMBER , https://doi.org/10.1109/ccgrid.2015.66
- Needham, Bhuiyan, Walker, Comput. Phys. Commun. 201, 95 (2016): https://doi.org/10.1016/j.cpc.2015.12.025
- Serial native figure: http://archive.ambermd.org/201409/0027.html ; offload figure: http://archive.ambermd.org/201703/0171.html
- Mixed precision: https://ambermd.org/intel/midpoint.htm
- Thesis: http://hdl.handle.net/10451/5451 ; paper: https://doi.org/10.1021/jp3034837 ; earlier work: https://pubs.acs.org/jpcbfk/article/114/39/12692/990525/
- Constant pH: Swails et al. 2014, https://theory.rutgers.edu/resources/pdfs/Swails_JChemTheoryComput_2014_v10_p1341.pdf ; tutorial: https://ambermd.org/tutorials/advanced/tutorial18/index.php
- Stochastic titration package: https://itqb.unl.pt/labs/molecular-simulation/in-house-software
- GROMACS on KNC: https://manual.gromacs.org/documentation/2019/install-guide/index.html ; forum: https://gromacs.bioexcel.eu/t/intel-xeon-phi/762
- SIMD for MD on KNC: Pennycook 2013, https://aiichironakano.github.io/cs653/Pennycook-PhiMD-IPDPS13.pdf ; Ma 2015, https://aiichironakano.github.io/cs653/Ma-XeonPhiOpt-HPCC15.pdf
- Flang runtime: https://llvm.googlesource.com/llvm-project.git/+/refs/heads/main/flang-rt/README.md
- FFTW options: https://www.fftw.org/fftw3_doc/Installation-on-Unix.html
- KNC toolchain: https://github.com/Lasimeri/Intel-Phi-LLVM
