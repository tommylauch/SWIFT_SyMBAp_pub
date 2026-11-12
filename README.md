# SWIFT with SyMBAp

This package is the [SWIFT](https://www.boulder.swri.edu/~hal/swift.html) N-body
subroutine package by Hal Levison and Martin Duncan (retrieved Jan 2019),
extended with **SyMBAp**, an OpenMP-parallel version of the SyMBA integrator by
Tommy Chi Ho Lau and Man Hoi Lee. SWIFT integrates a set of mutually interacting massive bodies,
optionally with test particles that feel the massive bodies but do not affect
them or each other. The original drivers and their input/output are unchanged;
SyMBAp reads the same files as SyMBA.

The pebble/gas accretion and migration model built on SyMBAp lives in the
`gas_giant` submodule (see `gas_giant/doc/model.tex`).

Contents: [Building](#building) · [Drivers](#drivers) · [Running](#running) ·
[SyMBAp](#symbap-parallel-symba) · [Input files](#input-files) ·
[Output files](#output-files) · [Status flags](#test-particle-status-flags) ·
[Oblateness](#oblateness-j2-j4) · [Unit numbers](#fortran-unit-numbers) ·
[Changelog](#symbap-changelog)

## Building

Requirements: a Fortran compiler (GNU Fortran is tested; GCC 8.5 is the oldest
version known to build this branch, GCC 12.3 is used in production), GNU Make,
`cpp`, and OpenMP for SyMBAp.

```bash
make            # libswift.a and all drivers in main/ (same as: make main)
make library    # libswift.a only
make tool       # libswift.a and the tools in tools/
make USE_OPENMP=0   # everything except swift_symba5p, without OpenMP
make clean
make help       # targets, flags and the programs that will be built
```

Edit the top of `Makefile` if needed:

| Variable   | Meaning |
|------------|---------|
| `FC`       | Fortran compiler (default `gfortran`) |
| `FFLAGS`   | compiler flags without `-c` (default `-O2 -ftree-vectorize -march=native`). When compiling on a different CPU than the one that runs the code (e.g. a cluster login node), replace `-march=native` with the target, e.g. `-march=znver2`. |
| `CPP`, `CPPFLAGS` | C preprocessor and its flags, see [FXDR](#fxdr) |
| `USE_OPENMP` | `1` builds SyMBAp with `-fopenmp` (added automatically for `symba5p/` and `swift_symba5p`) |

`swift_rmvs3` and `swift_rmvs4` (and the `rmvs3_step*`/`rmvs4` routines they use)
are compiled with `-mcmodel=medium` because their static arrays exceed 2 GB in
total; nothing else uses the medium code model.

### FXDR

SWIFT is standard Fortran 77 except for one optional feature. XDR is a
platform-independent binary format SWIFT can use for `bin.dat`; it requires the
[FXDR library](http://meteora.ucsd.edu/~pierce/fxdr_home_page.html). With FXDR
set `CPPFLAGS = -D_FXDR_AVAIL`, otherwise keep the default `-U_FXDR_AVAIL`.
Requesting XDR output without FXDR makes SWIFT stop with an error.

## Drivers

All drivers in `main/` take essentially the same input and produce the same type
of output. Each step starts from heliocentric positions and velocities and
advances them by `dt`; internal substeps (e.g. in Bulirsch–Stoer) are hidden from
the user. The calls look identical, so one integrator can be swapped for another.

| Program | Method |
|---------|--------|
| `swift_mvs` | Wisdom–Holman mapping (WHM; Wisdom & Holman 1991, AJ 102, 1528). Test particles are removed at close planetary encounters; arbitrarily close solar encounters are allowed. |
| `swift_rmvs3` | Regularized Mixed Variable Symplectic method (Levison & Duncan 1994, Icarus 108, 18): WHM extended to close approaches between test particles and planets. |
| `swift_rmvs4` | RMVS variant, includes close encounters. |
| `swift_bs` | Bulirsch–Stoer (Press et al. 1992, *Numerical Recipes*). Same I/O as WHM, plus a prompt for the tolerance EPS (1e-8 is common; can be appended to the input file). Handles close encounters but is about 10× slower than RMVS. |
| `swift_tu4` | Fourth-order T+U symplectic method (Candy & Rozmus 1991, J. Comput. Phys. 92, 230; Gladman, Duncan & Candy 1991, CeMDA 52, 221). Mostly used to check the other integrators. |
| `swift_helio` | Integrator in heliocentric coordinates (`helio/`). |
| `swift_skeel` | Includes close encounters (`skeel/`). |
| `swift_lyap` | Lyapunov exponents of test particles from shadow particles read from `shadow.in` (same format as `tp.in`), renormalised every `dtout`; writes log10 of the exponent to `lyap.out`. Talk to the authors before using it. |
| `swift_lyap2` | Lyapunov exponents from the difference equations; asks how often to output and for the initial phase-space separation; writes `lyap2.out`. Still being tested. |
| `swift_symba5` | SyMBA (Duncan, Levison & Lee 1998, AJ 116, 2067): symplectic integrator that handles close encounters between massive bodies. No test particles. |
| `swift_symba5p` | SyMBAp, the OpenMP-parallel SyMBA, see [below](#symbap-parallel-symba). |

## Running

The drivers prompt for their input on standard input, so the usual way is to
redirect a small file listing the answers:

```bash
cd example
../main/swift_mvs < mvs.in          # mvs.in: param.in, pl.in, tp.in
```

| Driver | Standard input, one entry per line |
|--------|------------------------------------|
| WHM, RMVS, TU4, helio, skeel, lyap | parameter file, planet file, test-particle file |
| `swift_bs` | as above, then EPS |
| `swift_symba5` | parameter file, planet file, `mtiny` (smallest mass to self-gravitate) |
| `swift_symba5p` | parameter file, planet file, `mtiny`, lower limit of the particle-to-thread ratio, max. no. of threads |

### Examples

`example/` holds the CHO initial conditions for the four giant planets with Pluto
as a test particle (`param.in`, `pl.in`, `tp.in`, `mvs.in`; dt = 1 yr = 365.25 d,
1 Myr). `example/symba/` holds the Sun and the four giant planets (no test
particle) for SyMBA (`symba.in`) and SyMBAp (`symbap.in`: `th_low` = 20,
`th_max` = 2). The planetary system is not a good Solar
System; do not use it for science.

## SyMBAp (parallel SyMBA)

`swift_symba5p` parallelises SyMBA with OpenMP: the triangular pair loops are
flattened and vectorised, and close-encounter pairs are grouped so that each group
is integrated in parallel. It reads the same parameter and planet files as
`swift_symba5`, plus two thread parameters (see [Running](#running)):

- **Lower limit of the particle-to-thread ratio** `th_low`: at most `nbod/th_low`
  threads are used, re-checked whenever bodies are removed.
- **Max. no. of threads** `th_max`.

The run starts with `max(min(nbod/th_low, th_max), 1)` threads. With more than one
thread, `symba5p_tune` times 10 000 steps at the current count and at one thread
fewer, keeps the faster one, and repeats after every dump and every change in the
number of bodies. The chosen count is stored in `threads.dat` in the run directory
and read back at the next tuning. Writing `-n` into `threads.dat` fixes the count
at `n` and switches the timing off (still capped by `th_low`); a positive value is
used as the starting point of the next tuning.

Notes:

- Compile with `-fopenmp` (done by the Makefile). gfortran puts local arrays on the
  stack with `-fopenmp`; with the pre-set array sizes this needs a large stack:
  `ulimit -s unlimited` (or a large value) and `export OMP_STACKSIZE=<size>`
  (e.g. `12m`) for the OpenMP threads.
- `OMP_NUM_THREADS` is not needed; the thread count comes from the input.
- On SIGTERM the run stops cleanly after the current step without an extra dump;
  restart from the last dump as usual.
- Errors occur if the maximum number of encounter groups or the maximum group size
  set in `symba5p/symba5p.inc` is too small for the simulation.
- Energy offsets due to mergers are not tracked by SyMBAp (flag L4 can be on).

## Input files

Units must make the gravitational constant 1. For example, `example/pl.in` uses AU
and days (1 yr = 365.25 d), which makes the solar mass about 2.96e-4; with AU and
solar masses instead, an orbit at 1 AU has a period of 2π.

### param.in

```
t0 tstop dt
dtout dtdump
L1 L2 L3 L4 L5 L6
rmin rmax rmaxu qmin lclose
binary_outputfile
status_flag_for_open_statements
```

- `t0`, `tstop`, `dt`: start time, stop time, time step.
- `dtout`: interval between outputs; `dtdump`: interval between dumps (setting
  `dtdump = dtout` makes restarts easiest).
- `L1`–`L6` (`T`/`F`):
  - **L1** include J2 and J4 of the central body (values go in `pl.in`).
  - **L2** check for bodies too close to the Sun or a planet, or too far from the
    Sun. Removed test particles are written to `discard.out` (time, position,
    velocity, `istat`/`rstat` flags and all planetary positions and velocities).
    If L2 is `F`, omit the `rmin ...` line so that the file has only three lines
    before the binary file name.
  - **L3** Jacobi integral of the test particles (only useful for one planet on a
    circular orbit), written to `jacobi.out`.
  - **L4** energy and angular momentum of the massive bodies, written to
    `energy.out`.
  - **L5** heliocentric elements a, e, i, Ω, ω, M of all bodies as REAL*4 into the
    binary output file (read it e.g. with `tools/follow`: positive numbers are test
    particles, negative numbers massive bodies, -3 being the third body).
  - **L6** like L5 but as XDR (see [FXDR](#fxdr)). L5 and L6 cannot both be `T`.
- `rmin`: heliocentric distance below which a test particle is stopped; `rmax`:
  distance beyond which it is stopped; `rmaxu`: distance beyond which an *unbound*
  particle is stopped; `qmin`: perihelion distance below which it is stopped. Use
  `-1` to skip a check. The elements are osculating heliocentric values, which are
  unrepresentative beyond a few hundred AU because of the solar reflex motion; use
  `rmax` of a few hundred AU if `qmin` is set.
- `lclose`: if `T`, a test particle closer than `rpl` to a planet is stopped (`rpl`
  in `pl.in`; with WHM use at least a Hill radius). In SyMBA/SyMBAp, bodies closer
  than the sum of their radii merge.
- `binary_outputfile`: name of the binary output file (written every `dtout`).
- `status_flag_for_open_statements`: `new` (stop if the file exists), `append`, or
  `unknown` (overwrite); used for all output files except `lyap.dat`.

### pl.in

First line: number of bodies including the Sun. Then for each body one line with
the mass, one with heliocentric x, y, z and one with vx, vy, vz. The Sun's position
and velocity must be 0.

- If L1 is `T`, the Sun's mass line also holds J2·R² and J4·R⁴ (R = radius of the
  central body), see [Oblateness](#oblateness-j2-j4).
- If L2 and `lclose` are `T`, each planet's mass line also holds its stopping radius
  `rpl`.
- For SyMBA/SyMBAp each planet's mass line holds `mass rhill` (Hill radius), or
  `mass rhill rpl` if `lclose` is `T`.

### tp.in

First line: number of test particles. Then for each particle heliocentric x, y, z;
vx, vy, vz; the 13 `istat` integers and the 13 `rstat` reals (see
[Status flags](#test-particle-status-flags)); start with all flags 0 unless
resuming. The code may change the flags as it runs.

## Output files

- Every `dtout`: the binary output file (L5/L6), `jacobi.out` (L3), `energy.out`
  (L4); SyMBA/SyMBAp also write the masses to `mass.<binary_outputfile>`.
- Every `dtdump`: `dump_pl.dat`, `dump_tp.dat` and `dump_param.dat`, with the same
  formats as `pl.in`, `tp.in` and `param.in`. `t0` in `dump_param.dat` is the dump
  time and the open status is changed to `append`, so the dumps restart a stopped
  run (increase `tstop` to extend it). Only the latest dump is kept.
- `discard.out` (test particles, L2) or `discard_mass.out` (SyMBA/SyMBAp:
  discarded and merged bodies). On restart with `append`, SyMBA/SyMBAp remove
  entries after the restart time from `discard_mass.out`.

Note for developers: for efficiency, some subroutines that only need to know
whether `istat(i,1)` is zero dimension `istat` as a vector rather than a 2D array;
this is noted in each such subroutine.

## Test-particle status flags

`istat` (integer) and `rstat` (real) are dimensioned `(NTPMAX, NSTAT)`, set in
`swift.inc`.

`istat(i,1) = 0`: particle active, and `istat(i,2)` is
- `0` if not in an encounter,
- `+n` if in the outer region of planet n,
- `-n` if in the inner region of planet n.

`istat(i,1) = 1`: particle inactive, and `istat(i,2)` gives the reason:

| `istat(i,2)` | Reason |
|---|---|
| -1 | Danby did not converge |
| -2 | a < 0 and r_helio > rmaxu |
| -3 | r_helio > rmax |
| -4 | q < qmin |
| 1 | r_helio < rmin |
| n | too close to planet n |

`istat(i,j>2)` is otherwise unused. In the RMVS routines `istat(i,3)` is the planet
of the last encounter and `istat(i,m>3)` the number of encounters (within a Hill
sphere) of particle i with planet k = m-2.

`rstat`:

- Active: `rstat(i,2)` is the perihelion distance if it is checked.
- Inactive: `rstat(i,1)` is the discard time. If `istat(i,2) = 1`, `rstat(i,2)` is
  0 if r < rmin, or the perihelion distance if q < qmin. If the particle was
  discarded near planet `istat(i,2)`, `rstat(i,2) > 0` is the closest approach
  predicted by the extrapolation and `rstat(i,2) < 0` the pericentre distance.
- RMVS: `rstat(i,3)` is the pericentre distance to a planet if the particle is in
  the inner encounter region and passed pericentre in the last step;
  `rstat(i,m>3)` is the smallest pericentre distance to planet k = m-2 (used only if
  q < r_Hill).
- LYAP2: `rstat(i,3..5)` and `rstat(i,6..8)` hold dx(1..3) and dv(1..3) of the
  difference equations.

## Oblateness (J2, J4)

The potential of the central body is

    V = -(GM/r) (1 - j2rp2 P2/r² - j4rp4 P4/r⁴)

with `j2rp2 = J2 Rp²`, `j4rp4 = J4 Rp⁴`, Rp the radius of the central body, and the
Legendre polynomials P_n evaluated at z²/r². Some authors use plus signs in this
sum, so check the signs of J2 and J4: in this convention J2 of an oblate body is
positive and advances the argument of perihelion (as for Mercury around an oblate
Sun). `j2rp2` and `j4rp4` must use the same length unit as the positions.

## Fortran unit numbers

| Routine | Unit |
|---|---|
| `io_dump_param`, `io_dump_pl`, `io_dump_tp` | 7 |
| `apple_init`, `io_init_param`, `io_init_pl`, `io_init_tp`, `io_lyap_init` | 7 |
| binary output file | 10 |
| `io_write_frame` | 20 |
| `io_jacobi_write` | 30 |
| `io_discard_write` | 40 |
| `io_lyap_init` | 50 |
| `io_energy_write` | 60 |

## SyMBAp changelog

- **Nov/20, 1.1**: fixed `symba5_getacch.f` passing `ielst` from `integer*2` to
  `integer`.
- **Feb/21, 1.4**: progress output shows particle numbers up to 6 digits.
- **Nov/21, 1.5**: thread count depends on the user-set lower limit of the
  particle-to-thread ratio.
- **Dec/21**: energy checking found incompatible with SyMBAp; turn off L4 until the
  next version.
- **Jan/22, 1.6**: energy offsets due to mergers are not tracked; L4 can be on again.
- **Jul/22**: `io_discard_cleanup` cleans `discard_mass.out` on resumption (SyMBA and
  SyMBAp).
- **Jul/23, 1.7**: automatic collapse of non-rectangular loop nests (required
  GCC 13+).
- **Dec/23, 1.8**: x, y, z arrays merged into `(3,n)` arrays; `-ftree-vectorize`
  added.
- **2025**: array versions of the MVS routines; SIGTERM handler; GNU Make replaces
  the `@make*` scripts; thread count tuned by timing (`symba5p_tune`).
- **Nov/26, 1.9**: pair loops rewritten for speed (results agree with 1.8 to
  round-off): mutual accelerations in `symba5p_pairacc` (vectorised row loop,
  reductions sized `nbod` instead of `NTPMAX`), encounter check inlined and
  vectorised in `symba5p_step_pl`, per-step reset of the group counters removed. No
  loop collapse is used any more, so GCC 13+ is no longer required (the
  `pre_gcc13_array` branch is retired; this is the `array` branch, formerly
  `gcc13_array`). `threads.dat` holding `-n` fixes the thread count. The thread
  input of version ≤ 1.6 (particle-to-thread ratio, then max. threads) is read
  again, so old input files work unchanged. `swift_rmvs3`/`swift_rmvs4` build with
  the medium code model.

## Authors and license

SWIFT: Hal Levison and Martin Duncan (user's manual June 1993). SyMBAp
parallelisation: Tommy Chi Ho Lau and Man Hoi Lee (2020–2026). See
`LICENSE.txt`.
