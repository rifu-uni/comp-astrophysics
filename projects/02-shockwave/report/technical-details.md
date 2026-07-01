# Technical Details: MHD Shock-Cloud Project

This document is the tracked implementation companion to the scientific report.
It collects the setup-selection rationale, PLUTO build organization, Python
solver design, comparison workflow, and known numerical limitations.

The scientific narrative stays in `report.tex`; this file is for the mechanics.

## 1. PLUTO Problem and Shared Initial Conditions

The project uses PLUTO's `Test_Problems/MHD/Shock_Cloud/` benchmark. The
upstream folder ships nine configurations that share the same physical initial
condition and runtime parameters, but differ in numerical scheme.

The initial condition is the Dai & Woodward shock-cloud setup. The shock is
initially at `x = 0.6`. The cloud has radius `RADIUS = 0.15` and center
`(0.8, 0.5)` in the 2D run. PLUTO keeps three vector components even though the
grid has two spatial dimensions.

| Region | rho | vx | p | By | Bz |
| --- | ---: | ---: | ---: | ---: | ---: |
| Post-shock, `x < 0.6` | `3.86859` | `0.0` | `167.345` | `+B_POST` | `-B_POST` |
| Pre-shock, `x >= 0.6` | `1.0` | `-11.2536` | `1.0` | `+B_PRE` | `+B_PRE` |
| Cloud | `10.0` | inherited | inherited | inherited | inherited |

with:

```text
B_POST = 2.1826182
B_PRE  = 0.56418958
```

The boundary conditions are:

| Boundary | PLUTO condition | Purpose |
| --- | --- | --- |
| `X1-beg`, `x = 0` | `outflow` | Let post-shock material leave with minimal reflection. |
| `X1-end`, `x = 1` | `userdef` | Reset the inflow to the pre-shock state. |
| `X2-beg`, `y = 0` | `outflow` | Let transverse material leave. |
| `X2-end`, `y = 1` | `outflow` | Let transverse material leave. |
| `X3-beg/end` | `outflow` | Present because PLUTO keeps one inactive z cell. |

The runtime grid and time settings are:

```text
X1-grid = 0.0 .. 1.0, 256 cells
X2-grid = 0.0 .. 1.0, 256 cells
X3-grid = 0.0 .. 1.0, 1 cell
CFL     = 0.4
tstop   = 0.06
first_dt = 1e-5
Solver  = roe
```

The `[Chombo Refinement]` block appears in the upstream `.ini` files and is
kept for fidelity, but it is inert in the standard non-Chombo build used here.

## 2. Upstream Shock_Cloud Setups

All upstream setups satisfy the "MHD, at least 2D" requirement at the physics
level, but only the low-cost 2D variants are appropriate for a course project
that also needs an independent Python reproduction.

| Setup | Dim | Static cells | Reconstruction | Time step | divB control | Entropy | Solver | AMR | Role |
| --- | --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `01` | 2D | `256 x 256 x 1` | `LINEAR` | `RK2` | `CT` | `NO` | `roe` | no | Canonical baseline. |
| `02` | 2D | `256 x 256 x 1` | `LINEAR` | `ChTr` | `CT` | `NO` | `roe` | no | Time-stepping comparison. |
| `03` | 2D | `256 x 256 x 1` | `LINEAR` | `HANCOCK` | `8W` | `NO` | `roe` | no | Divergence-control comparison. |
| `04` | 3D | `96 x 96 x 96` | `LINEAR` | `ChTr` | `CT` | `NO` | `roe` | no | 3D characteristic tracing. |
| `05` | 3D | `96 x 96 x 96` | `LINEAR` | `RK3` | `CT` | `NO` | `roe` | no | 3D third-order RK. |
| `06` | 3D | `96 x 96 x 96` | `LINEAR` | `ChTr` | `GLM` | `NO` | `roe` | no | 3D GLM divergence cleaning. |
| `07` | 3D | `64 x 64 x 64` | `WENO3_FD` | `RK3` | `GLM` | `NO` | `roe` | no | WENO3 + GLM. |
| `08` | 3D | `16 x 16 x 16` base | `LINEAR` | `HANCOCK` | `GLM` | `NO` | `hll` | yes | Chombo AMR setup. |
| `09` | 2D | `256 x 256 x 1` | `LINEAR` | `RK2` | `CT` | `SELECTIVE` | `roe` | no | Entropy-fix comparison. |

`ChTr` means `CHARACTERISTIC_TRACING`. `CT` means constrained transport.
`8W` means Powell's eight-wave formulation. `GLM` means Dedner-style
hyperbolic divergence cleaning. `SELECTIVE` enables the entropy fix in
troubled cells near shocks.

Setups `04` through `08` are skipped because a 3D Python reproduction would
substantially exceed the project scope. Setup `08` also needs a Chombo-enabled
PLUTO build.

## 3. Chosen Comparison Set

The project runs setups `01`, `02`, `09`, and a local compound setup `02+09`.
They all share the same grid, initial condition, runtime parameters, and
boundary conditions.

| Comparison | Isolated effect | What it probes |
| --- | --- | --- |
| `01` vs `02` | `RK2` vs `CHARACTERISTIC_TRACING` | Time integration and characteristic update effects. |
| `01` vs `09` | `ENTROPY_SWITCH NO` vs `SELECTIVE` | Harten-Hyman entropy fix near strong shocks. |
| `02` vs `02+09` | Entropy fix on the ChTr branch | Same entropy-fix probe after changing the time stepper. |
| `09` vs `02+09` | RK2 vs ChTr with entropy fix enabled | Time-stepping effect after stabilizing the shock. |

The local `02+09` setup is not a new physical problem. It combines setup `02`'s
characteristic tracing with setup `09`'s entropy fix and Van Leer limiter
constant. This gives a clean two-knob extension after the one-knob comparisons.

Setup `03` would be a useful divergence-control comparison, but it is left out
of the default set because mixing CT and eight-wave divergence control makes
the error analysis harder to interpret in a short report.

## 4. Tracked PLUTO Inputs

The tracked inputs are:

```text
inputs/
├── init.c
├── definitions_01.h
├── definitions_02.h
├── definitions_09.h
├── definitions_02+09.h
├── pluto_01.ini
├── pluto_02.ini
├── pluto_09.ini
└── pluto_02+09.ini
```

The `pluto_*.ini` files are intentionally almost identical. The numerical
differences are compile-time PLUTO choices in `definitions_*.h`:

| File | Time stepping | divB control | Entropy switch | Extra limiter |
| --- | --- | --- | --- | --- |
| `definitions_01.h` | `RK2` | `CONSTRAINED_TRANSPORT` | `NO` | none |
| `definitions_02.h` | `CHARACTERISTIC_TRACING` | `CONSTRAINED_TRANSPORT` | `NO` | none |
| `definitions_09.h` | `RK2` | `CONSTRAINED_TRANSPORT` | `SELECTIVE` | `VANLEER_LIM` |
| `definitions_02+09.h` | `CHARACTERISTIC_TRACING` | `CONSTRAINED_TRANSPORT` | `SELECTIVE` | `VANLEER_LIM` |

All four use:

```c
#define PHYSICS          MHD
#define DIMENSIONS       2
#define COMPONENTS       3
#define GEOMETRY         CARTESIAN
#define RECONSTRUCTION   LINEAR
#define EOS              IDEAL
```

## 5. Makefile Workflow

The Makefile expects PLUTO at:

```text
../../tools/pluto
```

Each `make pluto-NN` target stages a clean build directory:

```text
build/NN/
├── pluto.ini       # copied from inputs/pluto_NN.ini
├── definitions.h   # copied from inputs/definitions_NN.h
└── init.c          # shared initial condition
```

Then it writes a small noninteractive makefile stub, runs PLUTO's
`setup.py --auto-update`, builds the executable, runs `./pluto -i pluto.ini`,
copies `.dbl`, `dbl.out`, and `grid.out` into `output/pluto/NN/`, and writes
`output/pluto/NN/done` as a sentinel.

The sentinel files are deliberate. Plotting and comparison targets depend on
them so repeated `make all` invocations do not rerun finished simulations.

The full analysis target is:

```bash
make all
```

which expands to the four PLUTO runs/plots, both Python solver plots, the
PLUTO-vs-PLUTO comparison, all four PLUTO-vs-baseline-Python comparisons, the
baseline-vs-CT Python comparison, all four PLUTO-vs-CT comparisons, and density
animations.

Useful Python-side targets are:

| Target | Output |
| --- | --- |
| `make python` | `output/python/data.npz` from `python/solver.py`. |
| `make python-ct` | `output/python-ct/data.npz` from `python/solver_ct.py`. |
| `make plot-python` | `plots/python/summary.png` and per-variable panels. |
| `make plot-python-ct` | `plots/python-ct/summary.png` and per-variable panels. |
| `make compare-python-all` | `output/compare*` norms/plots for all PLUTO setups vs baseline Python. |
| `make compare-python-solvers` | `output/compare-python-solvers/` baseline-vs-CT norms and diff plots. |
| `make compare-ct-all` | `output/compare-ct-*` norms/plots for all PLUTO setups vs CT Python. |
| `make animate-python` | `plots/animations/python_rho.gif`. |
| `make animate-python-ct` | `plots/animations/python_ct_rho.gif`. |

## 6. Python Solver Overview

There are now two independent Python solvers. They share the same CLI/output
contract and the same physical initial/boundary conditions, so all existing
plotting, animation, and PLUTO comparison scripts can consume either `.npz`
file.

| Solver | Module | Role |
| --- | --- | --- |
| Baseline | `python/solver.py` | Cell-centered MUSCL + Rusanov/HLL + RK2 control case. |
| CT upgrade | `python/solver_ct.py` | Same fluid update, but with face-centered `Bx`/`By` and constrained transport. |

Shared baseline choices:

| Component | Implementation |
| --- | --- |
| Physics | Ideal MHD with eight conserved variables. |
| Geometry | 2D Cartesian grid with three velocity and magnetic components. |
| Grid | Cell-centered, default `256 x 256`, domain `[0,1] x [0,1]`. |
| Reconstruction | MUSCL piecewise-linear reconstruction. |
| Limiter | Minmod. |
| Riemann solver | HLL with Einfeldt/Rusanov wave speeds. |
| Time stepping | RK2 / Heun method. |
| Timestep | Adaptive CFL based on the maximum fast magnetosonic speed. |
| Boundary conditions | Same intent as PLUTO: outflow left/top/bottom and pre-shock user state at right. |
| Stabilization | Interface fallback plus density/pressure floors. |

The default commands used by the Makefile are:

```bash
pixi run python -m python.solver \
  --nx 256 --ny 256 \
  --tstop 0.06 \
  --cfl 0.3 \
  --snapshot-count 21 \
  --snapshot-vars rho \
  --out output/python/data.npz

pixi run python -m python.solver_ct \
  --nx 256 --ny 256 \
  --tstop 0.06 \
  --cfl 0.3 \
  --snapshot-count 21 \
  --snapshot-vars rho \
  --out output/python-ct/data.npz
```

Both solvers write:

```text
x1, x2, t,
rho, vx, vy, vz, prs, Bx, By, Bz,
snapshot_t, <var>_series   # when snapshots are requested
```

The CT solver additionally writes scalar diagnostics:

```text
divB_max, divB_mean,
pressure_floor_cells, pressure_floor_repairs,
nsteps, wall_s
```

`divB_max` and `divB_mean` are the discrete face-divergence diagnostics from
the CT mesh, not a cell-centered finite-difference estimate from the output
arrays.

## 7. Python Initial Condition

The Python initial condition mirrors `inputs/init.c`:

1. Cell centers are built on `[0,1] x [0,1]`.
2. Cells with `x < 0.6` receive the post-shock state.
3. Cells with `x >= 0.6` receive the pre-shock state.
4. The spherical cloud is represented as a circle in the 2D plane with center
   `(0.8, 0.5)` and radius `0.15`.
5. Cloud cells get `rho = 10.0` while inheriting the local pressure, velocity,
   and magnetic field state.

The primitive variables are:

```text
rho, vx, vy, vz, p, Bx, By, Bz
```

They are converted to conserved form using:

```text
E = p / (gamma - 1)
    + 0.5 * rho * (vx^2 + vy^2 + vz^2)
    + 0.5 * (Bx^2 + By^2 + Bz^2)
```

with `gamma = 5/3`.

## 8. Fluxes and Wave Speeds

The solver computes ideal-MHD conservative fluxes in either the x or y
direction. The HLL flux uses symmetric Einfeldt/Rusanov speeds:

```text
Smax = max(|u_n,L| + cf,L, |u_n,R| + cf,R)
SL   = -Smax
SR   = +Smax
```

where `u_n` is the velocity normal to the interface and `cf` is the fast
magnetosonic speed:

```text
cf^2 = 0.5 * (cs^2 + ca^2
       + sqrt((cs^2 + ca^2)^2 - 4 * cs^2 * Bn^2 / rho))
```

This is more diffusive than Roe or HLLD, but it is robust enough for this
strong-shock problem at moderate resolution.

## 9. MUSCL Reconstruction and Positivity Repair

For each direction, the solver:

1. Pads the conserved state with one ghost cell.
2. Computes one-sided differences.
3. Applies the minmod limiter.
4. Reconstructs left and right interface states.
5. Replaces nonphysical reconstructed states with first-order states.
6. Computes HLL fluxes.
7. Forms the finite-volume flux divergence.

The first-order fallback is applied interface-wise when either reconstructed
side would have nonfinite or floor-violating density/pressure.

After each RK stage, the solver enforces:

```text
RHO_FLOOR = 1e-3
PRS_FLOOR = 1e-5
```

If pressure falls below the floor, the total energy is reset to the minimum
energy consistent with the floored pressure plus the existing kinetic and
magnetic energy. This keeps the run finite, but it is not strictly
conservative. That limitation is part of the reported PLUTO-vs-Python error.

## 10. Boundary Conditions in Python

The ghost-cell fill is:

| Side | Python behavior |
| --- | --- |
| Left, `x = 0` | Copy the first interior column, matching outflow/zero-gradient. |
| Right, `x = 1` | Fill ghost cells with the pre-shock conservative state. |
| Bottom, `y = 0` | Copy the first interior row, matching outflow/zero-gradient. |
| Top, `y = 1` | Copy the last interior row, matching outflow/zero-gradient. |

The right boundary matters because the pre-shock inflow must remain stable as
the shock-cloud interaction evolves.

## 11. RK2 Update and Snapshots

The RK2 update is:

```text
k1 = -RHS(U)
U1 = floor(U + dt*k1)
k2 = -RHS(U1)
U  = floor(U + 0.5*dt*(k1 + k2))
```

The timestep is recomputed after each full step:

```text
dt = CFL * min(dx, dy) / max(|v| + cf)
```

The output `.npz` contains:

```text
x1, x2, t,
rho, vx, vy, vz, prs, Bx, By, Bz
```

When `--snapshot-count` is set, it also stores `snapshot_t` and
`<variable>_series` arrays for animation.

## 12. CT Solver Internals

`python/solver_ct.py` keeps the conserved fluid variables cell-centered but
stores the transverse magnetic field on a staggered mesh:

```text
Bx_face: shape (NX + 1, NY)   # x-faces
By_face: shape (NX, NY + 1)   # y-faces
```

For fluxes and output, the solver reconstructs cell-centered fields by face
averaging:

```text
Bx_cell[i,j] = 0.5 * (Bx_face[i,j] + Bx_face[i+1,j])
By_cell[i,j] = 0.5 * (By_face[i,j] + By_face[i,j+1])
```

The induction update uses the edge-centered electric field

```text
Ez = vx*By - vy*Bx
```

with local Lax-Friedrichs upwinding around the arithmetic corner average. The
face update is the discrete curl:

```text
dBx_face/dt = -(Ez[:, 1:] - Ez[:, :-1]) / dy
dBy_face/dt = +(Ez[1:, :] - Ez[:-1, :]) / dx
```

Because both face components are updated by the same edge-centered curl, the
discrete face divergence

```text
(Bx_face[1:, :] - Bx_face[:-1, :]) / dx
+ (By_face[:, 1:] - By_face[:, :-1]) / dy
```

is conserved to round-off. At the end of each RK stage, the face fields are
averaged back to cell centers and copied into the conservative state. The total
energy is corrected by the magnetic-energy difference during that projection,
so the CT synchronization preserves the cell's kinetic plus thermal energy
rather than artificially draining pressure.

The CT solver deliberately does not add HLLD, Roe, entropy fixes, AMR, or
higher-order reconstruction. This isolates divergence control as the upgrade.

### GLM fallback criteria

GLM cleaning was reserved as a fallback only if the CT implementation could
not reach the production run after conservative tuning. The fallback criteria
were:

- CT fails the `32 x 32`, `tstop=0.005` smoke test with finite positive fields.
- CT fails the `64 x 64`, `tstop=0.02`, `snapshot-count=5` smoke test.
- CT cannot reach `256 x 256`, `tstop=0.06` without runaway wave speeds after
  adding magnetic-energy synchronization and upwind EMF dissipation.

The final CT solver passes all three criteria, so no GLM implementation is
used in the deliverable.

## 13. Plotting and Comparison Scripts

PLUTO plots:

```bash
pixi run python scripts/plot_pluto.py --setup 01 --out plots/01
```

Python plots:

```bash
pixi run python -m python.postprocess \
  --in output/python/data.npz \
  --out plots/python

pixi run python -m python.postprocess \
  --in output/python-ct/data.npz \
  --out plots/python-ct
```

PLUTO-vs-PLUTO comparison:

```bash
pixi run python scripts/compare_schemes.py \
  --out output/compare-schemes
```

This compares all pairs among `01`, `02`, `09`, and `02+09`, writing:

```text
output/compare-schemes/norms.txt
output/compare-schemes/density.png
output/compare-schemes/diff_01_02.png
output/compare-schemes/diff_01_09.png
output/compare-schemes/diff_02_02_09.png
output/compare-schemes/diff_09_02_09.png
```

PLUTO-vs-Python overlays and norms:

```bash
pixi run python -m compare.overlay --setup 01 --out output/compare
pixi run python -m compare.norms   --setup 01 --out output/compare

pixi run python -m compare.overlay --setup 01 \
  --python-in output/python-ct/data.npz \
  --out output/compare-ct-01
pixi run python -m compare.norms --setup 01 \
  --python-in output/python-ct/data.npz \
  --out output/compare-ct-01
```

The comparison assumes PLUTO and Python are on the same final grid. If a shape
or grid mismatch appears, the plotting script warns instead of silently
interpolating.

Python-vs-Python comparison:

```bash
pixi run python -m compare.python_solvers \
  --baseline output/python/data.npz \
  --candidate output/python-ct/data.npz \
  --out output/compare-python-solvers
```

Animations:

```bash
pixi run python scripts/animate.py --source python \
  --python-in output/python/data.npz \
  --out plots/animations/python_rho.gif

pixi run python scripts/animate.py --source python \
  --python-in output/python-ct/data.npz \
  --title "Python CT solver" \
  --out plots/animations/python_ct_rho.gif
```

## 14. Current Numerical Results

The PLUTO-vs-PLUTO scheme spread is small compared with the PLUTO-vs-Python
solver-family gap.

Representative PLUTO-vs-PLUTO density and pressure norms:

| Pair | L1(rho) | L2(rho) | L1(p) | L2(p) | Isolated factor |
| --- | ---: | ---: | ---: | ---: | --- |
| `01-02` | `0.0636` | `0.332` | `1.09` | `3.26` | Time stepping. |
| `01-09` | `0.0977` | `0.461` | `1.78` | `4.95` | Entropy fix. |
| `02-02+09` | `0.0859` | `0.470` | `1.56` | `4.28` | Entropy fix on ChTr branch. |
| `09-02+09` | `0.0466` | `0.228` | `0.993` | `2.46` | Time stepping with entropy fix. |

Final Python diagnostics at `256 x 256`, `t=0.06`:

| Diagnostic | Baseline Python | CT Python |
| --- | ---: | ---: |
| Total mass | `3.9661` | `3.9653` |
| Pressure-floor cells | `188` | `176` |
| `Bx` range | `[-76.46, 76.46]` | `[-14.53, 14.53]` |
| Cell-centered `max|divB|` | `9.63e3` | `7.14e1` |
| Face-CT `max|divB|` | n/a | `6.80e-12` |
| Face-CT mean `|divB|` | n/a | `5.01e-13` |
| Steps | not saved | `841` |

The cell-centered divergence row is computed from the public output arrays.
The face-CT rows are the physically relevant CT diagnostics saved by
`solver_ct.py`.

Representative PLUTO-vs-Python density and pressure norms:

| PLUTO setup | Python solver | L1(rho) | L2(rho) | L1(p) | L2(p) |
| --- | --- | ---: | ---: | ---: | ---: |
| `01` | baseline | `0.611` | `2.39` | `15.95` | `42.66` |
| `02` | baseline | `0.607` | `2.39` | `15.51` | `42.12` |
| `09` | baseline | `0.625` | `2.46` | `16.16` | `42.85` |
| `02+09` | baseline | `0.621` | `2.45` | `15.98` | `42.71` |
| `01` | CT | `0.250` | `0.956` | `6.87` | `16.54` |
| `02` | CT | `0.247` | `0.946` | `6.62` | `15.72` |
| `09` | CT | `0.229` | `0.834` | `6.89` | `16.35` |
| `02+09` | CT | `0.232` | `0.838` | `6.70` | `15.99` |

Baseline-vs-CT Python norms:

| Variable | L1 | L2 | Linf |
| --- | ---: | ---: | ---: |
| `rho` | `0.570` | `2.19` | `45.3` |
| `Bx` | `1.84` | `4.94` | `78.7` |
| `By` | `1.74` | `4.39` | `65.6` |
| `prs` | `15.45` | `41.18` | `1134` |

The CT solver is not a PLUTO clone, but it substantially reduces the dominant
magnetic-divergence failure mode and brings the PLUTO-vs-Python density and
pressure norms down by roughly a factor of two to three.

## 15. Known Limitations

PLUTO limitations in this project:

- The resolution is only `256 x 256`, enough for morphology but not for a
  converged turbulent wake study.
- Only two numerical axes are explored by default: time stepping and entropy
  fix.
- The 3D and AMR upstream configurations are intentionally out of scope.

Python solver limitations:

- The baseline solver has no constrained transport, so discrete `nabla dot B`
  is not machine-zero.
- The CT solver controls face-divergence, but it still uses a simple Rusanov
  fluid update and a low-order upwind EMF rather than a full production CT
  Godunov scheme.
- No GLM or eight-wave divergence cleaning is implemented because CT succeeded.
- No HLLD/Roe contact and Alfven-wave resolution.
- Rusanov wave speeds are stable but strongly diffusive.
- Positivity repair changes total energy locally and is not strictly
  conservative.
- One-cell ghost boundary handling is much simpler than PLUTO's boundary
  machinery.
- The solver is a qualitative independent reproduction, not a calibrated
  production replacement.

These limitations are intentional: the Python solvers demonstrate the core
finite-volume pieces and the value of CT in a transparent way while PLUTO
remains the quantitative reference.

## 16. Corrections Preserved from the Setup Notes

The setup comparison in this file preserves several corrections from the
working notes:

1. The Chombo block is present in upstream `.ini` files, but only setup `08`
   is a Chombo AMR setup.
2. Setup `08` also changes the Riemann solver to `hll`; the other listed
   setups use `roe`.
3. Setup `03` is not "just another 2D setup"; it switches divergence control
   from CT to eight-wave, which changes the interpretation of magnetic errors.
4. The physical benchmark is Dai & Woodward's shock-cloud problem, not a new
   PLUTO-only configuration.
5. The project default is a controlled 2D numerical-method comparison, not a
   broad sweep through every upstream Shock_Cloud variant.

## 17. References for Implementation Choices

- Dai & Woodward (1994): original MHD shock-cloud benchmark.
- Toth (2000): divergence-control comparison in shock-capturing MHD.
- Dedner et al. (2002): GLM divergence cleaning.
- Mignone et al. (2007): PLUTO code paper.
- Toro (2009): approximate Riemann solvers and fast magnetosonic wave speeds.
