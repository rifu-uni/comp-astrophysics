# Project 2: MHD Shock-Cloud Interaction

This project reproduces the PLUTO `Test_Problems/MHD/Shock_Cloud/`
benchmark and compares four production PLUTO variants against two independent
Python finite-volume solvers.

The physical problem is the Dai & Woodward shock-cloud setup: a strong
magnetized shock interacts with an overdense cloud on a 2D Cartesian grid,
while the velocity and magnetic field retain all three vector components.
The project runs four PLUTO variants and two Python reference solvers:

| Run | Role |
| --- | --- |
| `01` | Baseline PLUTO setup: RK2, linear reconstruction, constrained transport, no entropy fix. |
| `02` | Same physical setup, but with characteristic tracing instead of RK2. |
| `09` | Same physical setup, but with the SELECTIVE Harten-Hyman entropy fix enabled. |
| `02+09` | Local compound setup combining characteristic tracing with the entropy fix. |
| `python` | Baseline educational solver: cell-centered MUSCL + Rusanov/HLL + RK2 with positivity repair. |
| `python-ct` | CT-upgraded Python solver: same fluid update, but with face-centered `Bx`/`By` and constrained-transport induction. |

The scientific discussion lives in `report/report.tex` and `report/report.pdf`.
Implementation details, setup-selection rationale, and Python solver notes live
in `report/technical-details.md`.

## What Gets Produced

`make all` regenerates the full analysis stack:

- PLUTO outputs for setups `01`, `02`, `09`, and `02+09`.
- Baseline and CT Python outputs at the same `256 x 256`, `tstop = 0.06` target.
- Per-variable summary plots for each PLUTO setup and for both Python solvers.
- PLUTO-vs-PLUTO numerical-scheme comparisons.
- PLUTO-vs-baseline-Python and PLUTO-vs-CT-Python overlays, difference maps,
  and L1/L2/Linf norms.
- Baseline-vs-CT Python overlays, difference maps, and L1/L2/Linf norms.
- Density animations for PLUTO and both Python solvers.

The main interpretation is:

- PLUTO setups differ by only a few percent, mostly near the strong shock.
- The baseline Python solver captures the broad shock envelope but is not a
  quantitative substitute for PLUTO.
- The CT Python solver materially improves magnetic divergence and reduces the
  strongest magnetic/pressure artifacts, but the remaining Rusanov diffusion,
  positivity repair, and simple boundary treatment still leave it below PLUTO.

## Repository Layout

```text
projects/02-shockwave/
├── README.md
├── Makefile
├── inputs/
│   ├── init.c
│   ├── definitions_01.h
│   ├── definitions_02.h
│   ├── definitions_09.h
│   ├── definitions_02+09.h
│   ├── pluto_01.ini
│   ├── pluto_02.ini
│   ├── pluto_09.ini
│   └── pluto_02+09.ini
├── python/
│   ├── solver.py
│   ├── solver_ct.py
│   └── postprocess.py
├── scripts/
│   ├── plot_pluto.py
│   ├── compare_schemes.py
│   └── animate.py
├── compare/
│   ├── overlay.py
│   ├── norms.py
│   └── python_solvers.py
├── report/
│   ├── report.tex
│   ├── report.pdf
│   ├── technical-details.md
│   └── external-figures/
├── build/       # generated PLUTO build directories, gitignored
├── output/      # generated simulation and comparison data, gitignored
└── plots/       # generated figures and animations, gitignored
```

`inputs/` contains the tracked PLUTO problem files used by the Makefile. The
generated directories `build/`, `output/`, and `plots/` are intentionally not
tracked.

## Requirements

- PLUTO source checkout at `../../tools/pluto` relative to this project.
- `pixi` environment with Python, NumPy, Matplotlib, Pillow, and pyPLUTO.
- A C compiler supported by the local PLUTO setup.

The Makefile uses `pixi run python` for Python-side tools and PLUTO's own
`setup.py` for the C build.

## Quick Start

```bash
# Build/run all PLUTO setups, run Python, plot, compare, and animate.
make all
```

Useful individual targets:

```bash
make help

make pluto          # default PLUTO setup, currently 01
make pluto-all      # setups 01, 02, 09, and 02+09
make pluto-01
make pluto-02
make pluto-09
make pluto-02+09

make python         # baseline independent Python MHD solver
make python-ct      # CT-upgraded Python MHD solver

make plot-all       # pyPLUTO plots for all PLUTO setups
make plot-python    # baseline Python solver plots
make plot-python-ct # CT Python solver plots

make compare-schemes # PLUTO-vs-PLUTO norms and density differences
make compare         # PLUTO 01 vs Python overlays and norms
make compare-02      # PLUTO 02 vs Python overlays and norms
make compare-09      # PLUTO 09 vs Python overlays and norms
make compare-02+09   # PLUTO 02+09 vs Python overlays and norms
make compare-python-all     # all four PLUTO setups vs baseline Python
make compare-python-solvers # baseline Python vs CT Python
make compare-ct-all         # all four PLUTO setups vs CT Python

make animate-python
make animate-python-ct
make animate-all     # density GIFs for PLUTO and both Python solvers
make clean           # remove build/, output/, and plots/
```

The full run is idempotent: each PLUTO setup and each Python solver write a
small sentinel file under `output/`, so downstream plot/comparison targets do
not rerun expensive simulations unless the sentinel is removed or the source
changes.

## Outputs

Important generated products:

```text
output/pluto/01/
output/pluto/02/
output/pluto/09/
output/pluto/02+09/
output/python/data.npz
output/python-ct/data.npz
output/compare-schemes/norms.txt
output/compare/norms.txt
output/compare-02/norms.txt
output/compare-09/norms.txt
output/compare-02+09/norms.txt
output/compare-python-solvers/norms.txt
output/compare-ct-01/norms.txt
output/compare-ct-02/norms.txt
output/compare-ct-09/norms.txt
output/compare-ct-02+09/norms.txt
plots/01/summary.png
plots/02/summary.png
plots/09/summary.png
plots/02+09/summary.png
plots/python/summary.png
plots/python-ct/summary.png
plots/animations/
```

The comparison figures use the same final time and grid for PLUTO and Python.
The baseline solver remains the transparent control; the CT solver is the
preferred upgrade path for isolating divergence control without changing the
Rusanov/MUSCL/RK2 fluid update.

## Documentation Map

- `report/report.tex`: scientific report, literature context, simulation
  description, results, and interpretation.
- `report/report.pdf`: compiled version of the scientific report.
- `report/technical-details.md`: implementation guide for the PLUTO setup
  selection, Makefile workflow, Python solver, comparison scripts, and known
  numerical limitations.
- `README.md`: this overview and run guide.

## PLUTO Fork Notes

This project expects a PLUTO checkout compatible with the course branch used
for the project. The local build recipe stages `definitions_*.h`, `pluto_*.ini`,
and `init.c` into `build/<setup>/`, runs PLUTO setup in noninteractive mode,
builds the executable, and copies the resulting `.dbl` outputs into
`output/pluto/<setup>/`.
