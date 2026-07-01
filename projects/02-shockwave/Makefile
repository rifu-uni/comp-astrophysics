# Project 2: MHD Shock-Cloud — Build Automation
# =============================================
# PLUTO_DIR — path to the PLUTO C source (currently at ../tools/pluto/;
# will become a git submodule after fork is pushed)
PU_DIR := $(realpath $(CURDIR)/../../tools/pluto)
ifeq ($(wildcard $(PU_DIR)/Src/main.c),)
  $(error PLUTO source not found at $(PU_DIR). Expected a PLUTO checkout.)
endif

WORKDIR  := $(CURDIR)/build
OUTDIR   := $(CURDIR)/output
PLOTDIR  := $(CURDIR)/plots
ANIMDIR  := $(PLOTDIR)/animations
MPLCONFIGDIR := $(CURDIR)/.mplconfig
export MPLCONFIGDIR

# Use the pixi environment's python (pyPLUTO, numpy, matplotlib, Pillow).
PYTHON   := pixi run python

# Setups used in the numerical-scheme comparison.
# See report/technical-details.md for the rationale.
SETUPS        := 01 02 09 02+09
DEFAULT_SETUP := 01

# ---- default ----
.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "Project 2 — MHD Shock Cloud (setups 01, 02, 09, 02+09)"
	@echo ""
	@echo "Setup comparison (01 vs 02 vs 09 vs 02+09): see report/technical-details.md"
	@echo ""
	@echo "Targets:"
	@echo "  make pluto          — build + run default setup ($(DEFAULT_SETUP))"
	@echo "  make pluto-01       — build + run setup 01 (RK2, LINEAR, CT, no entropy)"
	@echo "  make pluto-02       — build + run setup 02 (ChTr, LINEAR, CT)"
	@echo "  make pluto-09       — build + run setup 09 (RK2, LINEAR, CT, SELECTIVE entropy)"
	@echo "  make pluto-02+09    — build + run compound setup 02+09 (ChTr, LINEAR, CT, SELECTIVE entropy)"
	@echo "  make pluto-all      — build + run setups 01, 02, 09, 02+09"
	@echo ""
	@echo "  make plot-NN        — pyPLUTO plots for setup NN -> plots/NN/"
	@echo "  make plot-all       — pyPLUTO plots for all setups"
	@echo "  make animate-all    — density GIFs for PLUTO setups + Python solvers"
	@echo "  make animate-python — density GIF for the Python solver"
	@echo "  make animate-python-ct — density GIF for the CT Python solver"
	@echo ""
	@echo "  make python         — run the independent Python MHD solver"
	@echo "  make python-ct      — run the CT Python MHD solver"
	@echo "  make compare        — PLUTO 01 vs Python overlay + L1/L2 norms"
	@echo "  make compare-02+09  — PLUTO 02+09 vs Python overlay + L1/L2 norms"
	@echo "  make compare-python-all — PLUTO 01/02/09/02+09 vs baseline Python"
	@echo "  make compare-python-solvers — baseline Python vs CT Python"
	@echo "  make compare-ct-all — PLUTO 01/02/09/02+09 vs CT Python"
	@echo "  make compare-schemes — PLUTO-vs-PLUTO across 01/02/09/02+09"
	@echo ""
	@echo "  make all            — pluto + python + plots + comparisons + animations"
	@echo "  make clean          — wipe build/ output/ plots/"
	@echo ""
	@echo "PLUTO_DIR = $(PU_DIR)"

# ---- per-setup build + run ----
# pluto-NN targets recurse into _pluto-build with SETUP=NN so the
# recipe is defined in exactly one place. A sentinel file per setup
# ($(OUTDIR)/pluto/NN/done) gates the work; downstream targets depend
# on the sentinel, so re-running `make plot-all` or `make compare-schemes`
# does NOT re-trigger PLUTO.
#
# pluto-NN -> $(OUTDIR)/pluto/NN/done is wired through a static pattern
# rule (target-pattern in the rule header), so the `%` is substituted
# with NN for each target. Without a static-pattern header, make would
# interpret `$(OUTDIR)/pluto/%/done` as a literal target name and fail
# with "No rule to make target output/pluto/%/done" on a fresh checkout
# (where no sentinels exist yet).
.PHONY: $(addprefix pluto-,$(SETUPS)) pluto pluto-all _pluto-build
SENTINELS := $(addprefix $(OUTDIR)/pluto/,$(SETUPS))

pluto: pluto-$(DEFAULT_SETUP)
pluto-all: $(addprefix pluto-,$(SETUPS))

# Static pattern rule: for each pluto-NN target, declare the matching
# $(OUTDIR)/pluto/NN/done sentinel as its prereq. The recipe (@:) is a
# no-op so this rule only fires when the sentinel is missing (and the
# recipe below creates it).
$(addprefix pluto-,$(SETUPS)): pluto-%: $(OUTDIR)/pluto/%/done
	@:

# pluto_01.ini == pluto_02.ini == pluto_09.ini == pluto_02+09.ini modulo whitespace
# (the setups differ only in definitions_*.h). One file per setup is
# kept for parity with the upstream PLUTO Shock_Cloud/ folder.
$(OUTDIR)/pluto/%/done: inputs/pluto_%.ini inputs/definitions_%.h inputs/init.c
	@$(MAKE) _pluto-build SETUP=$*
	@touch $@

_pluto-build:
	@echo ">>> Building PLUTO for Shock_Cloud setup $(SETUP)..."
	@rm -rf $(WORKDIR)/$(SETUP)
	@mkdir -p $(WORKDIR)/$(SETUP) $(OUTDIR)/pluto/$(SETUP)
	@cp inputs/pluto_$(SETUP).ini $(WORKDIR)/$(SETUP)/pluto.ini
	@cp inputs/definitions_$(SETUP).h $(WORKDIR)/$(SETUP)/definitions.h
	@cp inputs/init.c $(WORKDIR)/$(SETUP)/
# Pre-stage a makefile stub with ARCH + PLUTO_DIR so PLUTO's setup.py
# bypasses the curses-based ARCH menu (which fails in non-TTY envs).
	@printf 'ARCH      = Darwin.gcc.defs\nPLUTO_DIR = %s\n' "$(PU_DIR)" > $(WORKDIR)/$(SETUP)/makefile
	@cd $(WORKDIR)/$(SETUP) && PLUTO_DIR=$(PU_DIR) python3 $(PU_DIR)/setup.py --auto-update
	@$(MAKE) -C $(WORKDIR)/$(SETUP) PLUTO_DIR=$(PU_DIR)
	@cd $(WORKDIR)/$(SETUP) && ./pluto -i pluto.ini
	@cp $(WORKDIR)/$(SETUP)/data.*.dbl $(WORKDIR)/$(SETUP)/dbl.out $(WORKDIR)/$(SETUP)/grid.out $(OUTDIR)/pluto/$(SETUP)/ 2>/dev/null || true
	@echo ">>> Setup $(SETUP) done: $(OUTDIR)/pluto/$(SETUP)/"

# ---- per-setup pyPLUTO plots ----
.PHONY: plot-01 plot-02 plot-09 plot-02+09 plot-all

# Use explicit targets (GNU make 3.81 on macOS is flaky with
# implicit `plot-%` pattern rules). Each target depends on its setup's
# sentinel, so it will re-run only if the underlying PLUTO data changed.
plot-01: $(OUTDIR)/pluto/01/done
	@echo ">>> pyPLUTO plots for setup 01 -> $(PLOTDIR)/01/"
	@mkdir -p $(PLOTDIR)/01
	@$(PYTHON) scripts/plot_pluto.py --setup 01 --out $(PLOTDIR)/01

plot-02: $(OUTDIR)/pluto/02/done
	@echo ">>> pyPLUTO plots for setup 02 -> $(PLOTDIR)/02/"
	@mkdir -p $(PLOTDIR)/02
	@$(PYTHON) scripts/plot_pluto.py --setup 02 --out $(PLOTDIR)/02

plot-09: $(OUTDIR)/pluto/09/done
	@echo ">>> pyPLUTO plots for setup 09 -> $(PLOTDIR)/09/"
	@mkdir -p $(PLOTDIR)/09
	@$(PYTHON) scripts/plot_pluto.py --setup 09 --out $(PLOTDIR)/09

plot-02+09: $(OUTDIR)/pluto/02+09/done
	@echo ">>> pyPLUTO plots for setup 02+09 -> $(PLOTDIR)/02+09/"
	@mkdir -p $(PLOTDIR)/02+09
	@$(PYTHON) scripts/plot_pluto.py --setup 02+09 --out $(PLOTDIR)/02+09

plot-all: plot-01 plot-02 plot-09 plot-02+09

# ---- PLUTO-vs-PLUTO comparison across setups ----
.PHONY: compare-schemes
compare-schemes: $(addsuffix /done,$(SENTINELS))
	@echo ">>> Cross-scheme comparison (01 vs 02 vs 09 vs 02+09)..."
	@$(PYTHON) scripts/compare_schemes.py --out $(OUTDIR)/compare-schemes
	@echo ">>> Cross-scheme output: $(OUTDIR)/compare-schemes/"

# ---- Python solvers and PLUTO-vs-Python comparison ----
.PHONY: python python-ct plot-python plot-python-ct compare compare-02 compare-09 compare-02+09 compare-python-all \
        compare-python-solvers compare-ct-01 compare-ct-02 compare-ct-09 \
        compare-ct-02+09 compare-ct-all

PY_SENTINEL := $(OUTDIR)/python/done
PY_CT_SENTINEL := $(OUTDIR)/python-ct/done

# The Python solver is independent of PLUTO.  It writes output/python/data.npz.
# Gated by a sentinel so re-runs are idempotent.
python: $(PY_SENTINEL)
	@echo ">>> Python solver output: $(OUTDIR)/python/"

$(PY_SENTINEL): python/solver.py Makefile
	@echo ">>> Running Python MHD solver (Rusanov + MUSCL + RK2 on 256x256)..."
	@mkdir -p $(OUTDIR)/python
	@$(PYTHON) -m python.solver --nx 256 --ny 256 --tstop 0.06 --cfl 0.3 --snapshot-count 21 --snapshot-vars rho --out $(OUTDIR)/python/data.npz
	@touch $@

python-ct: $(PY_CT_SENTINEL)
	@echo ">>> CT Python solver output: $(OUTDIR)/python-ct/"

$(PY_CT_SENTINEL): python/solver_ct.py python/solver.py Makefile
	@echo ">>> Running CT Python MHD solver (Rusanov + MUSCL + RK2 + CT on 256x256)..."
	@mkdir -p $(OUTDIR)/python-ct
	@$(PYTHON) -m python.solver_ct --nx 256 --ny 256 --tstop 0.06 --cfl 0.3 --snapshot-count 21 --snapshot-vars rho --out $(OUTDIR)/python-ct/data.npz
	@touch $@

plot-python: $(PY_SENTINEL)
	@echo ">>> Python postprocess plots -> $(PLOTDIR)/python/"
	@mkdir -p $(PLOTDIR)/python
	@$(PYTHON) -m python.postprocess --in $(OUTDIR)/python/data.npz --out $(PLOTDIR)/python

plot-python-ct: $(PY_CT_SENTINEL)
	@echo ">>> CT Python postprocess plots -> $(PLOTDIR)/python-ct/"
	@mkdir -p $(PLOTDIR)/python-ct
	@$(PYTHON) -m python.postprocess --in $(OUTDIR)/python-ct/data.npz --out $(PLOTDIR)/python-ct

compare-python-solvers: $(PY_SENTINEL) $(PY_CT_SENTINEL)
	@echo ">>> Python baseline vs CT comparison -> $(OUTDIR)/compare-python-solvers/"
	@mkdir -p $(OUTDIR)/compare-python-solvers
	@$(PYTHON) -m compare.python_solvers --baseline $(OUTDIR)/python/data.npz --candidate $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-python-solvers

# PLUTO-vs-Python comparison: overlay + L1/L2 norms.  Depends on PLUTO 01
# and the Python solver output.
compare: $(OUTDIR)/pluto/01/done $(PY_SENTINEL)
	@echo ">>> PLUTO-vs-Python comparison -> $(OUTDIR)/compare/"
	@mkdir -p $(OUTDIR)/compare
	@$(PYTHON) -m compare.overlay --setup 01 --out $(OUTDIR)/compare
	@$(PYTHON) -m compare.norms   --setup 01 --out $(OUTDIR)/compare

compare-02: $(OUTDIR)/pluto/02/done $(PY_SENTINEL)
	@echo ">>> PLUTO 02 vs Python comparison -> $(OUTDIR)/compare-02/"
	@mkdir -p $(OUTDIR)/compare-02
	@$(PYTHON) -m compare.overlay --setup 02 --out $(OUTDIR)/compare-02
	@$(PYTHON) -m compare.norms   --setup 02 --out $(OUTDIR)/compare-02

compare-09: $(OUTDIR)/pluto/09/done $(PY_SENTINEL)
	@echo ">>> PLUTO 09 vs Python comparison -> $(OUTDIR)/compare-09/"
	@mkdir -p $(OUTDIR)/compare-09
	@$(PYTHON) -m compare.overlay --setup 09 --out $(OUTDIR)/compare-09
	@$(PYTHON) -m compare.norms   --setup 09 --out $(OUTDIR)/compare-09

compare-02+09: $(OUTDIR)/pluto/02+09/done $(PY_SENTINEL)
	@echo ">>> PLUTO 02+09 vs Python comparison -> $(OUTDIR)/compare-02+09/"
	@mkdir -p $(OUTDIR)/compare-02+09
	@$(PYTHON) -m compare.overlay --setup 02+09 --out $(OUTDIR)/compare-02+09
	@$(PYTHON) -m compare.norms   --setup 02+09 --out $(OUTDIR)/compare-02+09

compare-python-all: compare compare-02 compare-09 compare-02+09

compare-ct-01: $(OUTDIR)/pluto/01/done $(PY_CT_SENTINEL)
	@echo ">>> PLUTO 01 vs CT Python comparison -> $(OUTDIR)/compare-ct-01/"
	@mkdir -p $(OUTDIR)/compare-ct-01
	@$(PYTHON) -m compare.overlay --setup 01 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-01
	@$(PYTHON) -m compare.norms   --setup 01 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-01

compare-ct-02: $(OUTDIR)/pluto/02/done $(PY_CT_SENTINEL)
	@echo ">>> PLUTO 02 vs CT Python comparison -> $(OUTDIR)/compare-ct-02/"
	@mkdir -p $(OUTDIR)/compare-ct-02
	@$(PYTHON) -m compare.overlay --setup 02 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-02
	@$(PYTHON) -m compare.norms   --setup 02 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-02

compare-ct-09: $(OUTDIR)/pluto/09/done $(PY_CT_SENTINEL)
	@echo ">>> PLUTO 09 vs CT Python comparison -> $(OUTDIR)/compare-ct-09/"
	@mkdir -p $(OUTDIR)/compare-ct-09
	@$(PYTHON) -m compare.overlay --setup 09 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-09
	@$(PYTHON) -m compare.norms   --setup 09 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-09

compare-ct-02+09: $(OUTDIR)/pluto/02+09/done $(PY_CT_SENTINEL)
	@echo ">>> PLUTO 02+09 vs CT Python comparison -> $(OUTDIR)/compare-ct-02+09/"
	@mkdir -p $(OUTDIR)/compare-ct-02+09
	@$(PYTHON) -m compare.overlay --setup 02+09 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-02+09
	@$(PYTHON) -m compare.norms   --setup 02+09 --python-in $(OUTDIR)/python-ct/data.npz --out $(OUTDIR)/compare-ct-02+09

compare-ct-all: compare-ct-01 compare-ct-02 compare-ct-09 compare-ct-02+09

# ---- animations ----
.PHONY: animate-01 animate-02 animate-09 animate-02+09 animate-pluto animate-python animate-python-ct animate-all animations

PLUTO_GIFS := $(ANIMDIR)/pluto_01_rho.gif $(ANIMDIR)/pluto_02_rho.gif $(ANIMDIR)/pluto_09_rho.gif $(ANIMDIR)/pluto_02+09_rho.gif
PYTHON_GIF := $(ANIMDIR)/python_rho.gif
PYTHON_CT_GIF := $(ANIMDIR)/python_ct_rho.gif

animate-01: $(ANIMDIR)/pluto_01_rho.gif
animate-02: $(ANIMDIR)/pluto_02_rho.gif
animate-09: $(ANIMDIR)/pluto_09_rho.gif
animate-02+09: $(ANIMDIR)/pluto_02+09_rho.gif
animate-pluto: $(PLUTO_GIFS)
animate-python: $(PYTHON_GIF)
animate-python-ct: $(PYTHON_CT_GIF)
animate-all animations: animate-pluto animate-python animate-python-ct

$(ANIMDIR)/pluto_01_rho.gif: $(OUTDIR)/pluto/01/done scripts/animate.py
	@echo ">>> PLUTO setup 01 density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source pluto --setup 01 --var rho --out $@

$(ANIMDIR)/pluto_02_rho.gif: $(OUTDIR)/pluto/02/done scripts/animate.py
	@echo ">>> PLUTO setup 02 density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source pluto --setup 02 --var rho --out $@

$(ANIMDIR)/pluto_09_rho.gif: $(OUTDIR)/pluto/09/done scripts/animate.py
	@echo ">>> PLUTO setup 09 density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source pluto --setup 09 --var rho --out $@

$(ANIMDIR)/pluto_02+09_rho.gif: $(OUTDIR)/pluto/02+09/done scripts/animate.py
	@echo ">>> PLUTO setup 02+09 density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source pluto --setup 02+09 --var rho --out $@

$(PYTHON_GIF): $(PY_SENTINEL) scripts/animate.py
	@echo ">>> Python solver density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source python --var rho --out $@

$(PYTHON_CT_GIF): $(PY_CT_SENTINEL) scripts/animate.py
	@echo ">>> CT Python solver density GIF -> $@"
	@mkdir -p $(ANIMDIR)
	@$(PYTHON) scripts/animate.py --source python --python-in $(OUTDIR)/python-ct/data.npz --var rho --title "Python CT solver" --out $@

# ---- all ----
.PHONY: all
all: plot-all plot-python plot-python-ct compare-schemes compare-python-all compare-python-solvers compare-ct-all animate-all
	@echo ">>> All done."

# ---- clean ----
.PHONY: clean
clean:
	@rm -rf $(WORKDIR) $(OUTDIR) $(PLOTDIR) $(MPLCONFIGDIR)
	@echo ">>> Cleaned."
