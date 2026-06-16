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

# ---- default ----
.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "Project 2 — MHD Shock Cloud"
	@echo ""
	@echo "Targets:"
	@echo "  make setup   — copy config files into build/ and generate PLUTO Makefile"
	@echo "  make pluto   — build and run the PLUTO simulation"
	@echo "  make plot    — pyPLUTO visualization -> plots/"
	@echo "  make python  — run the independent Python solver"
	@echo "  make compare — overlay + L1/L2 norms analysis"
	@echo "  make all     — setup -> pluto -> plot -> python -> compare"
	@echo "  make clean   — wipe build/ output/ plots/"
	@echo ""
	@echo "PLUTO_DIR = $(PU_DIR)"

# ---- setup working directory for PLUTO ----
.PHONY: setup
setup: $(WORKDIR)
	@echo ">>> Setting up PLUTO build in $(WORKDIR)"
	@cp inputs/pluto.ini $(WORKDIR)/
	@cp inputs/definitions.h $(WORKDIR)/
	@cp inputs/init.c $(WORKDIR)/
	@cd $(WORKDIR) && python3 $(PU_DIR)/setup.py --auto-update
	@echo ">>> Setup complete. Run 'make pluto' to build."

$(WORKDIR):
	@mkdir -p $(WORKDIR)

# ---- build and run PLUTO ----
.PHONY: pluto
pluto:
	@echo ">>> Building PLUTO for Shock_Cloud..."
	@$(MAKE) -C $(WORKDIR) PLUTO_DIR=$(PU_DIR)
	@echo ">>> Running PLUTO simulation..."
	@cd $(WORKDIR) && ./pluto -i pluto.ini
	@mkdir -p $(OUTDIR)/pluto
	@cp $(WORKDIR)/data.*.dbl $(WORKDIR)/dbl.out $(OUTDIR)/pluto/ 2>/dev/null || true
	@echo ">>> PLUTO output in $(OUTDIR)/pluto/"

# ---- pyPLUTO plots ----
.PHONY: plot
plot:
	@echo ">>> Generating pyPLUTO plots..."
	@PYTHONPATH="$(PU_DIR)/Tools/Python:$$PYTHONPATH" python3 scripts/plot_pluto.py
	@echo ">>> Plots in $(PLOTDIR)/"

# ---- Python reproduction ----
.PHONY: python
python:
	@echo ">>> Running Python MHD solver..."
	@python3 python/solver.py
	@echo ">>> Python output in $(OUTDIR)/python/"

# ---- comparison ----
.PHONY: compare
compare:
	@echo ">>> Running comparison analysis..."
	@python3 compare/overlay.py
	@python3 compare/norms.py
	@echo ">>> Comparison output in $(OUTDIR)/compare/"

# ---- all ----
.PHONY: all
all: setup pluto plot python compare
	@echo ">>> All done."

# ---- clean ----
.PHONY: clean
clean:
	@rm -rf $(WORKDIR) $(OUTDIR) $(PLOTDIR)
	@echo ">>> Cleaned."
