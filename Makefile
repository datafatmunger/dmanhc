PYTHON ?= .venv/bin/python3
SRC ?= src
BUILD_DIR := build
EXPERIMENTS ?= experiments

.DEFAULT_GOAL := dmanh

.PHONY: clean compile plots comparison dmanh dmanh-notebook dmanh-frequency-sweep readout-x-trace notebook-observable-audit non-abelian-flux

clean:
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)

# --- Symmetric double well (McGarry baseline) ---

compile: clean
	$(PYTHON) $(SRC)/compiler.py $(EXPERIMENTS)/double_well.toml

plots: compile
	$(PYTHON) $(SRC)/plots.py $(EXPERIMENTS)/double_well.toml
	$(PYTHON) $(SRC)/measure.py $(EXPERIMENTS)/double_well.toml

comparison: compile
	$(PYTHON) $(SRC)/measure.py $(EXPERIMENTS)/double_well.toml \
		--comparison-output build/double_well_measurement_comparison.png \
		--comparison-title 'Direct state access versus McGarry characteristic-function readout'

# --- DMANH+ (Phil's parameters) ---

dmanh:
	$(PYTHON) $(SRC)/compiler.py $(EXPERIMENTS)/dmanh.toml
	$(PYTHON) $(SRC)/plots.py $(EXPERIMENTS)/dmanh.toml
	$(PYTHON) $(SRC)/measure.py $(EXPERIMENTS)/dmanh.toml

dmanh-notebook:
	$(PYTHON) $(SRC)/compiler.py $(EXPERIMENTS)/dmanh.toml

# --- DMANH+ frequency sweep ---

dmanh-frequency-sweep:
	$(PYTHON) $(SRC)/compiler.py $(EXPERIMENTS)/dmanh_frequency_sweep.toml
	$(PYTHON) $(SRC)/frequency_sweep.py $(EXPERIMENTS)/dmanh_frequency_sweep.toml

readout-x-trace:
	$(PYTHON) $(SRC)/readout_x_trace.py

notebook-observable-audit:
	$(PYTHON) $(SRC)/notebook_observable_audit.py
