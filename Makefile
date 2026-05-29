PYTHON ?= .venv/bin/python3
SRC ?= src
BUILD_DIR := build
JAQAL ?= $(BUILD_DIR)/double_well.jaqal
DIRECT_PNG ?= $(BUILD_DIR)/double_well_oracle.png
HSIM_PNG ?= $(BUILD_DIR)/double_well_hsim.png
MEASUREMENT_PNG ?= $(BUILD_DIR)/double_well_measurement_panels.png
CHI_PNG ?= $(BUILD_DIR)/double_well_chi_slice_panels.png
COMPARISON_PNG ?= $(BUILD_DIR)/double_well_measurement_comparison.png
DIRECT_TITLE ?= Symmetric double well # from ideal truncated-Fock simulator
HSIM_TITLE ?= Exact $$H_{\mathrm{sim}}$$ versus compiled-gate dynamics
MEASUREMENT_TITLE ?= McGarry Eq. 33 readout from chi(beta)
CHI_TITLE ?= McGarry Fig. 7-style characteristic-function slice
COMPARISON_TITLE ?= Direct state access versus McGarry characteristic-function readout

DMANH_JAQAL ?= $(BUILD_DIR)/protocol_b_dmanh.jaqal
DMANH_DIRECT_PNG ?= $(BUILD_DIR)/protocol_b_dmanh.png
DMANH_DIRECT_TITLE ?= Protocol B DMANH # from ideal truncated-Fock simulator

.DEFAULT_GOAL := plots

.PHONY: clean compile plots comparison dmanh

clean:
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)

compile: clean
	$(PYTHON) $(SRC)/compiler.py --output $(JAQAL)

plots: compile
	$(PYTHON) $(SRC)/plots.py --jaqal $(JAQAL) --output $(DIRECT_PNG) --title '$(DIRECT_TITLE)' --hsim-output $(HSIM_PNG) --hsim-title '$(HSIM_TITLE)'
	$(PYTHON) $(SRC)/measure.py --jaqal $(JAQAL) --output $(MEASUREMENT_PNG) --title '$(MEASUREMENT_TITLE)' --chi-output $(CHI_PNG) --chi-title '$(CHI_TITLE)'

comparison: compile
	$(PYTHON) $(SRC)/measure.py --jaqal $(JAQAL) --output $(MEASUREMENT_PNG) --title '$(MEASUREMENT_TITLE)' --chi-output $(CHI_PNG) --chi-title '$(CHI_TITLE)' --comparison-output $(COMPARISON_PNG) --comparison-title '$(COMPARISON_TITLE)'

dmanh:
	$(PYTHON) $(SRC)/compiler.py \
		--output $(DMANH_JAQAL) \
		--max-time-ms 4 \
		--dt-us 200 \
		--delta-hz 754.95 \
		--alpha0 0.49365 \
		--vartheta 0.8 \
		--varphi 0 \
		--x-min 1.208 \
		--probe-qubit-index 0
	$(PYTHON) $(SRC)/plots.py \
		--jaqal $(DMANH_JAQAL) \
		--output $(DMANH_DIRECT_PNG) \
		--title '$(DMANH_DIRECT_TITLE)' \
		--no-hsim-output
