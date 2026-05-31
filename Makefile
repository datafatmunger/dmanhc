PYTHON ?= .venv/bin/python3
SRC ?= src
BUILD_DIR := build
JAQAL ?= $(BUILD_DIR)/double_well.jaqal
DIRECT_PNG ?= $(BUILD_DIR)/double_well_oracle.png
HSIM_PNG ?= $(BUILD_DIR)/double_well_hsim.png
MEASUREMENT_PNG ?= $(BUILD_DIR)/double_well_measurement_panels.png
CHI_PNG ?= $(BUILD_DIR)/double_well_chi_slice_panels.png
COMPARISON_PNG ?= $(BUILD_DIR)/double_well_measurement_comparison.png
DIRECT_TITLE ?= Symmetric double well
HSIM_TITLE ?= Exact $$H_{\mathrm{sim}}$$ versus compiled-gate dynamics
MEASUREMENT_TITLE ?= McGarry Eq. 33 readout from chi(beta)
CHI_TITLE ?= McGarry Fig. 7-style characteristic-function slice
COMPARISON_TITLE ?= Direct state access versus McGarry characteristic-function readout

VARTHETA ?= 0.8
ALPHA_PHASE_OFFSET ?= -1.5707963267948966

DOUBLE_WELL_STEPS ?= 20
DOUBLE_WELL_B_RAD_S ?= 4.0e3
DOUBLE_WELL_DELTA_RAD_S ?= 3.141592653589793e3
DOUBLE_WELL_ALPHA0 ?= 0.5235987755982989
DOUBLE_WELL_X_MIN ?= 1.5
DOUBLE_WELL_TIMES_MS ?= 0 2.00 4.00

DMANH_JAQAL ?= $(BUILD_DIR)/dmanh.jaqal
DMANH_DIRECT_PNG ?= $(BUILD_DIR)/dmanh.png
DMANH_DIRECT_TITLE ?= Phil DMANH
DMANH_STEPS ?= 49
DMANH_B_RAD_S ?= 5.09628e3
DMANH_DELTA_RAD_S ?= 1.29817e3
DMANH_ALPHA0 ?= 0.18512
DMANH_X_MIN ?= 1.208
DMANH_TIMES_MS ?= 0 4.081408 7.691885

.DEFAULT_GOAL := dmanh

.PHONY: clean compile plots comparison phil dmanh

clean:
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)

compile: clean
	$(PYTHON) $(SRC)/compiler.py \
		--output $(JAQAL) \
		--steps $(DOUBLE_WELL_STEPS) \
		--B-rad-s $(DOUBLE_WELL_B_RAD_S) \
		--delta-rad-s $(DOUBLE_WELL_DELTA_RAD_S) \
		--alpha0 $(DOUBLE_WELL_ALPHA0) \
		--vartheta $(VARTHETA) \
		--x-min $(DOUBLE_WELL_X_MIN) \
		--alpha-phase-offset $(ALPHA_PHASE_OFFSET)

plots: compile
	$(PYTHON) $(SRC)/plots.py --jaqal $(JAQAL) --output $(DIRECT_PNG) --title '$(DIRECT_TITLE)' --times-ms $(DOUBLE_WELL_TIMES_MS) --hsim-output $(HSIM_PNG) --hsim-title '$(HSIM_TITLE)'
	$(PYTHON) $(SRC)/measure.py --jaqal $(JAQAL) --output $(MEASUREMENT_PNG) --title '$(MEASUREMENT_TITLE)' --times-ms $(DOUBLE_WELL_TIMES_MS) --chi-output $(CHI_PNG) --chi-title '$(CHI_TITLE)'

comparison: compile
	$(PYTHON) $(SRC)/measure.py --jaqal $(JAQAL) --output $(MEASUREMENT_PNG) --title '$(MEASUREMENT_TITLE)' --times-ms $(DOUBLE_WELL_TIMES_MS) --chi-output $(CHI_PNG) --chi-title '$(CHI_TITLE)' --comparison-output $(COMPARISON_PNG) --comparison-title '$(COMPARISON_TITLE)'

dmanh:
	$(PYTHON) $(SRC)/compiler.py \
		--output $(DMANH_JAQAL) \
		--steps $(DMANH_STEPS) \
		--B-rad-s $(DMANH_B_RAD_S) \
		--delta-rad-s $(DMANH_DELTA_RAD_S) \
		--alpha0 $(DMANH_ALPHA0) \
		--vartheta $(VARTHETA) \
		--x-min $(DMANH_X_MIN) \
		--alpha-phase-offset $(ALPHA_PHASE_OFFSET)
	$(PYTHON) $(SRC)/plots.py \
		--jaqal $(DMANH_JAQAL) \
		--output $(DMANH_DIRECT_PNG) \
		--title '$(DMANH_DIRECT_TITLE)' \
		--times-ms $(DMANH_TIMES_MS) \
		--no-hsim-output
