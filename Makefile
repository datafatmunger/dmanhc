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
ALPHA_PHASE_OFFSET ?= 0 # -1.5707963267948966

DOUBLE_WELL_STEPS ?= 20
DOUBLE_WELL_B_RAD_S ?= 4.0e3
DOUBLE_WELL_DELTA_RAD_S ?= 3.141592653589793e3
DOUBLE_WELL_ALPHA0 ?= 0.5235987755982989
DOUBLE_WELL_X_MIN ?= 1.5
DOUBLE_WELL_TIMES_MS ?= 0 2.00 4.00

DMANH_JAQAL ?= $(BUILD_DIR)/dmanh.jaqal
DMANH_DIRECT_PNG ?= $(BUILD_DIR)/dmanh.png
DMANH_HSIM_PNG ?= $(BUILD_DIR)/dmanh_hsim.png
DMANH_MEASUREMENT_PNG ?= $(BUILD_DIR)/dmanh_measurement_panels.png
DMANH_CHI_PNG ?= $(BUILD_DIR)/dmanh_chi_slice_panels.png
DMANH_DIRECT_TITLE ?= DMANH+
DMANH_HSIM_TITLE ?= DMANH+ exact $$H_{\mathrm{sim}}$$ versus compiled-gate dynamics
DMANH_MEASUREMENT_TITLE ?= DMANH+ McGarry Eq. 33 readout from chi(beta)
DMANH_CHI_TITLE ?= DMANH+ characteristic-function slice
DMANH_STEPS ?= 49
DMANH_B_RAD_S ?= 5.09628e3
DMANH_DELTA_RAD_S ?= 1.29817e3
DMANH_ALPHA0 ?= 0.18512
DMANH_X_MIN ?= 1.25895
DMANH_TIMES_MS ?= 0 4.081408 7.691885
DMANH_HSIM_MAX_TIME_MS ?= 7.6918850612603702

DMANH_VARTHETA_1P6_JAQAL ?= $(BUILD_DIR)/dmanh_vartheta_1p6.jaqal
DMANH_VARTHETA_1P6_DIRECT_PNG ?= $(BUILD_DIR)/dmanh_vartheta_1p6.png
DMANH_VARTHETA_1P6_HSIM_PNG ?= $(BUILD_DIR)/dmanh_vartheta_1p6_hsim.png
DMANH_VARTHETA_1P6_MEASUREMENT_PNG ?= $(BUILD_DIR)/dmanh_vartheta_1p6_measurement_panels.png
DMANH_VARTHETA_1P6_CHI_PNG ?= $(BUILD_DIR)/dmanh_vartheta_1p6_chi_slice_panels.png
DMANH_VARTHETA_1P6_DIRECT_TITLE ?= DMANH+ vartheta=1.6
DMANH_VARTHETA_1P6_HSIM_TITLE ?= DMANH+ vartheta=1.6 exact $$H_{\mathrm{sim}}$$ versus compiled-gate dynamics
DMANH_VARTHETA_1P6_MEASUREMENT_TITLE ?= DMANH+ vartheta=1.6 McGarry Eq. 33 readout from chi(beta)
DMANH_VARTHETA_1P6_CHI_TITLE ?= DMANH+ vartheta=1.6 characteristic-function slice
DMANH_VARTHETA_1P6 ?= 1.6
DMANH_VARTHETA_1P6_STEPS ?= 25
DMANH_VARTHETA_1P6_TIMES_MS ?= 0 4.081408399852441 7.848862307408541
DMANH_VARTHETA_1P6_HSIM_MAX_TIME_MS ?= 7.848862307408541

DMANH_FREQUENCY_CSV ?= $(BUILD_DIR)/dmanh_frequency_sweep.csv
DMANH_FREQUENCY_SHIFT_PNG ?= $(BUILD_DIR)/dmanh_frequency_shift_vs_vartheta.png
DMANH_FREQUENCY_SPECTRA_PNG ?= $(BUILD_DIR)/dmanh_frequency_spectra.png
DMANH_FREQUENCY_TRACE_PNG ?= $(BUILD_DIR)/dmanh_frequency_trace_examples.png
DMANH_FREQUENCY_MAX_TIME_MS ?= 160
DMANH_FREQUENCY_VARTHETA_MIN ?= 0.1
DMANH_FREQUENCY_VARTHETA_MAX ?= 3.0
DMANH_FREQUENCY_VARTHETA_STEP ?= 0.1
DMANH_FREQUENCY_SELECTED_VARTHETA ?= 0.1 0.8 1.6 3.0
DMANH_FREQUENCY_TITLE ?= DMANH+ FFT peak shift versus compiled-gate timestep

.DEFAULT_GOAL := dmanh

DMANH_NOTEBOOK_DIR ?= $(BUILD_DIR)/notebook
DMANH_READOUT_BETAS ?= 0 -0.4 0 -0.2 0 0 0 0.2 0 0.4

.PHONY: clean compile plots comparison phil dmanh dmanh-notebook dmanh-vartheta-1p6 dmanh-frequency-sweep

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
		--hsim-output $(DMANH_HSIM_PNG) \
		--hsim-title '$(DMANH_HSIM_TITLE)' \
		--hsim-max-time-ms $(DMANH_HSIM_MAX_TIME_MS)
	$(PYTHON) $(SRC)/measure.py \
		--jaqal $(DMANH_JAQAL) \
		--output $(DMANH_MEASUREMENT_PNG) \
		--title '$(DMANH_MEASUREMENT_TITLE)' \
		--times-ms $(DMANH_TIMES_MS) \
		--chi-output $(DMANH_CHI_PNG) \
		--chi-title '$(DMANH_CHI_TITLE)'

dmanh-notebook:
	$(PYTHON) $(SRC)/compiler.py \
		--output $(DMANH_JAQAL) \
		--steps $(DMANH_STEPS) \
		--B-rad-s $(DMANH_B_RAD_S) \
		--delta-rad-s $(DMANH_DELTA_RAD_S) \
		--alpha0 $(DMANH_ALPHA0) \
		--vartheta $(VARTHETA) \
		--x-min $(DMANH_X_MIN) \
		--alpha-phase-offset $(ALPHA_PHASE_OFFSET) \
		--export-numpy $(DMANH_NOTEBOOK_DIR) \
		--readout-betas $(DMANH_READOUT_BETAS)

dmanh-vartheta-1p6:
	$(PYTHON) $(SRC)/compiler.py \
		--output $(DMANH_VARTHETA_1P6_JAQAL) \
		--steps $(DMANH_VARTHETA_1P6_STEPS) \
		--B-rad-s $(DMANH_B_RAD_S) \
		--delta-rad-s $(DMANH_DELTA_RAD_S) \
		--alpha0 $(DMANH_ALPHA0) \
		--vartheta $(DMANH_VARTHETA_1P6) \
		--x-min $(DMANH_X_MIN) \
		--alpha-phase-offset $(ALPHA_PHASE_OFFSET)
	$(PYTHON) $(SRC)/plots.py \
		--jaqal $(DMANH_VARTHETA_1P6_JAQAL) \
		--output $(DMANH_VARTHETA_1P6_DIRECT_PNG) \
		--title '$(DMANH_VARTHETA_1P6_DIRECT_TITLE)' \
		--times-ms $(DMANH_VARTHETA_1P6_TIMES_MS) \
		--hsim-output $(DMANH_VARTHETA_1P6_HSIM_PNG) \
		--hsim-title '$(DMANH_VARTHETA_1P6_HSIM_TITLE)' \
		--hsim-max-time-ms $(DMANH_VARTHETA_1P6_HSIM_MAX_TIME_MS)
	$(PYTHON) $(SRC)/measure.py \
		--jaqal $(DMANH_VARTHETA_1P6_JAQAL) \
		--output $(DMANH_VARTHETA_1P6_MEASUREMENT_PNG) \
		--title '$(DMANH_VARTHETA_1P6_MEASUREMENT_TITLE)' \
		--times-ms $(DMANH_VARTHETA_1P6_TIMES_MS) \
		--chi-output $(DMANH_VARTHETA_1P6_CHI_PNG) \
		--chi-title '$(DMANH_VARTHETA_1P6_CHI_TITLE)'

dmanh-frequency-sweep:
	$(PYTHON) $(SRC)/compiler.py \
		--output $(DMANH_JAQAL) \
		--steps $(DMANH_STEPS) \
		--B-rad-s $(DMANH_B_RAD_S) \
		--delta-rad-s $(DMANH_DELTA_RAD_S) \
		--alpha0 $(DMANH_ALPHA0) \
		--vartheta $(VARTHETA) \
		--x-min $(DMANH_X_MIN) \
		--alpha-phase-offset $(ALPHA_PHASE_OFFSET)
	$(PYTHON) $(SRC)/frequency_sweep.py \
		--jaqal $(DMANH_JAQAL) \
		--output-csv $(DMANH_FREQUENCY_CSV) \
		--shift-output $(DMANH_FREQUENCY_SHIFT_PNG) \
		--spectra-output $(DMANH_FREQUENCY_SPECTRA_PNG) \
		--trace-output $(DMANH_FREQUENCY_TRACE_PNG) \
		--max-time-ms $(DMANH_FREQUENCY_MAX_TIME_MS) \
		--vartheta-min $(DMANH_FREQUENCY_VARTHETA_MIN) \
		--vartheta-max $(DMANH_FREQUENCY_VARTHETA_MAX) \
		--vartheta-step $(DMANH_FREQUENCY_VARTHETA_STEP) \
		--selected-vartheta $(DMANH_FREQUENCY_SELECTED_VARTHETA) \
		--title '$(DMANH_FREQUENCY_TITLE)'
