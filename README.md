# VTAPnO_Locomotion_Analysis

Code and representative example datasets associated with the manuscript:

**A VTA-pontine GABA pathway biases locomotor direction via local and distal inhibition**

This repository contains analysis workflows used for:
- open-field locomotion direction classification from DeepLabCut tracking
- fiber photometry PSTH processing
- opto-tagging classification

## Repository structure

- `notebooks/`  
  Jupyter notebooks for behavior and fiber photometry analyses

- `matlab/`  
  MATLAB script for opto-tagging classification

- `examples/behavior/`  
  Representative example input files for the locomotion direction workflow

- `examples/photometry/`  
  Representative example input files for the fiber photometry workflow

- `examples/optotag/`  
  Representative example input files for the opto-tagging workflow

- `output_examples/`  
  Example outputs generated from the provided example data

## Included workflows

### 1. Open-field locomotion direction analysis
Notebook:
- `notebooks/Backward_Loc_min.ipynb`

This workflow uses DeepLabCut tracking data to classify movement direction as forward, backward, leftward, or rightward, compute peri-stimulus movement proportions, and fit a linear mixed-effects model.

Example input files are provided in:
- `examples/behavior/`

### 2. Fiber photometry PSTH analysis
Notebook:
- `notebooks/Fiber_Photometry_PSTHv2.ipynb`

This workflow processes Doric-exported signal files and event timestamps to compute ΔF/F, z-scored peri-stimulus traces, and summary PSTHs.

Example input files are provided in:
- `examples/photometry/`

Note: the example signal file is a temporally trimmed recording segment containing 3 representative stimulation events from a longer recording. The original sampling rate and signal structure were preserved.

### 3. Opto-tagging classification
MATLAB script:
- `matlab/optoTag_classification_unifiedCriteria_v2.m`

This workflow classifies optogenetically responsive units using spike-time data and stimulus timestamps, based on significance against a shuffled null distribution and a fidelity criterion.

Example input files are provided in:
- `examples/optotag/`

## Requirements

Install Python dependencies with:

```bash
pip install -r requirements.txt
