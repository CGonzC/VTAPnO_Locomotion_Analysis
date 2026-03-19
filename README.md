# VTAPnO_Locomotion_Analysis

Code and representative example datasets associated with the manuscript:

**A VTA-pontine GABA pathway biases locomotor direction via local and distal inhibition**

## Repository structure

- `notebooks/` – Python analysis workflows (behavior, photometry)  
- `matlab/` – Opto-tagging analysis script  
- `examples/` – Example input datasets for each workflow  
- `output_examples/` – Example outputs generated from the provided data  

## Workflows

### Behavior (DeepLabCut)
`notebooks/Backward_Loc_min.ipynb`  
Classifies locomotion direction and computes peri-stimulus movement metrics.

### Fiber photometry
`notebooks/Fiber_Photometry_PSTHv2.ipynb`  
Computes ΔF/F and peri-stimulus responses from photometry recordings.

### Opto-tagging (MATLAB)
`matlab/optoTag_classification_unifiedCriteria_v2.m`  
Identifies optogenetically responsive units using spike timing and shuffle-based significance.

## Requirements

Install dependencies with:

pip install -r requirements.txt

## Usage

Download or clone the repository and run:

- Python notebooks in `notebooks/`
- MATLAB script in `matlab/`

using the example files in `examples/`.

## Data availability

Example datasets are provided to demonstrate the analysis workflows.  
Full datasets used in the study are available from the corresponding author upon reasonable request.

## Code availability

All custom analysis code is provided in this repository.

## License

This code is released under the MIT License.
