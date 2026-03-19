# Opto-tagging example dataset

This folder contains representative example input data for the MATLAB-based opto-tagging classification workflow used in this study.

## Files
- Curated_clusters_spike_times_columns_sorted1.csv: example spike-times file
- Stim_B1.mat (or Stim_R1.mat): example stimulus file

## Expected file structure

### Curated_clusters_spike_times_columns_sorted1.csv
This CSV contains spike times in samples.

- The first row contains neuron IDs
- Each column corresponds to one neuron
- Each column contains spike times (in samples)
- Empty entries or zero-padding are ignored by the script

### Stim_B1.mat / Stim_R1.mat
The MATLAB stimulus file must contain a variable named:

selectedPoints

This variable must contain stimulus timestamps in samples.

## Used by
- matlab/optoTag_classification_unifiedCriteria_v2.m

## What the analysis does
The script:
- loads spike times and stimulus timestamps
- converts spike and stimulus times from samples to milliseconds
- computes trial-aligned post-stimulus spike counts
- applies one of two opto-tagging criteria:
  - 2ms_peak_anywhere
  - 6ms_fixed
- builds a shuffled null distribution using circular trial-wise shifts
- computes a 99.9th percentile significance threshold
- applies a fidelity criterion (default: at least 10% of trials with a spike in the criterion window)
- outputs a summary CSV with classification results
- generates quick aggregate plots for passing units

## User settings
The main analysis parameters are defined at the top of the MATLAB script, including:
- tagging criterion
- input filenames
- sampling rate
- bin size
- shuffling window
- number of shuffles
- fidelity threshold

## Running with different files
The script is designed to run on one spike-times file and one stimulus file at a time.

To analyze a different recording, modify the following lines in the MATLAB script:

spikeCsv = 'Curated_clusters_1.csv';
stimMat  = 'Stim_B1.mat';

Ensure that:
- the spike CSV follows the expected format
- the MAT file contains the variable 'selectedPoints'

## Output
The script generates:
- optotag_unified_<criterion>.csv: summary results table
- aggregate PSTH plot for passing units
- heatmap of concatenated trial-by-time data for passing units

## Notes
This example dataset is provided to allow reviewers and users to run the opto-tagging workflow on representative data. The full datasets used in the study are available from the corresponding author upon reasonable request.
