# Behavior example dataset

This folder contains representative example input data for the open-field locomotion direction analysis based on DeepLabCut tracking.

## Files
For each animal/session, the pipeline expects:
- `X.h5`: DeepLabCut tracking file
- `stim_X.csv`: stimulus frame file

where `X` is the animal/session index (for example, `1.h5` and `stim_1.csv`).

## Expected file structure

### `X.h5`
The H5 file must be a DeepLabCut output file containing tracked body-part coordinates.

This analysis uses:
- `Back`
- `Tailbase`

Coordinates are read from the DeepLabCut multi-index structure.

### `stim_X.csv`
The stimulus file must be a CSV containing stimulus frame numbers.

The code reads stimulus frames from the **first column**.

## Used by
- `notebooks/Backward_Loc_min.ipynb`

## What the analysis does
The pipeline:
- loads DeepLabCut tracking data
- smooths positional coordinates
- computes movement and orientation vectors
- classifies movement into:
  - Backward
  - Forward
  - Leftward
  - Rightward
- computes peri-stimulus movement proportions in 200 ms bins
- fits a linear mixed-effects model (LMM)
- generates summary plots and aggregated output

## Running with multiple animals

The analysis pipeline is designed to process multiple animals/sessions.

In the full workflow, files are named:
- `1.h5`, `2.h5`, ..., `N.h5`
- `stim_1.csv`, `stim_2.csv`, ..., `stim_N.csv`

In this repository, a reduced set of example files is provided to demonstrate the workflow.

To analyze a different number of animals, modify:

```python
for mouse_id in range(1, N+1):
