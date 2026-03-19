# Fiber photometry example dataset

This folder contains representative example input data for the fiber photometry PSTH workflow used in this study.

## Files
- `1_signal.csv`: example Doric signal file containing fluorescence traces
- `1_events.csv`: example event file used for peri-stimulus alignment

## Expected file structure

### `1_signal.csv`
The signal file is a CSV exported from Doric and must contain at least the first three columns in this order:
1. time
2. 470-nm signal
3. 405-nm isosbestic signal

The analysis notebook reads only these first three columns.

### `1_events.csv`
The event file is a CSV containing event timestamps.
The notebook reads timestamps from the second column.

## Used by
- `notebooks/Fiber_Photometry_PSTHv2.ipynb`

## What the notebook does
The notebook:
- computes ΔF/F from the 470-nm and 405-nm signals using a robust 405→470 fit
- aligns trials to event times
- computes peri-stimulus averages
- generates ΔF/F, standard z-score, and robust z-score outputs

## Running with multiple recordings

The analysis pipeline is designed to process multiple recordings (paired signal and event files).

In this repository, a single example file pair (`1_signal.csv`, `1_events.csv`) is provided for demonstration.

To analyze multiple recordings, modify the file indexing in the notebook:

```python
signal_files = [DATA_DIR / f"{i}_signal.csv" for i in range(1, N+1)]
event_files  = [DATA_DIR / f"{i}_events.csv" for i in range(1, N+1)]
dff_files    = [DATA_DIR / f"dff_{i}.csv" for i in range(1, N+1)]
