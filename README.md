# Axion SPK Analyzer

Browser analysis of Axion BioSystems spike files (`.spk`, written by AxIS / AxIS Navigator on Maestro systems), with Python twins of every analysis.

The parser is a re-implementation of Axion's official MATLAB **AxionFileLoader**
(`AxisFile.m`, `CombinedBlockVectorHeaderEntry.m`, `DiscontinuousBlockVectorHeaderEntry.m`,
`SpikeDataSet.m`, `ChannelArray.m`, the tag classes…). The SAMS repository ships a copy
under `src/AxionFileLoader`. All header CRC32 checksums are verified, as in the original loader.

## `mea_pipeline.html` (new, in progress)

This is a new page for our own analysis pipeline, independent of SAMS, built one panel at a time. Its styles are in `css/style.css`, which all the pages share. Its code is in three files in `js/`, loaded in this order:
- `js/spk_reader.js`: the `.spk` reader (`parseSpk`);
- `js/analysis.js`: every function that computes a result, with an index at the top. These are the recording preparation, activity counts, waveform measures, firing statistics, network bursts and comparison values. They take the recording and the parameters as arguments and never read the page or draw. The verification scripts run this file in Node;
- `js/mea_pipeline.js`: the page itself (state, drawing, controls, tabs).

They are plain scripts, not ES modules, so the page still works when opened straight from disk. Its `.spk` reader is the verified one (see [verification/](verification/README.md)): the page's parser is bit-identical to the Python reader on 7 recordings.

Tabs:
- **Open .spk…** accepts several files at once (drag and drop too). Each recording opens in its own tab, labelled with the file name, analysed with the default parameters; each tab keeps its own parameters, view and selection. The header stays the same for all tabs.
- **Multi-recording analysis** sets selected recordings side by side (up to 8). **Add .spk files…** there analyses files that aren't open yet with the default parameters, each in a new tab; a file already open is not read again. Each recording uses the parameters of its own tab, and any that differ from the defaults are listed. When a selected recording was analysed with parameters other than the defaults, a warning says so (and whether the recordings differ from each other), and offers to re-analyse all open recordings with the default Firing statistics and Network bursts parameters, after a confirmation. It compares network bursts per method, for one well or all wells: count, rate (60 · count / period, per min) and median duration (end − start; all wells' bursts pooled). These appear in a spider plot (one axis per ticked measure and method, up to 9, and one outline per recording; each axis scaled to its own maximum or shared within a measure) and in a table. A second spider plot and table compare the firing of each recording: network firing rate (all spikes of the recorded electrodes / period, spikes/s), mean electrode firing rate (mean ± SD over the active electrodes, Hz) and the number of active electrodes (`firingCompare` in `js/analysis.js`).

Current panels in each recording's tab:
- summary tiles;
- an activity replay, which drives the plate map:
  - controls: play/pause (Space), stop, speeds from 0.1× to 60×, and ← → to step one activity window (Shift: ten);
  - the timeline shows spikes/s over the whole recording; click or drag it to move the playhead;
  - at the playhead, each electrode shows the number of spikes it fired in the preceding activity window (10 ms – 1 s, default 100 ms). This is a plain count in the causal window (t − Δ, t], with no smoothing kernel.
  - the colour scale is fixed for the whole replay, saturating at the 99th percentile of non-zero counts;
  - the raster shows the playhead, and its window follows it;
- the plate map and the spike raster side by side in one panel (collapsed and expanded together; the plate takes only the width it needs, the raster the rest; stacked on windows narrower than 1000 px):
- a plate map coloured by firing rate (linear or log scale), spike count or median peak-to-peak amplitude, showing AxIS well treatments, with hover details and electrode selection by click or arrow keys;
- a spike raster with a population-rate strip:
  - zoom by dragging, ◀ ▶ to move, "Zoom out", "Whole recording", or typed start and end times (seconds or m:ss);
  - filter by well, and show electrodes with spikes or all recorded electrodes;
  - hover for per-electrode and per-bin details; click a row to select its electrode, linked with the plate map;
- waveforms per electrode, for the spikes in the raster window:
  - samples are drawn exactly as recorded, as straight lines between samples with no interpolation or smoothing;
  - alignment on the stored threshold crossing, or on the negative peak, positive peak or largest absolute sample. Alignment only shifts by whole samples.
  - the selected electrode is shown as individual spikes or as a density of all spikes, normalised per time column;
  - mean ± SD, with dots on the recorded samples, and the detector's ± threshold;
  - measures of the mean waveform: trough, peak, peak-to-peak, trough-to-next-peak time, noise SD and signal-to-noise;
  - a grid of every electrode of the well, laid out as on the plate, with a shared or per-electrode amplitude scale; click to select;
- firing statistics, over the whole recording or the raster window:
  - per electrode: mean rate, active (≥ 5 spikes/min), Fano factor (1 s bins), mean and median ISI, CV, CV2, LV, and the % of ISIs below a limit;
  - charts of the selected electrode: rate over time, log-ISI histogram and autocorrelogram;
  - per-well summaries (median, or mean ± SD over active electrodes), grouped by treatment;
  - CSV export, and plate-map colours by statistic;
  - an **i** button next to every statistic opens its definition, with equations and references;
  - a "detector limit" note: the shortest interval in the file, below which the ISI histogram and autocorrelogram are empty because of the detector's dead time, not the neurons;
  - `python/firing_stats.py` computes the same statistics in Python, and `verification/compare_firing_stats.py` checks they match;
- network bursts per well, on the pooled train of the well's active electrodes, over the period of Firing statistics, with three methods side by side:
  1. ISI<sub>N</sub> (Bakkum et al. 2014): N spikes within θ<sub>N</sub>. The defaults are N = 10 and θ<sub>N</sub> = 25 ms, with ≥ 10 spikes on ≥ 3 electrodes (MEA-NAP's rule). **Valley** sets θ<sub>N</sub> from the ISI<sub>N</sub> histogram;
  2. burst grouping (after Chiappalone 2006 / Pasquale 2010): single-electrode bursts (ISI ≤ 100 ms, ≥ 5 spikes) whose onsets are ≤ 100 ms apart, on ≥ 20 % of active electrodes;
  3. recruitment: runs of 50 ms bins in which ≥ 25 % of active electrodes fire. This is an active-electrode-count detector of the kind reviewed by Välkki et al. 2017 (Front Comput Neurosci 11:40); the exact rule and the 50 ms / 25 % defaults are this pipeline's own, with no published reference;
  - every parameter is editable, and every method requires ≥ 3 electrodes;
  - a timeline shows electrodes firing per bin, the raster of the well and one lane per method. Navigate with an overview strip (drag), typed times, or by clicking a network burst to zoom on it;
  - ISI<sub>N</sub> histogram, agreement matrix between the methods, a per-well/per-method statistics table, and CSV export (all network bursts, summary);
  - `python/network_bursts.py` is the Python twin, and `verification/compare_network_bursts.py` checks they match;
- recording details.

## HTML only: no installation

`mea_pipeline.html` needs the shared stylesheet `css/style.css` and the `js` folder next to it. It has no external scripts, libraries, fonts or server, and it works offline.
Double-click it to open it in any modern browser (Chrome, Edge, Firefox, Safari), then drag `.spk` files onto the page or use **Open .spk…**.
Files are parsed locally in the browser and nothing is uploaded. To copy the app elsewhere or send it to colleagues, include the `css` and `js` folders with it.

The former `spk_viewer.html` (file viewer with SAMS-style spike sorting) was removed on 2026-11-23. Its PCA and spike sorting, the next analysis step of the pipeline, are documented in [docs/pca_spike_sorting.md](docs/pca_spike_sorting.md), with the code kept for reference in `docs/pca_sorting_reference.js`.

## Optional Python tools

They are in the `python/` folder. To import them from a script in another folder, add `python/` to `PYTHONPATH` (or to `sys.path`).

| File | What it does |
|---|---|
| `python/axion_spk.py` | Python/numpy reader for scripted analysis, plus a CLI (summary, CSV and NPZ export). |
| `python/axion_sort.py` | SAMS-style spike sorting (PCA, clustering, template merge) in Python/numpy, for scripts and batch runs. It is the twin of the former viewer's sorter; see [docs/pca_spike_sorting.md](docs/pca_spike_sorting.md). |
| `python/firing_stats.py` | Python twin of the pipeline's firing statistics (checked by `verification/compare_firing_stats.py`). |
| `python/network_bursts.py` | Python twin of the pipeline's network bursts (checked by `verification/compare_network_bursts.py`). |

```bash
python3 python/axion_spk.py path/to/file.spk --csv spikes.csv --npz spikes.npz
python3 python/axion_sort.py path/to/file.spk --units units.csv --csv spikes_with_units.csv
python3 python/axion_sort.py path/to/file.spk --electrodes A6_11 B4_12 --merge-thr 1.2   # subset, custom parameters
```

```python
from axion_spk import read_spk
f = read_spk("file.spk")
f.summary()
f.spike_times_s        # trigger time of each spike (s)
f.electrode_labels()   # e.g. 'C3_33' (well C3, electrode column 3, row 3)
f.waveforms_uv         # (n_spikes, n_samples) in µV
f.channels, f.header, f.metadata, f.tags

res = f.sort()               # or f.sort(electrodes=["A6_11"], merge_thr=1.2)
res.unit                     # unit of every spike: -1 not sorted, 0 unsorted, 1.. unit
res.electrodes["A6_11"].units[0].mean_uV   # template, plus n, rate_hz, snr, isi_cv, …
res.write_units_csv("units.csv")           # one row per unit
```

## Plates and electrode grids

The plate layout comes from each file's channel map, so any Axion plate works. Tested on:
- the 96-well SAMS sample: CytoView MEA 96, 3×3 electrodes per well;
- six recordings on CytoView MEA 6: 6 wells, 8×8 = 64 electrodes per well, 384 channels.

For dense grids the views adapt:
- the raster lists only electrodes with spikes, with a switch to show every recorded electrode;
- the raster grows taller with the number of rows, and when rows are thin it labels wells;
- the waveform grid widens each well block to fit 8×8 panels;

The plate is named as in AxIS (the "Plate Type" recording setting). `SpikeFile.plate_name` gives the same name in Python.

## Format notes (Spike_v1)

Each spike is one block of `BlockHeaderSize + 2*NumSamplesPerBlock` bytes, starting at `DataRegionStart`:
`int64 startingSample, uint8 channel, uint8 chip, int32 triggerSample, float64 noiseSD (V), float64 thresholdMultiplier, [reserved], int16 data[n]`.
Spike time = `(startingSample + triggerSample) / SamplingFrequency`, and voltage = `data * VoltageScale` (V).
The pair `(chip, channel)` is looked up in the ChannelArray to get the well and electrode position.

Only header version 1.x files are supported, which covers all modern AxIS files. Legacy v0.x files are not.
