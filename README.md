# Axion SPK Analyzer

A viewer for Axion BioSystems spike files (`.spk`, written by AxIS / AxIS Navigator on Maestro systems).

The parser is a re-implementation of Axion's official MATLAB **AxionFileLoader**
(`AxisFile.m`, `CombinedBlockVectorHeaderEntry.m`, `DiscontinuousBlockVectorHeaderEntry.m`,
`SpikeDataSet.m`, `ChannelArray.m`, the tag classes…). The SAMS repository ships a copy
under `src/AxionFileLoader`. All header CRC32 checksums are verified, as in the original loader.

## `mea_pipeline.html` (new, in progress)

This is a new page for our own analysis pipeline, independent of SAMS, built one panel at a time. Its styles are in `css/style.css`, which all the pages share. Its `.spk` reader is the verified one (see [verification/](verification/README.md)): the page's parser is bit-identical to the Python reader on 7 recordings.

Current panels:
- summary tiles;
- an activity replay, which drives the plate map:
  - controls: play/pause (Space), stop, speeds from 0.1× to 60×, and ← → to step one activity window (Shift: ten);
  - the timeline shows spikes/s over the whole recording; click or drag it to move the playhead;
  - at the playhead, each electrode shows the number of spikes it fired in the preceding activity window (10 ms – 1 s, default 100 ms). This is a plain count in the causal window (t − Δ, t], with no smoothing kernel.
  - the colour scale is fixed for the whole replay, saturating at the 99th percentile of non-zero counts;
  - the raster shows the playhead, and its window follows it;
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
  3. recruitment: runs of 50 ms bins in which ≥ 25 % of active electrodes fire;
  - every parameter is editable, and every method requires ≥ 3 electrodes;
  - a timeline shows electrodes firing per bin, the raster of the well and one lane per method. Navigate with an overview strip (drag), typed times, or by clicking a network burst to zoom on it;
  - ISI<sub>N</sub> histogram, agreement matrix between the methods, a per-well/per-method statistics table, and CSV export (all network bursts, summary);
  - `python/network_bursts.py` is the Python twin, and `verification/compare_network_bursts.py` checks they match;
- recording details.

## HTML only: no installation

`spk_viewer.html` needs only the shared stylesheet `css/style.css` next to it. It has no external scripts, libraries, fonts or server, and it works offline.
Double-click it to open it in any modern browser (Chrome, Edge, Firefox, Safari), then drag a `.spk` file onto the page or use **Open .spk…**.
The file is parsed locally in the browser and nothing is uploaded. To copy the viewer elsewhere or send it to colleagues, include the `css` folder with it.

### Loading from a URL

Paste a web address into the box at the top and press **Load URL**, or drop a link onto the page. The download shows a progress readout, and once the file loads the page address becomes a shareable link:
`spk_viewer.html?src=https://…/recording.spk`. Opening that link loads the file directly.

- The server must allow cross-origin downloads (CORS). GitHub raw, Zenodo, Figshare, OSF and CORS-enabled S3 / Google Cloud buckets all do.
  Many other sites don't. The viewer then explains why it couldn't load the file, and you can download the file and open it locally instead.
- GitHub `…/blob/…` page links and Dropbox share links are converted to their direct-download form automatically.
- If loading fails (404, CORS, or a file that isn't a .spk), the error is shown and the file that's already open stays on screen.

Everything the Python tools do is also available in the page. That includes the full spike-table export
(**Export all spikes (CSV)**), which has the same columns and values as `axion_spk.py --csv`.

## Optional Python tools

They are in the `python/` folder. To import them from a script in another folder, add `python/` to `PYTHONPATH` (or to `sys.path`).

| File | What it does |
|---|---|
| `python/axion_spk.py` | Python/numpy reader for scripted analysis, plus a CLI (summary, CSV and NPZ export). |
| `python/axion_sort.py` | The viewer's spike sorter in Python/numpy, for scripts and batch runs. It gives the same units, spike for spike, as the page. |
| `python/view_spk.py` | Opens a given `.spk` file straight in the viewer through a small local web server, with no drag and drop. |

```bash
python3 python/axion_spk.py path/to/file.spk --csv spikes.csv --npz spikes.npz
python3 python/view_spk.py path/to/file.spk
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
res.write_units_csv("units.csv")           # same file as the viewer's "Export units (CSV)"
```

## Plates and electrode grids

The plate layout comes from each file's channel map, so any Axion plate works. Tested on:
- the 96-well SAMS sample: CytoView MEA 96, 3×3 electrodes per well;
- six recordings on CytoView MEA 6: 6 wells, 8×8 = 64 electrodes per well, 384 channels.

For dense grids the views adapt:
- the raster lists only electrodes with spikes, with a switch to show every recorded electrode;
- the raster grows taller with the number of rows, and when rows are thin it labels wells;
- the waveform grid widens each well block to fit 8×8 panels;
- the electrode table hides electrodes without spikes, with a checkbox to show them.

The plate is named as in AxIS (the "Plate Type" recording setting). `SpikeFile.plate_name` gives the same name in Python.

## What the viewer shows

- **Summary tiles:** spike count, active and recorded electrodes, duration, sampling rate, plate type.
- **Activity replay:** play/pause (Space), stop, speeds from 0.1× to 60×, and ←/→ to step (Shift for bigger steps). The bar shows population spike rate over the whole recording; click or drag on it to seek. During replay the plate map shows each electrode's instantaneous firing rate, smoothed with an exponential kernel whose τ grows with speed (100 ms at 1×). A playhead also runs across the raster, and a zoomed raster window pages along with it.
- **Signal trace:** a full-width trace of the selected electrode, or of every electrode of its well stacked in lanes, over the same time window as the raster. A `.spk` file does not contain the continuous signal, only a 3.04 ms cut-out around each detected spike (the continuous signal is in AxIS `.raw` files). The trace draws each cut-out at its true time, coloured by unit; the dotted baseline marks stretches that were not saved.
  - **Drag** to zoom; the zoom is shared with the raster and the charts. Zoomed in far enough, individual samples are shown.
  - **Y scale** is automatic (shared across lanes) or fixed at ±25 to ±500 µV.
  - **Hover** a cut-out for its unit, spike time and the value at that sample; **click** it to open it in the Spike browser.
- **Plate map:** the mean firing rate of every electrode, laid out on the real well/electrode grid. Click an electrode to select it.
- **Spike raster:** every recorded electrode, grouped by well. Drag to zoom in time; click a row to select that electrode.
- **Selected electrode:** waveform overlay with the mean, firing rate over time, log-scale ISI histogram (with the % of intervals under 2 ms), and trigger amplitude over time. All four follow the zoom window.
- **Waveforms by electrode:** one small panel per recorded electrode, laid out on each well's real electrode grid. Each panel overlays up to N spikes from the current window (25 to 1000), aligned on their trigger, with the mean in blue. The y-axis can be shared, for comparing amplitudes, or scaled per electrode, for seeing shapes. Hover a panel for peak-to-peak and min/max of the mean; click it to select that electrode.
- **Electrodes table:** sortable, and click a row to select it. **File info:** header fields, metadata tags, and the AxIS recording and spike-detector settings.
- **Spike browser:** every spike of the selected electrode across the whole recording (time × trigger amplitude, peak-to-peak, PC1 or PC2), coloured by unit once sorted.
  - **Click** a spike to see its waveform over its unit's mean ± SD, with time, unit, amplitudes, noise and the intervals to the neighbouring spikes.
  - The waveform plot is labelled with the spike's unit and time (m:ss.sss and seconds).
  - **Drag** across time, or type a start and end time in the **Range** boxes (seconds like `385.2`, or `m:ss` like `6:25.2`; Enter applies, empty boxes clear), to overlay every spike in that range, coloured by unit and counted per unit. The boxes follow a drag selection.
  - In a range, **hover** an overlaid waveform to highlight it and see its unit and time; **click** it to open that spike.
  - **◀ Prev / Next ▶** step through spikes. Unticking a unit in the legend hides it and skips it when stepping.
  - **Show in raster** zooms the raster and the charts to the selected spike (±0.5 s) or range.
- **Spike sorting:** sort the selected electrode or the whole plate. Each unit gets a colour, and the panel shows the PCA feature space in 2D (PC1–PC3 pairs) and in a rotatable 3D view (drag to rotate, **Spin**, double-click to reset; all three axes share one scale so the cloud keeps its true shape), unit templates (mean ± SD) and a unit table (spikes, rate, peak-to-peak, SNR, ISI CV, refractory removals). You can tick units to **merge** them or **mark them as noise**. After sorting, the raster colours the selected electrode's spikes by unit, the waveform grid draws one mean per unit, and both CSV exports gain a `unit` column. **Export units (CSV)** writes one row per unit. See [Spike sorting](#spike-sorting) below.
- **CSV export:** the selected electrode's spikes in the current window, or every spike in the file (time, electrode, trigger amplitude, waveform min and max, noise SD, threshold multiplier).
- A light/dark theme.

## Spike sorting

The sorter is a JavaScript port of the automatic pipeline of **SAMS** (Ren et al. 2025, bioRxiv [10.1101/2025.02.08.637245](https://doi.org/10.1101/2025.02.08.637245); `src/process_electrode.m`). It uses the same default parameters and sorts each electrode independently:

1. PCA of the spike waveforms.
2. Initial clusters on PC1–PC2. K = 2…5 is chosen by the Davies–Bouldin index.
3. Waveform outliers are removed: spikes whose RMSE to the cluster mean is more than 3 SD from the mean RMSE.
4. Under-sorting fix: Hartigan's dip test is run on PC1 and PC2 of each cluster. Multimodal clusters (p < 0.05) are clustered again, and clusters of 10 spikes or fewer are dropped.
5. Template merge: clusters whose min–max normalised mean waveforms are within a DTW distance of 1.5, and whose peaks and troughs match within 30 % of the height, are merged. This repeats until nothing changes.
6. Final clean-up of each unit:
   - iterative RMSE outlier removal;
   - removal of the second spike of every pair closer than 1.5 ms;
   - units below 0.1 Hz are dropped.

   Electrodes below 0.1 Hz (or under 20 spikes) are not sorted.

Spikes removed along the way are **unsorted** (unit 0). The panel flags "possible multi-unit" when two final templates are close by DTW but were kept apart because their peaks or troughs differ.

**Differences from SAMS:**
- k-means (k-means++ seeding, best of 3) replaces MATLAB's spectral clustering.
- The dip test runs on at most 2,000 spikes per cluster (200-sample bootstrap).
- SAMS's waveform-lag alignment before merging is not ported.
- SAMS's re-classification of outlier spikes is not ported.
- Results are deterministic (seeded) for a given file and parameters.

**Check against the SAMS output shipped with the sample file** (`sample files/SAMS Output Files/spike_sorting.xlsx`):
- 37 of the 40 electrodes SAMS sorted get the same number of units.
- 100 units in total here, against 101 in SAMS.
- Where the unit count matches, unit sizes differ by 2.9 % on average, and most electrodes match spike for spike (e.g. A6_11: 6341 + 2051).
- The differences are on the most crowded electrodes: A6_33, B4_31, C3_33 and F5_23.
- Sorting the whole plate (181,548 spikes, 47 electrodes) takes about 4 s in the browser and about 9 s with `axion_sort.py`.

**Python and browser give identical results.** `axion_sort.py` uses the same random generator (mulberry32, seeded per electrode), so on the sample file every spike gets the same unit, and the units CSV is byte-for-byte the same as the viewer's.

## Format notes (Spike_v1)

Each spike is one block of `BlockHeaderSize + 2*NumSamplesPerBlock` bytes, starting at `DataRegionStart`:
`int64 startingSample, uint8 channel, uint8 chip, int32 triggerSample, float64 noiseSD (V), float64 thresholdMultiplier, [reserved], int16 data[n]`.
Spike time = `(startingSample + triggerSample) / SamplingFrequency`, and voltage = `data * VoltageScale` (V).
The pair `(chip, channel)` is looked up in the ChannelArray to get the well and electrode position.

Only header version 1.x files are supported, which covers all modern AxIS files. Legacy v0.x files are not.
