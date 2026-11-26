# PCA and spike sorting: what the former SPK viewer did

`spk_viewer.html` was removed from the web app on 2026-11-23. Its PCA and spike sorting are the starting point for the next step of the MEA pipeline, so this note records what they did and how they were checked.

Where the code is now:
- [`pca_sorting_reference.js`](pca_sorting_reference.js): the viewer's PCA, sorting, curation, drawing and CSV code, copied verbatim (lines 817–1365). No page loads it. It still reads the old page's globals (`F`, `S`, `SS`, `$`, `C`, `fmt`, `setup`, `frame`, `niceTicks`, `pctl`, `chartHover`, `download`, `uV`). The algorithm parts (`eigSym`, `pcaScores`, `kmeans`, `daviesBouldin`, `clusterBest`, `dipStat`, `dipTest`, `rmseOutliers`, `meanWave`, `dtw`, `mergePlan`) are pure functions. Only `sortElectrode` and `unitStats` read `F`, and they would need the recording `R` as an argument, as was done for `js/analysis.js`.
- [`python/axion_sort.py`](../python/axion_sort.py): the Python/numpy twin of the same sorter. It is still in the repository and does not depend on the viewer.
- The whole page: `git show 6fc708d:spk_viewer.html > spk_viewer.html`, the last committed version, identical to the one removed.

## PCA as implemented

Each electrode is handled on its own, using every spike of the whole recording; there was no time-window option.
- **Input:** the waveforms exactly as stored. That is `NumSamplesPerBlock` samples per spike (38 on the test files, 3.04 ms), in µV (sample × `VoltageScale`), aligned on the file's stored trigger. There is no re-alignment, interpolation, filtering or per-spike normalisation.
- **Covariance:** the sample covariance of those vectors, centred on the mean waveform and divided by n − 1.
- **Eigen-decomposition:** cyclic Jacobi (`eigSym`), with eigenvalues sorted in descending order. Python uses `np.linalg.eigh`.
- **Sign convention:** each of the first 3 eigenvectors is flipped so that its largest-|loading| element is positive. This makes the scores reproducible.
- **Output:** scores on PC1–PC3, `(x − mean) · v_k`, and the explained variance of each, `λ_k / Σ max(λ, 0)`. These are kept per electrode in `res.pcs` (n × 3) and `res.explained`.
- PCA runs only as part of sorting. An electrode with fewer than `max(20, duration × minRate)` spikes was neither sorted nor PCA'd.

## Sorting pipeline (port of SAMS, Ren et al. 2025, `src/process_electrode.m`)

1. PCA, as above.
2. Initial clusters on PC1–PC2. K-means is run for K = 2…maxK (k-means++ seeding, best of 3 runs by within-cluster sum of squares, mulberry32 RNG seeded `0x5A15 + channel`), and K is chosen by the lowest Davies–Bouldin index. SAMS uses spectral clustering here.
3. Outliers removed per cluster: spikes whose waveform RMSE to the cluster mean is more than `sdCut` SD from the mean RMSE.
4. Hartigan's dip test (port of AS 217) on PC1 and PC2 of each cluster. Its p-value comes from a 200-sample uniform bootstrap, on at most 2,000 spikes. A cluster with p < 0.05 on either axis is clustered again as in step 2, and clusters of 10 spikes or fewer are dropped.
5. Template merge, repeated until nothing changes. All templates are min–max normalised together, and two are merged when their DTW distance (sum of |a − b| along the path, as in MATLAB `dtw`) is below `mergeThr` and their peaks and troughs lie within (1 − overlap) × height. Merges are applied as connected components.
6. Clean-up of each unit:
   - iterative RMSE outlier removal;
   - the second spike of every pair closer than the refractory period is removed;
   - units with rate ≤ `minRate` are dropped.

Spikes removed along the way are "unsorted" (unit 0). Units are numbered by size. A "possible multi-unit" flag is raised when two final templates are within `mergeThr` by DTW but their peaks or troughs differ.

Defaults (SAMS): `mergeThr` 1.5, `overlap` 0.7, refractory 1.5 ms, `sdCut` 3, `minRate` 0.1 Hz, `maxK` 5.

Not ported from SAMS:
- the waveform-lag alignment before merging;
- the re-classification of outlier spikes.

## Unit statistics and exports

Each unit had:
- n, rate;
- mean and SD template;
- peak-to-peak, and peak (the extremum with the largest absolute value);
- noise (mean of the detector noise SD over its spikes);
- SNR = peak-to-peak / noise;
- ISI CV, minimum ISI, and the number of refractory removals.

**Export units (CSV)** wrote one row per unit with these columns: `electrode, well, treatment, unit, spikes, rate_hz, peak_to_peak_uV, peak_uV, snr, isi_cv, min_isi_ms, refractory_removed, unsorted_on_electrode`. `axion_sort.py --units` writes the same file.

## Views

- **Feature space 2D:** PC1×PC2, PC1×PC3 or PC2×PC3, coloured by unit.
  - Display only: at most 6,000 points, evenly subsampled, and axes clipped to the 0.2–99.8 % percentiles.
  - Axis labels show the % of variance.
  - Hovering shows the unit and spike time.
- **Feature space 3D:** PC1×PC2×PC3, orthographic, drag to rotate, Spin, double-click to reset. All three axes share one scale, so the cloud keeps its true shape. Near points are drawn slightly larger.
- **Unit templates:** mean ± SD per unit.
- **Unit table:** hover a row to highlight that unit in every view. Tick units to **Merge** them or **Mark as noise**; the result is flagged "edited by hand".
- **Spike browser:** the Y axis could be PC1 or PC2 once the electrode was sorted.

## How it was checked

- **Against SAMS's own output** for the 96-well sample (`sample files/SAMS Output Files/spike_sorting.xlsx`):
  - 37 of 40 electrodes get the same number of units, with 100 units in total against 101;
  - where the counts match, unit sizes differ by 2.9 % on average, and most electrodes match spike for spike (e.g. A6_11: 6341 + 2051);
  - the differences are on the most crowded electrodes (A6_33, B4_31, C3_33, F5_23).
- **Browser vs Python:** the same unit for every spike on the sample file, and a byte-identical units CSV.
- **Speed:** the whole plate (181,548 spikes, 47 electrodes) took about 4 s in the browser and 9 s in Python.
- **Not checked:** the PCA scores on their own were never compared with an independent reference such as `sklearn.decomposition.PCA`. They were checked only through the identical sorting.

## Points to decide for the pipeline's PCA step

These choices were made in the viewer but never written into a definition box:
- Which samples go into PCA: all stored samples, or a window around the trigger.
- Alignment: the stored trigger, or the same alignment options as the Waveforms panel.
- Covariance or correlation, i.e. whether each sample is scaled.
- Whole recording, or the Firing-statistics period.
- Per electrode, as before, or pooled per well.
- How many PCs to keep, and whether to report explained variance per electrode.
- Display subsampling: the 6,000-point subsample was only for drawing. The pipeline's rule of "no hidden smoothing" means any such subsample must be stated on the plot.
