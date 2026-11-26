# Verification of the .spk reader

Checks that our readers return exactly what Axion's own software reads from a `.spk` file.

| Check | Reference | Script |
|---|---|---|
| Python reader vs Axion **Data Export Tool** | `.nex` exported by Axion from the same `.spk` | `compare_with_nex.py recording.spk recording.nex` |
| Browser parser vs Python reader | — | `node dump_browser_parser.js .. recording.spk /tmp/out` (runs `js/spk_reader.js`; `..` is the app folder) then `compare_browser_parser.py recording.spk /tmp/out` |
| Python reader vs Axion's **MATLAB loader** (AxionFileLoader) | `.mat` written by Axion's loader | in MATLAB: `export_with_axion_loader('<AxionFileLoader folder>', 'recording.spk', 'axion_ref.mat')`, then `compare_with_axion_loader.py recording.spk axion_ref.mat` |

`nexread.py` is a minimal reader of the documented NeuroExplorer (.nex) format, used by `compare_with_nex.py`.

## Results (2026-11-11)

**Against Axion's Data Export Tool (.nex)**:

| Recording | Spikes | Electrodes | Result |
|---|---|---|---|
| `2-3-22_SampleElec_15min` (SAMS sample, 96-well) | 181,548 | 47 | all spikes match |
| `hNEU_CYC_SAG_DIV98_PTX_APV_DNQX_031126` (CytoView MEA 6) | 151,980 | 55 | all spikes match |

For each recording:
- the electrodes with spikes and the spike count of every electrode are the same;
- the timestamps are identical, sample for sample. NEX stores the first sample of each waveform (`startingSample`).
- every waveform sample is the same within the export's own rounding. The `.nex` stores millivolts, rescaled per electrode, and its largest difference from our values is 0.0018 µV, 30 times smaller than the 0.055 µV resolution of the `.spk` file.

The export's scaling equals the file's `VoltageScale`, sign included, which confirms that voltage = sample × `VoltageScale`.

**Browser parser vs Python reader**, on 7 recordings: the SAMS sample and the six CytoView files. This was run on the parser of the former `spk_viewer.html`, from which `js/spk_reader.js` is taken.
- These fields are bit-identical for every spike: time, trigger offset, electrode, every waveform sample and the threshold multiplier.
- The detector noise SD differs only because the page keeps it in single precision (relative error ≈ 1e-7).
- The header values are also equal: n, sampling frequency, VoltageScale, duration, plate type and channel count.

**Against Axion's MATLAB loader**: not run yet, because MATLAB isn't installed on this Mac. The comparison script was tested on a reference built in the same layout: it passes correct data and catches a single changed sample. The MATLAB export script itself is untested.

On 2026-11-23 the check was re-run on `js/spk_reader.js` itself with the test file `DIV98_PTX` (154,934 spikes): identical in every field. The noise SD and threshold multiplier are now equal too, because this reader keeps them in double precision.

## Firing statistics (mea_pipeline.html vs firing_stats.py)

`node dump_page_stats.js .. recording.spk /tmp/page.json [t0 t1]`, then `python3 compare_firing_stats.py recording.spk /tmp/page.json [t0 t1]`. The Node script runs the page's own parser and statistics code (`js/spk_reader.js`, `js/analysis.js`); `..` is the app folder.

The two implementations were checked on 2026-11-11 against each other:

| Recording / period | Electrodes compared | Result |
|---|---|---|
| test file `DIV98_PTX`, whole recording | 52 | identical |
| test file `DIV98_PTX`, 120.5–250.25 s | 52 | identical |
| SAMS 96-well sample, whole recording | 47 | identical |

- Spike counts, rates, active flags, interval counts, medians and short-ISI percentages are exact.
- The Fano factor, mean ISI, CV, CV2 and LV agree to ≤ 2·10⁻¹⁴ relative, which is floating-point summation order.

**Autocorrelogram and ISI histogram:** checked against brute-force counts in Python (all spike pairs within ±100 ms; whole-sample intervals), with zero difference. Intervals that fall exactly on a bin edge (e.g. 10 ms = 125 samples) are put in the bin that starts at that value, as the definition says. numpy's `histogram` on intervals computed in seconds can put them one bin lower, because of rounding.

Spike time convention: our spike time is `(startingSample + triggerSample) / SamplingFrequency`. This is the same formula as `SpikeDataSet.LoadAllSpikes` in Axion's loader. A waveform's `Start` in that loader is `startingSample / SamplingFrequency`, as in the `.nex` export.

## Network bursts (mea_pipeline.html vs network_bursts.py)

`node dump_page_nb.js .. recording.spk /tmp/nb.json [t0 t1]`, then `python3 compare_network_bursts.py recording.spk /tmp/nb.json [t0 t1]`. The Node script runs the page's own parser and the STATS and NB blocks of `js/analysis.js`, with the panel's default parameters; `..` is the app folder.

Checked on 2026-11-11:

| Recording / period | Wells | Network bursts compared (ISI_N / grouping / recruitment) | Result |
|---|---|---|---|
| test file `DIV98_PTX`, whole recording | 2 | 684 / 66 / 376 | identical |
| test file `DIV98_PTX`, 120.5–250.25 s | 2 | 125 / 16 / 57 | identical |
| SAMS 96-well sample, whole recording | 6 | 981 / 536 / 1052 | identical |

- Every network burst has the same start, end, spike count and electrode count (exact equality).
- The ISI_N histograms are identical, and the proposed valley threshold is the same.
- Summaries and agreement percentages agree to ≤ 5·10⁻¹⁶ relative.
