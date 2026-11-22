"""Compare axion_spk.read_spk with the output of Axion's own MATLAB loader
(saved by export_with_axion_loader.m), spike by spike and sample by sample.
   python compare_with_axion_loader.py recording.spk axion_ref.mat"""
import os, sys
import numpy as np
from scipy.io import loadmat
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from axion_spk import read_spk

spk, mat = sys.argv[1], sys.argv[2]
sf = read_spk(spk); s = sf.spikes; fs = sf.sampling_frequency
m = {k: np.asarray(v).squeeze() for k, v in loadmat(mat).items() if not k.startswith("__")}
ch = sf.channels
lut = sf._channel_lut()
idx = lut[s["chip"].astype(np.int64) * 256 + s["channel"].astype(np.int64)]
pos = np.stack([ch["WellRow"][idx], ch["WellColumn"][idx], ch["ElectrodeColumn"][idx], ch["ElectrodeRow"][idx]], 1).astype(float)
ref_pos = np.stack([m["wellRow"], m["wellCol"], m["elecCol"], m["elecRow"]], 1).astype(float)
n = len(s)
print(os.path.basename(spk), "- ours:", n, "spikes; Axion loader:", m["startS"].size, "spikes")
checks = {
    "spike count": m["startS"].size == n,
    "sampling frequency": float(m["fs"]) == fs,
    "voltage scale": float(m["voltageScale"]) == sf.voltage_scale,
    "electrode of every spike": np.array_equal(pos, ref_pos),
    "waveform start (s)": np.array_equal(s["startingSample"] / fs, m["startS"]),
    "trigger offset": np.array_equal(s["triggerSample"].astype(float), m["triggerSample"]),
    "noise SD": np.array_equal(s["standardDeviation"], m["stdDev"]),
    "threshold multiplier": np.array_equal(s["thresholdMultiplier"], m["thrMult"]),
    "raw samples": np.array_equal(s["data"], np.atleast_2d(m["data"]).T.astype(np.int16)),
    "voltages (V)": np.array_equal(s["data"].astype(np.float64) * sf.voltage_scale, np.atleast_2d(m["volts"]).T),
    "spike times (LoadAllSpikes)": np.array_equal(sf.spike_times_s, m["allTimes"]),
    "chip / channel": np.array_equal(s["chip"].astype(float), m["allAchk"]) and np.array_equal(s["channel"].astype(float), m["allChannel"]),
}
for k, v in checks.items():
    print(f"  {k:28s} {'identical' if v else 'DIFFERENT'}")
print("  RESULT:", "IDENTICAL to Axion's loader" if all(checks.values()) else "DIFFERENCES FOUND")
