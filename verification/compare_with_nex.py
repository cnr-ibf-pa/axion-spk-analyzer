"""Compare axion_spk.read_spk with an Axion Data Export Tool .nex of the same recording.
   python compare_with_nex.py recording.spk recording.nex"""
import sys, numpy as np
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from axion_spk import read_spk
from nexread import read_nex

spk_path, nex_path = sys.argv[1], sys.argv[2]
sf = read_spk(spk_path)
h, nv = read_nex(nex_path)
fs = sf.sampling_frequency
labels = sf.electrode_labels()
start = sf.spikes["startingSample"].astype(np.int64)
trig = sf.spikes["triggerSample"].astype(np.int64)
volts = sf.spikes["data"].astype(np.float64) * sf.voltage_scale          # V, exactly as Axion's loader (data * VoltageScale)

neur = {v["name"]: v for v in nv if v["type"] == "neuron"}
wave = {v["name"][:-3]: v for v in nv if v["type"] == "waveform"}
ours = {l: np.nonzero(labels == l)[0] for l in np.unique(labels)}
print(f"file: {spk_path.split('/')[-1]}  spikes: {len(labels)}  NEX clock {h['freq']} Hz, .spk fs {fs} Hz")

nex_active = {k for k, v in neur.items() if v["count"] > 0}
print("electrodes with spikes: ours", len(ours), " NEX", len(nex_active), " same set:", set(ours) == nex_active)
bad = []
tot = dict(spikes=0, ts_trig=0, ts_start=0, wf_ok=0, maxerr_steps=0.0, maxerr_uV=0.0)
for lab, idx in sorted(ours.items()):
    o = idx[np.argsort(start[idx] + trig[idx], kind="stable")]
    n_ts = neur[lab]["ts"]; w = wave[lab]
    if len(n_ts) != len(o):
        bad.append(f"{lab}: count ours {len(o)} vs NEX {len(n_ts)}"); continue
    tot["spikes"] += len(o)
    tot["ts_trig"] += int(np.sum(n_ts == start[o] + trig[o]))
    tot["ts_start"] += int(np.sum(n_ts == start[o]))
    if not np.array_equal(w["ts"], n_ts):
        bad.append(f"{lab}: waveform timestamps differ from neuron timestamps")
    nexV = (w["wf_raw"].astype(np.float64) * w["adtomv"] + w["mvoffset"]) / 1e3   # NEX stores millivolts: raw*ADtoMV + offset
    err = np.abs(nexV - volts[o])
    steps = err / (w["adtomv"] / 1e3)
    tot["maxerr_steps"] = max(tot["maxerr_steps"], float(steps.max()))
    tot["maxerr_uV"] = max(tot["maxerr_uV"], float(err.max() * 1e6))
    tot["wf_ok"] += int(np.sum(np.all(steps <= 0.5 + 1e-9, axis=1)))
print("compared spikes:", tot["spikes"])
print("timestamps equal to start+trigger sample:", tot["ts_trig"], " equal to start sample:", tot["ts_start"])
print("waveforms equal within NEX re-quantisation (<= 0.5 step):", tot["wf_ok"])
print(f"largest difference: {tot['maxerr_steps']:.3f} NEX steps = {tot['maxerr_uV']:.4f} uV  (.spk resolution {abs(sf.voltage_scale)*1e6:.4f} uV)")
print("problems:", bad or "none")
