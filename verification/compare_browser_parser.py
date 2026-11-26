"""Compare the browser page's parser (dumped by dump_browser_parser.js) with axion_spk.read_spk, spike by spike.
   node dump_browser_parser.js .. recording.spk /tmp/out
   python compare_browser_parser.py recording.spk /tmp/out"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from axion_spk import read_spk

spk, pre = sys.argv[1], sys.argv[2]
sf = read_spk(spk); s = sf.spikes; h = json.load(open(pre + ".hdr.json"))
n, ns = h["n"], h["ns"]
order = np.argsort(sf.spike_times_s, kind="stable") if h["sortedOnLoad"] else np.arange(len(s))   # page re-sorts unordered files
js = dict(times=np.fromfile(pre + ".times", "<f8"), sd=np.fromfile(pre + ".sd", "<f8"), thr=np.fromfile(pre + ".thr", "<f8"),
          trig=np.fromfile(pre + ".trig", "<i4"), wf=np.fromfile(pre + ".wf", "<i2").reshape(n, ns),
          lab=np.array(open(pre + ".lab").read().split("\n")))
py = dict(times=sf.spike_times_s[order], sd=s["standardDeviation"][order], thr=s["thresholdMultiplier"][order],
          trig=s["triggerSample"][order].astype(np.int32), wf=s["data"][order], lab=sf.electrode_labels()[order])
print(os.path.basename(spk))
print("  header:", {k: v for k, v in dict(n=n == len(s), fs=h["fs"] == sf.sampling_frequency, VoltageScale=h["scale"] == sf.voltage_scale,
      duration=h["duration"] == sf.header["Duration"], plate=h["plateType"] == sf.plate_type, channels=h["nch"] == len(sf.channels)).items()})
ok = True
for k in ("times", "trig", "lab", "wf", "sd", "thr"):
    a, b = js[k], py[k]
    eq = np.array_equal(a, b)
    note = ""
    if not eq and k in ("sd", "thr") and np.array_equal(a, b.astype(np.float32).astype(np.float64)):
        eq, note = True, " (equal after float32 rounding: the page stores this field in single precision)"
    ok &= eq
    print(f"  {k:6s} identical: {eq}{note}" + ("" if eq else f"  max diff {np.max(np.abs(a.astype(float) - b.astype(float)))}"))
print("  RESULT:", "IDENTICAL" if ok else "DIFFERENT")
