"""Compare the firing statistics computed by mea_pipeline.html (dumped by dump_page_stats.js)
with firing_stats.py, electrode by electrode.
   node dump_page_stats.js ../mea_pipeline.html recording.spk /tmp/page.json [t0 t1]
   python compare_firing_stats.py recording.spk /tmp/page.json [t0 t1]"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from axion_spk import read_spk
from firing_stats import firing_stats

spk, js = sys.argv[1], sys.argv[2]
t0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
t1 = float(sys.argv[4]) if len(sys.argv) > 4 else None
page = json.load(open(js))
el, _ = firing_stats(read_spk(spk), t0, t1)
pairs = [("N", "N"), ("rate", "rate"), ("active", "active"), ("fano", "fano"), ("nIsi", "n_isi"), ("meanIsi", "mean_isi"),
         ("medIsi", "median_isi"), ("cv", "cv"), ("cv2", "cv2"), ("lv", "lv"), ("shortPct", "short_pct")]
worst = {p: 0.0 for p, _ in pairs}; bad = []
for lab, s in page.items():
    r = el[lab]
    for pk, rk in pairs:
        a, b = s[pk], r[rk]
        if a is None or (isinstance(b, float) and math.isnan(b)):
            if not (a is None and isinstance(b, float) and math.isnan(b)): bad.append(f"{lab} {pk}: page {a} python {b}")
            continue
        rel = abs(float(a) - float(b)) / max(1e-300, abs(float(b))) if b else abs(float(a))
        worst[pk] = max(worst[pk], rel)
        if rel > 1e-9: bad.append(f"{lab} {pk}: page {a} python {b}")
print(f"{os.path.basename(spk)}: {len(page)} electrodes with spikes compared (period {t0} – {t1 or 'end'} s)")
for k, v in worst.items():
    print(f"  {k:9s} largest relative difference {v:.2e}")
print("  RESULT:", "IDENTICAL (within 1e-9)" if not bad else f"{len(bad)} DIFFERENCES: " + "; ".join(bad[:5]))
