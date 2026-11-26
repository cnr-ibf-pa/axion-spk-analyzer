"""Compare the network bursts detected by mea_pipeline.html (dumped by dump_page_nb.js) with network_bursts.py.
   node dump_page_nb.js .. recording.spk /tmp/nb.json [t0 t1]
   python compare_network_bursts.py recording.spk /tmp/nb.json [t0 t1]"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from axion_spk import read_spk
from network_bursts import METHODS, network_bursts

spk, js = sys.argv[1], sys.argv[2]
t0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
t1 = float(sys.argv[4]) if len(sys.argv) > 4 else None
page = json.load(open(js))
sf = read_spk(spk)
py = network_bursts(sf, t0, t1)
bad, worst, nb_total = [], 0.0, {m: 0 for m in METHODS}
def close(a, b, what):
    global worst
    if a is None or (isinstance(b, float) and math.isnan(b)):
        if not (a is None and isinstance(b, float) and math.isnan(b)): bad.append(f"{what}: page {a} python {b}")
        return
    rel = abs(a - b) / max(1e-300, abs(b)) if b else abs(a)
    worst = max(worst, rel)
    if rel > 1e-9: bad.append(f"{what}: page {a} python {b}")
py = {w: r for w, r in py.items() if r["active"]}
if set(page) != set(py): bad.append(f"wells differ: page {sorted(page)} python {sorted(py)}")
for w in sorted(set(page) & set(py)):
    p, q = page[w], py[w]
    if (p["active"], p["spikes"]) != (q["active"], q["spikes"]): bad.append(f"{w}: active/spikes page {p['active']}/{p['spikes']} python {q['active']}/{q['spikes']}")
    for m in METHODS:
        A, B = p["bursts"][m], q["bursts"][m]
        nb_total[m] += len(B)
        if len(A) != len(B): bad.append(f"{w} {m}: page {len(A)} python {len(B)} network bursts"); continue
        for i, (x, y) in enumerate(zip(A, B)):
            if x["start"] != y["start"] or x["end"] != y["end"] or x["spikes"] != y["spikes"] or x["electrodes"] != y["electrodes"]:
                bad.append(f"{w} {m} #{i + 1}: page {x} python {y}"); break
        for k, kk in [("n", "n"), ("perMin", "per_min"), ("durMedian", "dur_median"), ("durMean", "dur_mean"), ("ibiMean", "ibi_mean"), ("ibiCv", "ibi_cv"),
                      ("spikesMedian", "spikes_median"), ("electrodesMedian", "electrodes_median"), ("pctSpikes", "pct_spikes")]:
            close(p["summary"][m][k], q["summary"][m][kk], f"{w} {m} {k}")
        for m2 in METHODS: close(p["agreement"][m][m2], q["agreement"][m][m2], f"{w} agreement {m}->{m2}")
    # ISI_N histogram (exact counts) and the proposed valley threshold
    if p["hist"] != [float(x) for x in q["isin_hist"]]: bad.append(f"{w}: ISI_N histograms differ")
    close(p["valley"], q["isin_valley"], f"{w} ISI_N valley")
print(f"{os.path.basename(spk)}: {len(py)} wells with active electrodes (period {t0} – {t1 or 'end'} s)")
print("  network bursts compared:", ", ".join(f"{m} {n}" for m, n in nb_total.items()))
print(f"  summaries and agreement: largest relative difference {worst:.2e}")
print("  RESULT:", "IDENTICAL (bursts exact, summaries within 1e-9)" if not bad else f"{len(bad)} DIFFERENCES: " + "; ".join(bad[:5]))
