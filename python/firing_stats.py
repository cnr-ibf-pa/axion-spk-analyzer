"""
firing_stats.py - per-electrode and per-well firing statistics of an Axion .spk file,
with the same definitions as the "Firing statistics" panel of mea_pipeline.html.

Definitions (period [a, a + T]; t_1 < ... < t_N spike times of one electrode in it):
  rate      r = N / T                                   (spikes/s)
  active    60 r >= theta                               (theta in spikes/min, default 5)
  Fano      F = s^2 / mean(n_k), n_k counts in K = floor(T/D) whole bins [a+(k-1)D, a+kD),
            s^2 the sample variance (K - 1)
  ISIs      tau_i = t_{i+1} - t_i, n = N - 1 intervals; ISI statistics need n >= min_isi
  CV        sample SD of tau / mean of tau
  CV2       mean over i of 2 |tau_{i+1} - tau_i| / (tau_{i+1} + tau_i)       (Holt et al. 1996)
  LV        3/(n-1) * sum ((tau_i - tau_{i+1}) / (tau_i + tau_{i+1}))^2       (Shinomoto et al. 2003)
  short     100 * #{tau_i < L} / n
Well summaries: median, mean and sample SD over the well's active electrodes.

    python firing_stats.py recording.spk [--csv electrodes.csv] [--wells wells.csv] [--active 5] [--bin 1] ...
"""
from __future__ import annotations

import argparse
import math
import sys

import numpy as np

DEFAULTS = dict(active_per_min=5.0, fano_bin=1.0, short_isi=1.5e-3, min_isi=20)


def electrode_stats(t: np.ndarray, a: float, T: float, active_per_min=5.0, fano_bin=1.0, short_isi=1.5e-3, min_isi=20) -> dict:
    t = np.asarray(t, dtype=np.float64)
    N = len(t)
    rate = N / T
    out = dict(N=N, rate=rate, active=rate * 60 >= active_per_min)
    K = math.floor(T / fano_bin + 1e-9)
    k = np.floor((t - a) / fano_bin).astype(np.int64)
    n = np.bincount(k[(k >= 0) & (k < K)], minlength=K).astype(np.float64) if K > 0 else np.zeros(0)
    m = n.mean() if K else 0.0
    out["fano"] = float(n.var(ddof=1) / m) if K >= 2 and m > 0 else math.nan
    tau = np.diff(t)
    out["n_isi"] = len(tau)
    if len(tau) >= min_isi:
        mu = tau.mean()
        out["mean_isi"] = float(mu)
        out["median_isi"] = float(np.median(tau))
        out["cv"] = float(tau.std(ddof=1) / mu)
        s = tau[:-1] + tau[1:]
        ok = s > 0
        out["cv2"] = float(np.mean(2 * np.abs(tau[1:] - tau[:-1])[ok] / s[ok])) if ok.any() else math.nan
        out["lv"] = float(3 * np.mean(((tau[:-1] - tau[1:])[ok] / s[ok]) ** 2)) if ok.any() else math.nan
        out["short_pct"] = float(100 * np.mean(tau < short_isi))
    else:
        for key in ("mean_isi", "median_isi", "cv", "cv2", "lv", "short_pct"):
            out[key] = math.nan
    return out


def firing_stats(sf, t0: float = 0.0, t1: float | None = None, **params):
    """Statistics of every recorded electrode of a SpikeFile (axion_spk.read_spk) over [t0, t1]
    (default: the whole recording). Returns (electrodes, wells): dicts keyed by label / well name."""
    prm = {**DEFAULTS, **params}
    duration = sf.header.get("Duration") or float(sf.spike_times_s.max())
    if t1 is None:
        t1 = duration
    a, T = t0, t1 - t0
    times = sf.spike_times_s
    labels = sf.electrode_labels()
    # recorded electrodes: the spike header's channel list, plus any electrode with spikes
    lut = {(int(c["ChannelAchk"]), int(c["ChannelIndex"])): i for i, c in enumerate(sf.channels)}
    recorded = {lut[x] for x in sf.channel_ids if x in lut}
    lab_all = [sf._label(c) for c in sf.channels]
    recorded |= {lab_all.index(l) for l in set(labels) if l != "?"}
    order = np.argsort(times, kind="stable")
    el = {}
    for c in sorted(recorded):
        lab = lab_all[c]
        t = times[order][(labels[order] == lab)]
        t = t[(t >= a) & (t <= a + T)]
        el[lab] = electrode_stats(t, a, T, **prm)
    treat = {}
    for tg in sf.tags:
        if "TreatmentWhat" in tg and tg.get("WellRow"):
            treat[f"{chr(64 + tg['WellRow'])}{tg['WellColumn']}"] = tg["TreatmentWhat"]
    wells = {}
    for lab, s in el.items():
        wells.setdefault(lab.split("_")[0], []).append((lab, s))
    wsum = {}
    for w, items in wells.items():
        act = [s for _, s in items if s["active"]]
        row = dict(treatment=treat.get(w, ""), recorded=len(items), active=len(act), well_rate=sum(s["N"] for _, s in items) / T)
        for key in ("rate", "cv", "cv2", "lv", "median_isi", "fano"):
            v = np.array([s[key] for s in act if np.isfinite(s[key])])
            row[key] = dict(n=len(v), median=float(np.median(v)) if len(v) else math.nan, mean=float(v.mean()) if len(v) else math.nan,
                            sd=float(v.std(ddof=1)) if len(v) > 1 else math.nan)
        wsum[w] = row
    return el, wsum


def _main(argv=None):
    from axion_spk import read_spk
    ap = argparse.ArgumentParser(description="Firing statistics of an Axion .spk file (same definitions as mea_pipeline.html)")
    ap.add_argument("spk")
    ap.add_argument("--start", type=float, default=0.0, help="period start (s)")
    ap.add_argument("--end", type=float, default=None, help="period end (s), default: end of recording")
    ap.add_argument("--active", type=float, default=5.0, help="active threshold, spikes/min (default 5)")
    ap.add_argument("--bin", type=float, default=1.0, help="rate / Fano bin width, s (default 1)")
    ap.add_argument("--short", type=float, default=1.5, help="short-ISI limit, ms (default 1.5)")
    ap.add_argument("--min-isi", type=int, default=20, help="minimum intervals for ISI statistics (default 20)")
    ap.add_argument("--csv", help="write per-electrode statistics to this CSV")
    ap.add_argument("--wells", help="write per-well statistics to this CSV")
    a = ap.parse_args(argv)
    sf = read_spk(a.spk)
    el, wells = firing_stats(sf, a.start, a.end, active_per_min=a.active, fano_bin=a.bin, short_isi=a.short / 1000, min_isi=a.min_isi)
    act = sum(s["active"] for s in el.values())
    print(f"{len(el)} recorded electrodes, {act} active")
    for w, r in sorted(wells.items()):
        print(f"  {w} {r['treatment']:8s} active {r['active']}/{r['recorded']}  well rate {r['well_rate']:.2f}/s  "
              f"median rate {r['rate']['median']:.3f} Hz  median CV {r['cv']['median']:.2f}  median LV {r['lv']['median']:.2f}")
    f = lambda v, k=1: "" if not np.isfinite(v) else f"{v * k:.6f}"
    if a.csv:
        with open(a.csv, "w") as fh:
            fh.write("electrode,active,spikes,rate_hz,fano,median_isi_ms,mean_isi_ms,cv,cv2,lv,short_isi_pct,n_isi\n")
            for lab, s in el.items():
                fh.write(",".join([lab, str(int(s["active"])), str(s["N"]), f(s["rate"]), f(s["fano"]), f(s["median_isi"], 1000), f(s["mean_isi"], 1000),
                                   f(s["cv"]), f(s["cv2"]), f(s["lv"]), f(s["short_pct"]), str(s["n_isi"])]) + "\n")
        print(f"wrote {a.csv}")
    if a.wells:
        keys = ("rate", "cv", "cv2", "lv", "median_isi", "fano")
        with open(a.wells, "w") as fh:
            fh.write("well,treatment,active,recorded,well_rate," + ",".join(f"{k}_{s}" for k in keys for s in ("median", "mean", "sd", "n")) + "\n")
            for w, r in sorted(wells.items()):
                fh.write(",".join([w, r["treatment"], str(r["active"]), str(r["recorded"]), f(r["well_rate"])]
                                  + [x for k in keys for x in (f(r[k]["median"]), f(r[k]["mean"]), f(r[k]["sd"]), str(r[k]["n"]))]) + "\n")
        print(f"wrote {a.wells}")


if __name__ == "__main__":
    _main()
