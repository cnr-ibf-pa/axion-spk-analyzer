"""
network_bursts.py - network bursts of each well of an Axion .spk file with three methods, with the
same definitions as the "Network bursts" panel of mea_pipeline.html.

Common input (per well, period [a, a + T]):
  E   the well's active electrodes (firing_stats: 60 r >= theta), A = |E|
  s_1 <= ... <= s_M  the pooled spike train of E (ties ordered by electrode index)
  m   minimum electrodes in a network burst (default 3)

1. ISI_N (Bakkum et al. 2014): window i = spikes i .. i+N-1 qualifies when s_{i+N-1} - s_i <= theta_N.
   Qualifying windows that share at least one spike are merged; a network burst is the run of spikes they
   cover, kept if it has >= min_spikes spikes on >= m electrodes.
2. Burst grouping (after Chiappalone et al. 2006, Pasquale et al. 2010): single-electrode bursts are
   maximal runs of >= n_b spikes of one electrode with every ISI <= ISI_max. Their onsets (all electrodes)
   are sorted; consecutive onsets <= tau apart form a group, kept if >= max(m, ceil(f A)) distinct
   electrodes have a burst in it. Network burst = [first onset, last end of its bursts].
3. Recruitment: K = floor(T / b) bins [a + k b, a + (k+1) b); bin k qualifies when >= max(m, ceil(f A))
   electrodes of E fire in it. Runs of consecutive qualifying bins are network bursts, from their first
   to their last spike.
Then, for every method, network bursts whose intervals overlap are merged.

Per network burst [start, end]: spikes and electrodes of E with start <= t <= end.

    python network_bursts.py recording.spk [--csv bursts.csv] [--summary summary.csv] [--isin-thr 25] ...
"""
from __future__ import annotations

import argparse
import math

import numpy as np

from firing_stats import DEFAULTS as FS_DEFAULTS, firing_stats

DEFAULTS = dict(min_el=3, isin_n=10, isin_thr=0.025, isin_min_spikes=10,
                burst_isi=0.1, burst_min_spikes=5, onset_gap=0.1, group_pct=20.0,
                rec_bin=0.05, rec_pct=25.0)
METHODS = ("isin", "group", "recruit")


def required(pct, A, m):
    return max(m, math.ceil(pct / 100 * A - 1e-9))


def pooled(trains):
    """trains: list of (electrode index, ascending times). Returns times s and electrodes e, sorted by (t, e)."""
    if not trains:
        return np.zeros(0), np.zeros(0, np.int64)
    t = np.concatenate([x for _, x in trains])
    e = np.concatenate([np.full(len(x), c, np.int64) for c, x in trains])
    o = np.lexsort((e, t))
    return t[o], e[o]


def isin_values(s, N):
    return s[N - 1:] - s[:len(s) - N + 1] if len(s) >= N else np.zeros(0)


def nb_isin(s, e, N, thr, min_spikes, m):
    d = isin_values(s, N)
    out, i0, iend = [], -1, -1
    def close():
        if i0 >= 0 and iend - i0 + 1 >= min_spikes and len(np.unique(e[i0:iend + 1])) >= m:
            out.append((s[i0], s[iend]))
    for i in np.flatnonzero(d <= thr):
        if i0 >= 0 and i <= iend:
            iend = i + N - 1
        else:
            close()
            i0, iend = i, i + N - 1
    close()
    return out


def single_bursts(t, isi_max, n_min):
    """Maximal runs of >= n_min consecutive spikes with every ISI <= isi_max: list of (onset, end)."""
    out, j = [], 0
    n = len(t)
    for i in range(1, n + 1):
        if i == n or t[i] - t[i - 1] > isi_max:
            if i - j >= n_min:
                out.append((t[j], t[i - 1]))
            j = i
    return out


def nb_group(trains, isi_max, n_min, gap, req):
    b = [(on, c, end) for c, t in trains for on, end in single_bursts(t, isi_max, n_min)]
    b.sort(key=lambda x: (x[0], x[1]))
    out, g = [], []
    def close():
        if g and len({c for _, c, _ in g}) >= req:
            out.append((g[0][0], max(x[2] for x in g)))
    for x in b:
        if g and x[0] - g[-1][0] > gap:
            close(); g = []
        g.append(x)
    close()
    return out


def recruitment(trains, a, T, b):
    K = math.floor(T / b + 1e-9)
    rec = np.zeros(K, np.int64)
    for _, t in trains:
        k = np.floor((t - a) / b).astype(np.int64)
        rec[np.unique(k[(k >= 0) & (k < K)])] += 1
    return rec


def nb_recruit(s, a, T, b, trains, req):
    rec = recruitment(trains, a, T, b)
    q = rec >= req
    k = np.floor((s - a) / b).astype(np.int64)
    out, r = [], 0
    K = len(q)
    while r < K:
        if not q[r]:
            r += 1; continue
        r1 = r
        while r1 + 1 < K and q[r1 + 1]:
            r1 += 1
        inside = np.flatnonzero((k >= r) & (k <= r1))
        if len(inside):
            out.append((s[inside[0]], s[inside[-1]]))
        r = r1 + 1
    return out


def merge(nb):
    out = []
    for st, en in sorted(nb):
        if out and st <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], en))
        else:
            out.append((st, en))
    return out


def describe(nb, s, e):
    rows = []
    for st, en in nb:
        i0, i1 = np.searchsorted(s, st, "left"), np.searchsorted(s, en, "right")
        rows.append(dict(start=float(st), end=float(en), duration=float(en - st), spikes=int(i1 - i0), electrodes=int(len(np.unique(e[i0:i1])))))
    return rows


def summary(rows, M, T):
    n = len(rows)
    dur = np.array([r["duration"] for r in rows]); sp = np.array([r["spikes"] for r in rows]); ne = np.array([r["electrodes"] for r in rows])
    ibi = np.diff([r["start"] for r in rows])
    med = lambda x: float(np.median(x)) if len(x) else math.nan
    mibi = float(ibi.mean()) if len(ibi) else math.nan
    return dict(n=n, per_min=60 * n / T, dur_median=med(dur), dur_mean=float(dur.mean()) if n else math.nan,
                ibi_mean=mibi, ibi_cv=float(ibi.std(ddof=1) / mibi) if len(ibi) > 1 else math.nan,
                spikes_median=med(sp), electrodes_median=med(ne), pct_spikes=100 * float(sp.sum()) / M if M else math.nan)


def agreement(x, y):
    """Percentage of the network bursts in x that overlap at least one in y."""
    if not x:
        return math.nan
    ys = np.array([b[0] for b in y]); ye = np.array([b[1] for b in y])
    hit = sum(1 for st, en in x if np.any((ys <= en) & (ye >= st)))
    return 100 * hit / len(x)


def isin_valley(s, N, per_decade=10, lo=1e-4, hi=100.0):
    """Threshold proposed from the histogram of log10 ISI_N (per_decade bins from lo to hi, no smoothing):
    p = fullest bin; for each later bin v, depth(v) = min(h[p], max_{k>v} h[k]) - h[v]; the deepest v wins
    (first one on ties) and the threshold is its upper edge. Returns (threshold or nan, histogram)."""
    d = isin_values(s, N)
    nb = round(math.log10(hi / lo) * per_decade)
    h = np.zeros(nb, np.int64)
    with np.errstate(divide="ignore"):
        k = np.floor(np.log10(d / lo) * per_decade + 1e-9)
    ok = (d > 0) & (k >= 0) & (k < nb)
    np.add.at(h, k[ok].astype(np.int64), 1)
    if not h.any():
        return math.nan, h
    p = int(np.argmax(h))
    best, bv = 0, -1
    for v in range(p + 1, nb - 1):
        depth = min(h[p], h[v + 1:].max()) - h[v]
        if depth > best:
            best, bv = depth, v
    return (lo * 10 ** ((bv + 1) / per_decade) if bv >= 0 else math.nan), h


def well_network_bursts(trains, a, T, **params):
    p = {**DEFAULTS, **params}
    A = len(trains)
    s, e = pooled(trains)
    nb = dict(isin=nb_isin(s, e, p["isin_n"], p["isin_thr"], p["isin_min_spikes"], p["min_el"]),
              group=nb_group(trains, p["burst_isi"], p["burst_min_spikes"], p["onset_gap"], required(p["group_pct"], A, p["min_el"])),
              recruit=nb_recruit(s, a, T, p["rec_bin"], trains, required(p["rec_pct"], A, p["min_el"])))
    nb = {k: merge(v) for k, v in nb.items()}
    rows = {k: describe(v, s, e) for k, v in nb.items()}
    valley, hist = isin_valley(s, p["isin_n"])
    return dict(active=A, spikes=len(s), bursts=rows, isin_valley=valley, isin_hist=hist, summary={k: summary(rows[k], len(s), T) for k in METHODS},
                agreement={x: {y: agreement(nb[x], nb[y]) for y in METHODS} for x in METHODS})


def network_bursts(sf, t0=0.0, t1=None, active_per_min=FS_DEFAULTS["active_per_min"], **params):
    """Network bursts of every recorded well of a SpikeFile (axion_spk.read_spk) over [t0, t1]."""
    duration = sf.header.get("Duration") or float(sf.spike_times_s.max())
    t1 = duration if t1 is None else t1
    a, T = t0, t1 - t0
    el, _ = firing_stats(sf, t0, t1, active_per_min=active_per_min)
    lab_all = [sf._label(c) for c in sf.channels]
    labels, times = sf.electrode_labels(), sf.spike_times_s
    wells = {}
    for lab in el:
        wells.setdefault(lab.split("_")[0], []).append(lab)
    out = {}
    for w, labs in sorted(wells.items()):
        trains = []
        for lab in sorted(labs, key=lab_all.index):
            if not el[lab]["active"]:
                continue
            t = np.sort(times[labels == lab])
            trains.append((lab_all.index(lab), t[(t >= a) & (t <= a + T)]))
        out[w] = well_network_bursts(trains, a, T, **params)
    return out


def _main(argv=None):
    from axion_spk import read_spk
    ap = argparse.ArgumentParser(description="Network bursts of an Axion .spk file (same definitions as mea_pipeline.html)")
    ap.add_argument("spk")
    ap.add_argument("--start", type=float, default=0.0); ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--active", type=float, default=5.0, help="active electrode threshold, spikes/min")
    ap.add_argument("--min-el", type=int, default=3)
    ap.add_argument("--isin-n", type=int, default=10); ap.add_argument("--isin-thr", type=float, default=25, help="ms")
    ap.add_argument("--isin-min-spikes", type=int, default=10)
    ap.add_argument("--burst-isi", type=float, default=100, help="ms"); ap.add_argument("--burst-min-spikes", type=int, default=5)
    ap.add_argument("--onset-gap", type=float, default=100, help="ms"); ap.add_argument("--group-pct", type=float, default=20)
    ap.add_argument("--rec-bin", type=float, default=50, help="ms"); ap.add_argument("--rec-pct", type=float, default=25)
    ap.add_argument("--csv", help="write every network burst to this CSV")
    ap.add_argument("--summary", help="write per-well, per-method summaries to this CSV")
    a = ap.parse_args(argv)
    res = network_bursts(read_spk(a.spk), a.start, a.end, active_per_min=a.active, min_el=a.min_el, isin_n=a.isin_n,
                         isin_thr=a.isin_thr / 1000, isin_min_spikes=a.isin_min_spikes, burst_isi=a.burst_isi / 1000,
                         burst_min_spikes=a.burst_min_spikes, onset_gap=a.onset_gap / 1000, group_pct=a.group_pct,
                         rec_bin=a.rec_bin / 1000, rec_pct=a.rec_pct)
    for w, r in res.items():
        print(f"{w}: {r['active']} active electrodes, {r['spikes']} spikes")
        for k in METHODS:
            s = r["summary"][k]
            print(f"  {k:8s} {s['n']:5d} network bursts  {s['per_min']:.2f}/min  median duration {s['dur_median'] * 1000:.0f} ms  {s['pct_spikes']:.1f} % of spikes")
    f = lambda v, k=1: "" if not np.isfinite(v) else f"{v * k:.6f}"
    if a.csv:
        with open(a.csv, "w") as fh:
            fh.write("well,method,index,start_s,end_s,duration_ms,spikes,electrodes\n")
            for w, r in res.items():
                for k in METHODS:
                    for i, b in enumerate(r["bursts"][k], 1):
                        fh.write(f"{w},{k},{i},{b['start']:.6f},{b['end']:.6f},{b['duration'] * 1000:.3f},{b['spikes']},{b['electrodes']}\n")
        print(f"wrote {a.csv}")
    if a.summary:
        with open(a.summary, "w") as fh:
            fh.write("well,method,active,n,per_min,dur_median_ms,dur_mean_ms,ibi_mean_s,ibi_cv,spikes_median,electrodes_median,pct_spikes\n")
            for w, r in res.items():
                for k in METHODS:
                    s = r["summary"][k]
                    fh.write(",".join([w, k, str(r["active"]), str(s["n"]), f(s["per_min"]), f(s["dur_median"], 1000), f(s["dur_mean"], 1000),
                                       f(s["ibi_mean"]), f(s["ibi_cv"]), f(s["spikes_median"]), f(s["electrodes_median"]), f(s["pct_spikes"])]) + "\n")
        print(f"wrote {a.summary}")


if __name__ == "__main__":
    _main()
