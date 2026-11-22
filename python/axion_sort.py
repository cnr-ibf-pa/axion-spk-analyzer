"""
axion_sort.py - SAMS-style automatic spike sorting for Axion .spk files (numpy only).

A Python port of the sorter in spk_viewer.html, itself a port of the automatic
pipeline of SAMS (Ren et al. 2025, bioRxiv doi:10.1101/2025.02.08.637245,
src/process_electrode.m). It uses the same algorithm, defaults and random number
generator as the browser, so both give the same units for the same file.

Per electrode:
  1. PCA of the waveforms
  2. initial clusters on PC1-PC2: k-means (k-means++ seeding, best of 3) for
     K = 2..max_k, K chosen by the Davies-Bouldin index   [SAMS: spectral clustering]
  3. waveform outliers removed (RMSE to the cluster mean > sd_cut SD from the mean RMSE)
  4. Hartigan's dip test on PC1 and PC2 of each cluster; multimodal ones (p < 0.05)
     are clustered again; clusters of <= 10 spikes are dropped
  5. clusters with similar mean waveforms (DTW on min-max normalised templates <
     merge_thr, peaks and troughs within (1 - overlap) x height) are merged
  6. per unit: iterative RMSE outlier removal, second spike of every pair closer than
     the refractory period removed, units slower than min_rate dropped
Spikes removed along the way are "unsorted" (unit 0).

Usage
-----
    from axion_spk import read_spk
    from axion_sort import sort_spikes
    sf = read_spk("recording.spk")
    res = sort_spikes(sf)                      # every electrode with spikes
    res = sort_spikes(sf, electrodes=["A6_11"])
    res.unit          # unit of every spike: -1 not sorted, 0 unsorted, 1.. unit number
    res.electrodes    # per-electrode results (units with templates and statistics)
    res.write_units_csv("units.csv")

CLI
---
    python axion_sort.py recording.spk --units units.csv --csv spikes.csv
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass, field, asdict

import numpy as np

DEFAULTS = dict(merge_thr=1.5, overlap=0.7, refractory_ms=1.5, sd_cut=3.0, min_rate=0.1, max_k=5)


# --------------------------------------------------------------------------- RNG
class Mulberry32:
    """mulberry32, identical to the browser's mulberry() so results match the viewer."""
    def __init__(self, seed: int):
        self.a = seed & 0xFFFFFFFF

    def __call__(self) -> float:
        self.a = (self.a + 0x6D2B79F5) & 0xFFFFFFFF
        t = self.a
        t = _imul(t ^ (t >> 15), t | 1)
        t ^= (t + _imul(t ^ (t >> 7), t | 61)) & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0


def _imul(a, b):
    return (a * b) & 0xFFFFFFFF


# --------------------------------------------------------------------------- building blocks
def pca_scores(X: np.ndarray, nc: int = 3):
    mu = X.mean(axis=0)
    Z = X - mu
    cov = Z.T @ Z / max(1, len(X) - 1)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    tot = np.clip(vals, 0, None).sum() or 1.0
    V = vecs[:, :nc].copy()
    for c in range(V.shape[1]):  # sign convention: largest-magnitude loading positive
        if V[np.argmax(np.abs(V[:, c])), c] < 0:
            V[:, c] = -V[:, c]
    return Z @ V, np.clip(vals[:nc], 0, None) / tot


def kmeans(P: np.ndarray, K: int, rnd, reps: int = 3, max_it: int = 100):
    n = len(P)
    best = None
    for _ in range(reps):
        C = np.empty((K, P.shape[1]))
        C[0] = P[int(rnd() * n)]
        D = np.full(n, np.inf)
        for k in range(1, K):
            D = np.minimum(D, ((P - C[k - 1]) ** 2).sum(1))
            cs = np.cumsum(D)
            r = rnd() * cs[-1]
            i = min(int(np.searchsorted(cs, r, side="left")), n - 1)
            C[k] = P[i]
        lab = np.zeros(n, dtype=np.int64)
        sse = 0.0
        for it in range(max_it):
            d = ((P[:, None, :] - C[None, :, :]) ** 2).sum(2)
            new = d.argmin(1)
            sse = float(d[np.arange(n), new].sum())
            changed = it == 0 or bool((new != lab).any())
            lab = new
            if not changed:
                break
            cnt = np.bincount(lab, minlength=K)
            for k in range(K):
                if cnt[k]:
                    C[k] = P[lab == k].sum(0) / cnt[k]
                else:  # re-seed an empty cluster
                    C[k] = P[int(rnd() * n)]
        if best is None or sse < best[2]:
            best = (lab, C.copy(), sse)
    return best


def davies_bouldin(P, lab, C, K):
    s = np.zeros(K)
    for k in range(K):
        m = lab == k
        if not m.any():
            return math.inf
        s[k] = np.sqrt(((P[m] - C[k]) ** 2).sum(1)).mean()
    db = 0.0
    for i in range(K):
        db += max((s[i] + s[j]) / (math.sqrt(((C[i] - C[j]) ** 2).sum()) or 1e-12) for j in range(K) if j != i)
    return db / K


def cluster_best(P, kmin, kmax, rnd):
    best = None
    for K in range(kmin, min(kmax, len(P) - 1) + 1):
        lab, C, _ = kmeans(P, K, rnd)
        db = davies_bouldin(P, lab, C, K)
        if best is None or db < best[1]:
            best = (lab, db, K)
    return (best[0], best[2]) if best else (np.zeros(len(P), dtype=np.int64), 1)


def dip_stat(xs) -> float:
    """Hartigan's dip statistic of a sorted sample (AS 217, as SAMS' HartigansDipTest.m)."""
    N = len(xs)
    if N < 4 or not xs[-1] > xs[0]:
        return 0.0
    x = [0.0] + [float(v) for v in xs]
    min_neg, any_pos, pos_after = math.inf, False, False
    sg = [0] * N
    for i in range(1, N - 1):
        dd = x[i + 2] - 2 * x[i + 1] + x[i]
        sg[i] = -1 if dd > 0 else (1 if dd < 0 else 0)
        if sg[i] < 0 and i < min_neg:
            min_neg = i
    for i in range(1, N - 1):
        if sg[i] > 0:
            any_pos = True
            if i > min_neg:
                pos_after = True
    if not any_pos or min_neg == math.inf or not pos_after:
        return 0.0
    mn = [0] * (N + 1); mj = [0] * (N + 1); gcm = [0] * (N + 2); lcm = [0] * (N + 2)
    fn = float(N)
    low, high, dip = 1, N, 1.0 / fn
    mn[1] = 1
    for j in range(2, N + 1):
        mn[j] = j - 1
        mnj = mn[j]; mnmnj = mn[mnj]; a = mnj - mnmnj; b = j - mnj
        while not (mnj == 1 or (x[j] - x[mnj]) * a < (x[mnj] - x[mnmnj]) * b):
            mn[j] = mnmnj; mnj = mn[j]; mnmnj = mn[mnj]; a = mnj - mnmnj; b = j - mnj
    mj[N] = N
    for jk in range(1, N):
        k = N - jk
        mj[k] = k + 1
        mjk = mj[k]; mjmjk = mj[mjk]; a = mjk - mjmjk; b = k - mjk
        while not (mjk == N or (x[k] - x[mjk]) * a < (x[mjk] - x[mjmjk]) * b):
            mj[k] = mjmjk; mjk = mj[k]; mjmjk = mj[mjk]; a = mjk - mjmjk; b = k - mjk
    for _ in range(4 * N):
        ic = 1; gcm[1] = high; igcm1 = gcm[ic]; ic += 1; gcm[ic] = mn[igcm1]
        while gcm[ic] > low:
            igcm1 = gcm[ic]; ic += 1; gcm[ic] = mn[igcm1]
        icx = ic
        ic = 1; lcm[1] = low; lcm1 = lcm[ic]; ic += 1; lcm[ic] = mj[lcm1]
        while lcm[ic] < high:
            lcm1 = lcm[ic]; ic += 1; lcm[ic] = mj[lcm1]
        icv = ic
        ig, ih, ix, iv, d = icx, icv, icx - 1, 2, 0.0
        if icx == 2 and icv == 2:
            d = 1.0 / fn
        else:
            for _ in range(4 * N):
                igcmx = gcm[ix]; lcmiv = lcm[iv]
                if not igcmx > lcmiv:
                    lcmiv1 = lcm[iv - 1]; a = lcmiv - lcmiv1; b = igcmx - lcmiv1 - 1
                    dx = (x[igcmx] - x[lcmiv1]) * a / (fn * (x[lcmiv] - x[lcmiv1])) - b / fn
                    ix -= 1
                    if not dx < d:
                        d = dx; ig = ix + 1; ih = iv
                else:
                    igcm = gcm[ix]; ig1 = gcm[ix + 1]; a = lcmiv - ig1 + 1; b = igcm - ig1
                    dx = a / fn - ((x[lcmiv] - x[ig1]) * b) / (fn * (x[igcm] - x[ig1]))
                    iv += 1
                    if not dx < d:
                        d = dx; ig = ix + 1; ih = iv - 1
                if ix < 1:
                    ix = 1
                if iv > icv:
                    iv = icv
                if gcm[ix] == lcm[iv]:
                    break
        if d < dip:
            break
        dl = 0.0
        if ig != icx:
            for j in range(ig, icx):
                temp = 1.0 / fn; jb = gcm[j + 1]; je = gcm[j]
                if je - jb > 1 and x[je] != x[jb]:
                    cst = (je - jb) / (fn * (x[je] - x[jb]))
                    for jr in range(jb, je + 1):
                        t = (jr - jb + 1) / fn - (x[jr] - x[jb]) * cst
                        if t > temp:
                            temp = t
                if dl < temp:
                    dl = temp
        du = 0.0
        if ih != icv:
            for k in range(ih, icv):
                temp = 1.0 / fn; kb = lcm[k]; ke = lcm[k + 1]
                if ke - kb > 1 and x[ke] != x[kb]:
                    cst = (ke - kb) / (fn * (x[ke] - x[kb]))
                    for kr in range(kb, ke + 1):
                        t = (x[kr] - x[kb]) * cst - (kr - kb - 1) / fn
                        if t > temp:
                            temp = t
                if du < temp:
                    du = temp
        dip = max(dip, dl, du)
        low = gcm[ig]; high = lcm[ih]
    return 0.5 * dip


_DIP_NULL: dict = {}   # sample size -> [dips of uniform samples, rng state]


def dip_multimodal(vals, rnd, nboot: int = 200, cap: int = 2000, alpha: float = 0.05) -> bool:
    """True when the dip test rejects unimodality (p < alpha). Same subsampling and
    bootstrap null as the viewer; the null is only extended as far as the decision needs."""
    v = np.asarray(vals, dtype=np.float64).copy()
    if len(v) > cap:
        for i in range(cap):
            j = i + int(rnd() * (len(v) - i))
            v[i], v[j] = v[j], v[i]
        v = v[:cap]
    v.sort()
    n = len(v)
    dip = dip_stat(v)
    if n not in _DIP_NULL:
        _DIP_NULL[n] = [[], Mulberry32(n * 7919)]
    nul, r = _DIP_NULL[n]
    limit = int(math.ceil(alpha * nboot))      # p < alpha  <=>  fewer than `limit` null dips >= dip
    count = 0
    for b in range(nboot):
        if b == len(nul):
            u = np.array([r() for _ in range(n)])
            u.sort()
            nul.append(dip_stat(u))
        if nul[b] >= dip:
            count += 1
            if count >= limit:
                return False
    return True


def rmse_outliers(W, g, sd_cut):
    if len(g) < 3:
        return g, 0
    A = W[g].astype(np.float64)
    r = np.sqrt(((A - A.mean(0)) ** 2).mean(1))
    mu, sd = r.mean(), r.std(ddof=1)
    keep = g[np.abs(r - mu) <= sd_cut * sd]
    return keep, len(g) - len(keep)


def dtw(a, b) -> float:
    m = len(b)
    prev = [math.inf] * (m + 1); prev[0] = 0.0
    for i in range(1, len(a) + 1):
        cur = [math.inf] * (m + 1)
        ai = a[i - 1]
        for j in range(1, m + 1):
            cur[j] = abs(ai - b[j - 1]) + min(prev[j], cur[j - 1], prev[j - 1])
        prev = cur
    return prev[m]


def merge_plan(T, prm):
    lo = min(float(t.min()) for t in T); hi = max(float(t.max()) for t in T)
    N = [((t - lo) / ((hi - lo) or 1)).tolist() for t in T]
    tol = 1 - prm["overlap"]
    parent = list(range(len(T)))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    close_not_merged = False
    for i in range(len(T)):
        for j in range(i + 1, len(T)):
            d = dtw(N[i], N[j])
            pk1, tr1, pk2, tr2 = max(N[i]), min(N[i]), max(N[j]), min(N[j])
            h1 = pk1 - tr1
            same = abs(pk1 - pk2) < h1 * tol and abs(tr1 - tr2) < h1 * tol
            if d < prm["merge_thr"] and same:
                parent[root(j)] = root(i)
            elif d < prm["merge_thr"]:
                close_not_merged = True
    comps: dict = {}
    for i in range(len(T)):
        comps.setdefault(root(i), []).append(i)
    return list(comps.values()), close_not_merged


# --------------------------------------------------------------------------- results
@dataclass
class Unit:
    unit: int
    spikes: np.ndarray          # indices into the electrode's (time-ordered) spikes
    n: int
    rate_hz: float
    mean_uV: np.ndarray
    sd_uV: np.ndarray
    peak_to_peak_uV: float
    peak_uV: float
    noise_uV: float
    snr: float
    isi_cv: float
    min_isi_s: float
    refractory_removed: int


@dataclass
class ElectrodeSort:
    channel: int
    label: str
    n: int
    spike_index: np.ndarray     # global spike indices, time-ordered
    labels: np.ndarray          # unit per spike_index (0 = unsorted)
    units: list = field(default_factory=list)
    skipped: str = ""
    initial_k: int = 0
    dip_splits: int = 0
    merges: int = 0
    dropped_slow: int = 0
    possible_mua: bool = False
    pcs: np.ndarray | None = None
    explained: np.ndarray | None = None
    seconds: float = 0.0


def _unit_stats(k, g, W, times, noise_v, duration, viol):
    A = W[g].astype(np.float64)
    mean = A.mean(0)
    sd = np.sqrt(np.clip((A * A).mean(0) - mean ** 2, 0, None))
    mx, mn = float(mean.max()), float(mean.min())
    isi = np.diff(times[g])
    if len(isi) > 1:
        mi = isi.mean()
        cv = math.sqrt(max(0.0, (float((isi ** 2).sum()) - len(isi) * mi * mi) / (len(isi) - 1))) / mi
    else:
        cv = math.nan
    noise = float(noise_v[g].mean() * 1e6)
    return Unit(k, g, len(g), len(g) / duration, mean, sd, mx - mn, mx if abs(mx) > abs(mn) else mn, noise,
                (mx - mn) / noise, cv, float(isi.min()) if len(isi) else math.nan, viol)


def sort_electrode(W, times, noise_v, duration, channel=0, label="", spike_index=None, **params) -> ElectrodeSort:
    """Sort one electrode.
    W: (n, samples) waveforms in uV (float32), time-ordered; times: spike times (s);
    noise_v: detector noise SD per spike (V); duration: recording length (s);
    channel: channel-array index (seeds the random generator, as in the viewer)."""
    prm = {**DEFAULTS, **params}
    t0 = time.perf_counter()
    n = len(W)
    res = ElectrodeSort(channel, label, n, np.arange(n) if spike_index is None else spike_index, np.zeros(n, dtype=np.int16))
    min_n = max(20, math.ceil(duration * prm["min_rate"]))
    if n < min_n:
        res.skipped = f"{n} spikes; at least {min_n} needed"
        return res
    scores, explained = pca_scores(W.astype(np.float64), 3)
    res.pcs, res.explained = scores.astype(np.float32), explained
    rnd = Mulberry32(0x5A15 + channel)

    def split(g):
        lab, K = cluster_best(scores[g, :2], 2, prm["max_k"], rnd)
        return [g[lab == k] for k in range(K) if (lab == k).any()]

    groups = split(np.arange(n))
    res.initial_k = len(groups)
    groups = [rmse_outliers(W, g, prm["sd_cut"])[0] for g in groups]
    nxt = []
    for g in groups:
        if len(g) <= 10:
            continue
        m1 = dip_multimodal(scores[g, 0], rnd)
        m2 = dip_multimodal(scores[g, 1], rnd)   # always run both: each draws from the generator, as in the viewer
        if m1 or m2:
            sub = split(g)
            if len(sub) > 1:
                res.dip_splits += 1
                nxt.extend(sub)
                continue
        nxt.append(g)
    groups = nxt
    for _ in range(10):
        if len(groups) < 2:
            break
        comps, _c = merge_plan([W[g].astype(np.float64).mean(0) for g in groups], prm)
        if len(comps) == len(groups):
            break
        res.merges += len(groups) - len(comps)
        groups = [np.concatenate([groups[i] for i in comp]) for comp in comps]
    refr = prm["refractory_ms"] / 1000
    kept = []
    for g in groups:
        g = np.sort(g)
        for _ in range(W.shape[1]):
            if not len(g):
                break
            g, out = rmse_outliers(W, g, prm["sd_cut"])
            if not out:
                break
        viol = 0
        if len(g) > 1:
            drop = np.concatenate([[False], np.diff(times[g]) < refr])
            viol = int(drop.sum())
            g = g[~drop]
        if not len(g) or len(g) / duration <= prm["min_rate"]:
            res.dropped_slow += 1
            continue
        kept.append((g, viol))
    kept.sort(key=lambda u: -len(u[0]))
    for k, (g, viol) in enumerate(kept, 1):
        res.units.append(_unit_stats(k, g, W, times, noise_v, duration, viol))
        res.labels[g] = k
    if len(res.units) > 1:
        res.possible_mua = merge_plan([u.mean_uV for u in res.units], prm)[1]
    res.seconds = time.perf_counter() - t0
    return res


@dataclass
class SortResult:
    file: str
    params: dict
    unit: np.ndarray                    # per spike: -1 not sorted, 0 unsorted, 1.. unit
    electrodes: dict                    # label -> ElectrodeSort
    treatments: dict = field(default_factory=dict)

    def units_table(self):
        rows = []
        for lab, e in self.electrodes.items():
            uns = e.n - sum(u.n for u in e.units)
            for u in e.units:
                rows.append(dict(electrode=lab, well=lab.split("_")[0], treatment=self.treatments.get(lab.split("_")[0], ""),
                                 unit=u.unit, spikes=u.n, rate_hz=u.rate_hz, peak_to_peak_uV=u.peak_to_peak_uV,
                                 peak_uV=u.peak_uV, snr=u.snr, isi_cv=u.isi_cv, min_isi_ms=u.min_isi_s * 1000,
                                 refractory_removed=u.refractory_removed, unsorted_on_electrode=uns))
        return rows

    def write_units_csv(self, path):
        """Same columns and formatting as the viewer's 'Export units (CSV)'."""
        f = lambda v, d: "" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"
        with open(path, "w") as fh:
            fh.write("electrode,well,treatment,unit,spikes,rate_hz,peak_to_peak_uV,peak_uV,snr,isi_cv,min_isi_ms,refractory_removed,unsorted_on_electrode\n")
            for r in self.units_table():
                tr = '"' + r["treatment"].replace('"', '""') + '"' if r["treatment"] else ""
                fh.write(",".join([r["electrode"], r["well"], tr, str(r["unit"]), str(r["spikes"]), f(r["rate_hz"], 4),
                                   f(r["peak_to_peak_uV"], 2), f(r["peak_uV"], 2), f(r["snr"], 2), f(r["isi_cv"], 3),
                                   f(r["min_isi_ms"], 3), str(r["refractory_removed"]), str(r["unsorted_on_electrode"])]) + "\n")


def sort_spikes(sf, electrodes=None, progress=None, **params) -> SortResult:
    """Sort every electrode with spikes (or only the given labels, e.g. ["A6_11"])."""
    prm = {**DEFAULTS, **params}
    s = sf.spikes
    times = sf.spike_times_s
    lut = sf._channel_lut()
    ch = lut[s["chip"].astype(np.int64) * 256 + s["channel"].astype(np.int64)]
    if len(times) and np.any(np.diff(times) < 0):
        order = np.argsort(times, kind="stable")
    else:
        order = np.arange(len(times))
    duration = sf.header.get("Duration") or (float(times.max()) if len(times) else 1.0)
    scale = float(sf.voltage_scale) * 1e6
    labels = [sf._label(c) for c in sf.channels]
    unit = np.full(len(times), -1, dtype=np.int16)
    wanted = None if electrodes is None else set(electrodes)
    out = {}
    chans = [c for c in np.unique(ch[ch >= 0]) if wanted is None or labels[c] in wanted]
    chans.sort(key=lambda c: tuple(int(sf.channels[c][k]) for k in ("WellRow", "WellColumn", "ElectrodeRow", "ElectrodeColumn")))
    for k, c in enumerate(chans):
        idx = order[ch[order] == c]
        idx = idx[np.argsort(times[idx], kind="stable")]
        W = (s["data"][idx].astype(np.float64) * scale).astype(np.float32)
        if progress:
            progress(k, len(chans), labels[c], len(idx))
        e = sort_electrode(W, times[idx], s["standardDeviation"][idx], duration, channel=int(c), label=labels[c], spike_index=idx, **prm)
        unit[idx] = e.labels
        out[labels[c]] = e
    treat = {}
    for t in sf.tags:
        if "TreatmentWhat" in t and t.get("WellRow"):
            treat[f"{chr(64 + t['WellRow'])}{t['WellColumn']}"] = t["TreatmentWhat"]
    return SortResult(sf.path, prm, unit, out, treat)


def _main(argv=None):
    from axion_spk import read_spk
    ap = argparse.ArgumentParser(description="SAMS-style automatic spike sorting of an Axion .spk file")
    ap.add_argument("spk")
    ap.add_argument("--electrodes", nargs="+", help="only these electrodes, e.g. A6_11 B4_12")
    ap.add_argument("--units", help="write one row per unit to this CSV")
    ap.add_argument("--csv", help="write the spike table with a unit column to this CSV")
    ap.add_argument("--npz", help="write times, electrode, unit and unit templates to this .npz")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=type(v), default=v, help=f"default {v}")
    a = ap.parse_args(argv)
    sf = read_spk(a.spk)
    prm = {k: getattr(a, k) for k in DEFAULTS}
    t0 = time.perf_counter()
    prog = lambda k, n, lab, ns: print(f"\r[{k + 1}/{n}] {lab} ({ns} spikes)      ", end="", file=sys.stderr, flush=True)
    res = sort_spikes(sf, electrodes=a.electrodes, progress=prog, **prm)
    print(file=sys.stderr)
    nu = sum(len(e.units) for e in res.electrodes.values())
    print(f"{len(res.electrodes)} electrodes, {nu} units in {time.perf_counter() - t0:.1f} s")
    for lab, e in res.electrodes.items():
        if e.skipped:
            print(f"  {lab}: not sorted ({e.skipped})")
        else:
            uns = e.n - sum(u.n for u in e.units)
            print(f"  {lab}: {e.n} spikes -> {len(e.units)} units {[u.n for u in e.units]}, {uns} unsorted"
                  + ("  [possible multi-unit]" if e.possible_mua else ""))
    if a.units:
        res.write_units_csv(a.units); print(f"wrote {a.units}")
    if a.csv or a.npz:
        labels = sf.electrode_labels()
    if a.csv:
        wf = sf.waveforms_uv
        trig = sf.spikes["triggerSample"]
        amp = wf[np.arange(len(wf)), np.clip(trig, 0, wf.shape[1] - 1)]
        with open(a.csv, "w") as fh:
            fh.write("time_s,electrode,trigger_amplitude_uV,min_uV,max_uV,std_uV,threshold_multiplier,unit\n")
            for t, lab, am, mn, mx, sd, th, u in zip(sf.spike_times_s, labels, amp, wf.min(1), wf.max(1),
                                                     sf.spikes["standardDeviation"] * 1e6, sf.spikes["thresholdMultiplier"], res.unit):
                fh.write(f"{t:.6f},{lab},{am:.3f},{mn:.3f},{mx:.3f},{sd:.3f},{th:.3f},{'' if u < 0 else u}\n")
        print(f"wrote {a.csv}")
    if a.npz:
        tl, tu, tm = [], [], []
        for lab, e in res.electrodes.items():
            for u in e.units:
                tl.append(lab); tu.append(u.unit); tm.append(u.mean_uV)
        np.savez_compressed(a.npz, times_s=sf.spike_times_s, electrode=labels, unit=res.unit,
                            template_electrode=np.array(tl), template_unit=np.array(tu),
                            templates_uV=np.array(tm) if tm else np.zeros((0, 0)))
        print(f"wrote {a.npz}")


if __name__ == "__main__":
    _main()
