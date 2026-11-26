"use strict";
/* =====================================================================
 * MEA pipeline: analysis functions
 *
 * Everything that computes a result from the recording is here; the page
 * (js/mea_pipeline.js) only reads the controls, calls these functions and
 * draws what they return. None of them reads page state or touches the DOM:
 * the recording R, spike times and parameters are passed in. Times are in s,
 * voltages in µV. The definitions are written out in the page's info boxes.
 *
 * Index
 *   General       median, quantileSorted, sampleSD, lowerBound, upperBound
 *   Recording     parseSettings, prepare
 *   Replay        windowCounts
 *   Waveforms     readSpike, alignIndex, waveStats, waveMeasures
 *   Firing stats  binCounts, electrodeStats, isiHistogram, autocorrelogram   (STATS block)
 *   Net. bursts   NB_METHODS, nbRequired, nbPooled, nbIsinValues, nbIsin, nbSingle, nbGroup,
 *                 nbRecruitment, nbRecruit, nbMerge, nbLower, nbUpper, nbDescribe, nbSummary,
 *                 nbAgreement, nbIsinHist, nbValley, nbWell                   (NB block)
 *   Comparison    nbCompare, firingCompare
 *
 * The STATS and NB blocks are also run in Node by verification/dump_page_stats.js
 * and verification/dump_page_nb.js, and checked against python/firing_stats.py
 * and python/network_bursts.py.
 * ===================================================================== */

/* ---------------------------------------------------------------------
 * General
 * --------------------------------------------------------------------- */
function median(a) { if (!a.length) return NaN; const b = Float64Array.from(a).sort(), m = b.length >> 1; return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2; }
function lowerBound(arr, v) { let lo = 0, hi = arr.length; while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < v) lo = m + 1; else hi = m; } return lo; }
function upperBound(arr, v) { let lo = 0, hi = arr.length; while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] <= v) lo = m + 1; else hi = m; } return lo; }

/* ---------------------------------------------------------------------
 * Recording: per-electrode data derived from the parsed file (R = parseSpk(...))
 * --------------------------------------------------------------------- */
function parseSettings(desc) { // AxIS "Description": CSV-like sections of key,value lines
  const sections = []; let cur = { title: "General", rows: [] };
  for (const raw of desc.split(/\r?\n/)) {
    const line = raw.trim(); if (!line) continue;
    const i = line.indexOf(","), k = i < 0 ? line : line.slice(0, i).trim(), v = i < 0 ? "" : line.slice(i + 1).trim();
    if (!k && !v) continue;
    if (k && !v) { if (cur.rows.length) sections.push(cur); cur = { title: k, rows: [] }; continue; }
    cur.rows.push([k, v]);
  }
  if (cur.rows.length) sections.push(cur);
  return sections;
}
function prepare(R) {
  const h = R.header, nC = R.channels.length, ns = h.NumSamplesPerBlock, uV = h.VoltageScale * 1e6;
  R.duration = h.Duration || (R.n ? R.times[R.n - 1] : 1);
  // spikes of each electrode, in time order
  const counts = new Int32Array(nC);
  for (let i = 0; i < R.n; i++) if (R.ch[i] >= 0) counts[R.ch[i]]++;
  R.perCh = Array.from(counts, c => new Int32Array(c));
  const fill = new Int32Array(nC);
  for (let i = 0; i < R.n; i++) { const c = R.ch[i]; if (c >= 0) R.perCh[c][fill[c]++] = i; }
  R.counts = counts; R.rate = Array.from(counts, c => c / R.duration);
  // recorded electrodes: the channel list of the spike header (else those with spikes)
  const byId = new Map(R.channels.map(c => [c.achk * 256 + c.idx, c.i]));
  R.recorded = new Set(R.channelIds.map(([a, c]) => byId.get(a * 256 + c)).filter(x => x !== undefined));
  for (let c = 0; c < nC; c++) if (counts[c] > 0) R.recorded.add(c);
  // peak-to-peak amplitude of every spike (µV) and per-electrode medians
  R.ptp = new Float32Array(R.n);
  for (let i = 0; i < R.n; i++) {
    const o = R.waveOffset(i); let mn = Infinity, mx = -Infinity;
    for (let k = 0; k < ns; k++) { const v = R.dv.getInt16(o + 2 * k, true) * uV; if (v < mn) mn = v; if (v > mx) mx = v; }
    R.ptp[i] = mx - mn;
  }
  R.medPtp = R.perCh.map(idx => median(Array.from(idx, i => R.ptp[i])));
  R.medNoise = R.perCh.map(idx => median(Array.from(idx, i => R.sd[i] * 1e6)));
  R.thrMult = R.perCh.map(idx => idx.length ? R.thr[idx[0]] : NaN);
  // plate geometry from the channel map
  R.geo = { WR: 0, WC: 0, ER: 0, EC: 0 };
  for (const c of R.channels) { R.geo.WR = Math.max(R.geo.WR, c.wr); R.geo.WC = Math.max(R.geo.WC, c.wc); R.geo.ER = Math.max(R.geo.ER, c.er); R.geo.EC = Math.max(R.geo.EC, c.ec); }
  R.wellNames = []; for (let r = 1; r <= R.geo.WR; r++) for (let c = 1; c <= R.geo.WC; c++) R.wellNames.push(String.fromCharCode(64 + r) + c);
  const key = c => [c.wr, c.wc, c.er, c.ec];
  R.rows = [...R.recorded].sort((a, b) => { const A = key(R.channels[a]), B = key(R.channels[b]); for (let k = 0; k < 4; k++) if (A[k] !== B[k]) return A[k] - B[k]; return 0; });
  R.recordedWells = new Set([...R.recorded].map(c => R.channels[c].well));
  R.settings = parseSettings(h.Description || "");
  R.plateName = (() => { for (const s of R.settings) for (const [k, v] of s.rows) if (k === "Plate Type" && v) return v;
    return (PLATES[R.plateType] || `0x${(R.plateType >>> 0).toString(16)}`).replace(/([a-z])([A-Z])/g, "$1 $2"); })();
  // shortest interval between two spikes of the same electrode, over the whole recording: below it, nothing can be measured
  R.minIsiAll = Infinity;
  for (const idx of R.perCh) for (let k = 1; k < idx.length; k++) { const d = R.times[idx[k]] - R.times[idx[k - 1]]; if (d < R.minIsiAll) R.minIsiAll = d; }
  R.overlap = (() => { for (const s of R.settings) for (const [k, v] of s.rows) if (k === "Allow Sample Overlap") return v; return null; })();
  R.preTrigMs = R.n ? median(Array.from(R.trig.subarray(0, Math.min(R.n, 1000)))) / R.fs * 1000 : NaN;
}

/* ---------------------------------------------------------------------
 * Activity replay: spikes per electrode in the causal window (t − win, t]
 * --------------------------------------------------------------------- */
function windowCounts(R, t, win, out) { // spikes per electrode in (t - win, t]
  out.fill(0);
  const a = upperBound(R.times, t - win), b = upperBound(R.times, t);
  for (let i = a; i < b; i++) { const c = R.ch[i]; if (c >= 0) out[c]++; }
  return b - a;
}

/* ---------------------------------------------------------------------
 * Waveforms: samples exactly as stored, aligned by whole samples
 * --------------------------------------------------------------------- */
function readSpike(R, i) { const ns = R.header.NumSamplesPerBlock, o = R.waveOffset(i), uV = R.header.VoltageScale * 1e6, a = new Float64Array(ns); for (let k = 0; k < ns; k++) a[k] = R.dv.getInt16(o + 2 * k, true) * uV; return a; }
function alignIndex(R, w, i, mode) { // sample used as time zero for this spike
  if (mode === "trig") return Math.min(Math.max(R.trig[i], 0), w.length - 1);
  let best = 0; for (let k = 1; k < w.length; k++) if (mode === "neg" ? w[k] < w[best] : mode === "pos" ? w[k] > w[best] : Math.abs(w[k]) > Math.abs(w[best])) best = k;
  return best;
}
// waveforms + per-lag statistics of a set of spikes (lag = sample index - alignment index)
function waveStats(R, idx, mode) {
  const ns = R.header.NumSamplesPerBlock, waves = [], al = [];
  for (const i of idx) { const w = readSpike(R, i); waves.push(w); al.push(alignIndex(R, w, i, mode)); }
  if (!waves.length) return null;
  const lagMin = -Math.max(...al), lagMax = ns - 1 - Math.min(...al), L = lagMax - lagMin + 1;
  const sum = new Float64Array(L), sq = new Float64Array(L), cnt = new Int32Array(L);
  waves.forEach((w, q) => { for (let k = 0; k < ns; k++) { const j = k - al[q] - lagMin; sum[j] += w[k]; sq[j] += w[k] * w[k]; cnt[j]++; } });
  const mean = new Float64Array(L), sd = new Float64Array(L), ok = new Uint8Array(L);
  for (let j = 0; j < L; j++) { if (!cnt[j]) continue; mean[j] = sum[j] / cnt[j]; sd[j] = Math.sqrt(Math.max(0, sq[j] / cnt[j] - mean[j] ** 2)); ok[j] = cnt[j] >= waves.length / 2; }
  return { waves, al, lagMin, lagMax, mean, sd, ok, n: waves.length };
}
function waveMeasures(R, st, c) { // measures of the mean waveform (only lags present in >= half of the spikes)
  let tr = -1, pk = -1;
  for (let j = 0; j < st.mean.length; j++) { if (!st.ok[j]) continue; if (tr < 0 || st.mean[j] < st.mean[tr]) tr = j; if (pk < 0 || st.mean[j] > st.mean[pk]) pk = j; }
  let after = -1; for (let j = tr + 1; j < st.mean.length; j++) if (st.ok[j] && (after < 0 || st.mean[j] > st.mean[after])) after = j; // repolarisation peak after the trough
  const noise = R.medNoise[c];
  return { trough: st.mean[tr], peak: st.mean[pk], ptp: st.mean[pk] - st.mean[tr], troughToPeakMs: after > 0 ? (after - tr) / R.fs * 1000 : NaN, snr: (st.mean[pk] - st.mean[tr]) / noise, noise };
}

/* ---------------------------------------------------------------------
 * Firing statistics
 * --------------------------------------------------------------------- */
/* STATS-BEGIN: firing statistics, pure functions (also run in Node by verification/compare_firing_stats.py)
 * t: spike times (s) of one electrode inside the analysis period [a, a + T], ascending. */
function quantileSorted(s, q) { // linear interpolation between order statistics (same as numpy's default)
  if (!s.length) return NaN; const h = (s.length - 1) * q, i = Math.floor(h); return i + 1 < s.length ? s[i] + (h - i) * (s[i + 1] - s[i]) : s[i];
}
function sampleSD(x, mean) { if (x.length < 2) return NaN; let v = 0; for (const a of x) v += (a - mean) ** 2; return Math.sqrt(v / (x.length - 1)); }
function binCounts(t, a, T, bin) { // counts in K = floor(T / bin) whole bins [a + k·bin, a + (k+1)·bin)
  const K = Math.floor(T / bin + 1e-9), n = new Float64Array(K);
  for (const x of t) { const k = Math.floor((x - a) / bin); if (k >= 0 && k < K) n[k]++; }
  return n;
}
function electrodeStats(t, a, T, prm) {
  const N = t.length, rate = N / T, out = { N, rate, active: rate * 60 >= prm.activePerMin };
  // Fano factor of counts in bins of prm.fanoBin seconds (sample variance / mean)
  const n = binCounts(t, a, T, prm.fanoBin); let m = 0; for (const v of n) m += v; m /= n.length || 1;
  out.fano = n.length >= 2 && m > 0 ? sampleSD(n, m) ** 2 / m : NaN;
  // inter-spike intervals
  const tau = new Float64Array(Math.max(0, N - 1)); for (let i = 0; i + 1 < N; i++) tau[i] = t[i + 1] - t[i];
  out.nIsi = tau.length;
  if (tau.length >= prm.minIsi) {
    let s = 0; for (const x of tau) s += x; const mu = s / tau.length;
    const srt = Float64Array.from(tau).sort();
    out.meanIsi = mu; out.medIsi = quantileSorted(srt, .5); out.cv = sampleSD(tau, mu) / mu;
    let c2 = 0, lv = 0, nPair = 0;
    for (let i = 0; i + 1 < tau.length; i++) { const sum = tau[i] + tau[i + 1]; if (sum <= 0) continue; c2 += 2 * Math.abs(tau[i + 1] - tau[i]) / sum; lv += ((tau[i] - tau[i + 1]) / sum) ** 2; nPair++; }
    out.cv2 = nPair ? c2 / nPair : NaN; out.lv = nPair ? 3 * lv / nPair : NaN;
    let sh = 0; for (const x of tau) if (x < prm.shortIsi) sh++; out.shortPct = 100 * sh / tau.length;
  } else { out.meanIsi = out.medIsi = out.cv = out.cv2 = out.lv = out.shortPct = NaN; }
  return out;
}
function isiHistogram(t, perDecade, lo, hi) { // fraction of ISIs per log-spaced bin; edges lo·10^(k/perDecade)
  const nb = Math.round(Math.log10(hi / lo) * perDecade), h = new Float64Array(nb); let n = 0, out = 0;
  for (let i = 0; i + 1 < t.length; i++) { const d = t[i + 1] - t[i]; n++; const k = Math.floor(Math.log10(d / lo) * perDecade + 1e-9); if (d > 0 && k >= 0 && k < nb) h[k]++; else out++; }
  if (n) for (let k = 0; k < nb; k++) h[k] /= n;
  return { h, nb, n, outside: out, edge: k => lo * 10 ** (k / perDecade) };
}
function autocorrelogram(t, W, b) { // conditional rate (spikes/s) at lag τ given a spike at 0, bins [k·b, (k+1)·b) for k = -K..K-1, self pairs excluded
  const K = Math.round(W / b), c = new Float64Array(2 * K), N = t.length;
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const d = t[j] - t[i]; if (d >= W) break;
    const k = Math.floor(d / b + 1e-9); if (k < K) { c[K + k]++; c[K - 1 - k]++; }   // lag +d in bin k, lag -d mirrored
  }
  for (let k = 0; k < 2 * K; k++) c[k] = N ? c[k] / (N * b) : 0;
  return { c, K, b };
}
/* STATS-END */

/* ---------------------------------------------------------------------
 * Network bursts
 * --------------------------------------------------------------------- */
/* NB-BEGIN: network bursts, pure functions (also run in Node by verification/compare_network_bursts.py)
 * trains: [{ c, t }] for the well's active electrodes in ascending electrode index; t: ascending spike times (s) in the period [a, a + T].
 * A network burst is [start, end] (s). Same definitions as network_bursts.py. */
const NB_METHODS = ["isin", "group", "recruit"];
function nbRequired(pct, A, m) { return Math.max(m, Math.ceil(pct / 100 * A - 1e-9)); }
function nbPooled(trains) { // pooled spike train sorted by time, ties by electrode index
  let M = 0; for (const x of trains) M += x.t.length;
  const t = new Float64Array(M), e = new Int32Array(M); let k = 0;
  for (const x of trains) for (const v of x.t) { t[k] = v; e[k++] = x.c; }
  const o = new Uint32Array(M); for (let i = 0; i < M; i++) o[i] = i;
  o.sort((p, q) => t[p] - t[q] || e[p] - e[q]);
  return { s: Float64Array.from(o, i => t[i]), e: Int32Array.from(o, i => e[i]) };
}
function nbIsinValues(s, N) { const n = Math.max(0, s.length - N + 1), d = new Float64Array(n); for (let i = 0; i < n; i++) d[i] = s[i + N - 1] - s[i]; return d; }
function nbIsin(s, e, N, thr, minSpikes, m) { // windows of N spikes spanning <= thr, merged when they share a spike
  const d = nbIsinValues(s, N), out = []; let i0 = -1, iend = -1;
  const close = () => { if (i0 < 0 || iend - i0 + 1 < minSpikes) return; const u = new Set(); for (let i = i0; i <= iend; i++) u.add(e[i]); if (u.size >= m) out.push([s[i0], s[iend]]); };
  for (let i = 0; i < d.length; i++) { if (!(d[i] <= thr)) continue; if (i0 >= 0 && i <= iend) iend = i + N - 1; else { close(); i0 = i; iend = i + N - 1; } }
  close(); return out;
}
function nbSingle(t, isiMax, nMin) { // maximal runs of >= nMin spikes of one electrode with every ISI <= isiMax
  const out = [], n = t.length; let j = 0;
  for (let i = 1; i <= n; i++) if (i === n || t[i] - t[i - 1] > isiMax) { if (i - j >= nMin) out.push([t[j], t[i - 1]]); j = i; }
  return out;
}
function nbGroup(trains, isiMax, nMin, gap, req) { // burst onsets <= gap apart grouped; kept with >= req electrodes
  const b = []; for (const x of trains) for (const [on, end] of nbSingle(x.t, isiMax, nMin)) b.push([on, x.c, end]);
  b.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
  const out = []; let g = [];
  const close = () => { if (!g.length || new Set(g.map(x => x[1])).size < req) return; let en = -Infinity; for (const x of g) if (x[2] > en) en = x[2]; out.push([g[0][0], en]); };
  for (const x of b) { if (g.length && x[0] - g[g.length - 1][0] > gap) { close(); g = []; } g.push(x); }
  close(); return out;
}
function nbRecruitment(trains, a, T, b) { // electrodes with >= 1 spike in each whole bin [a + k·b, a + (k+1)·b)
  const K = Math.max(0, Math.floor(T / b + 1e-9)), rec = new Int32Array(K);
  for (const x of trains) { let last = -1; for (const v of x.t) { const k = Math.floor((v - a) / b); if (k >= 0 && k < K && k !== last) { rec[k]++; last = k; } } }
  return rec;
}
function nbRecruit(s, a, T, b, trains, req) { // runs of bins with >= req electrodes firing, from their first to their last spike
  const rec = nbRecruitment(trains, a, T, b), K = rec.length, out = [], kOf = j => Math.floor((s[j] - a) / b); let i = 0;
  for (let r = 0; r < K;) {
    if (rec[r] < req) { r++; continue; }
    let r1 = r; while (r1 + 1 < K && rec[r1 + 1] >= req) r1++;
    while (i < s.length && kOf(i) < r) i++;
    let j = i; while (j < s.length && kOf(j) <= r1) j++;
    if (j > i) out.push([s[i], s[j - 1]]);
    i = j; r = r1 + 1;
  }
  return out;
}
function nbMerge(nb) { // network bursts whose intervals overlap become one
  const out = []; for (const [st, en] of nb.slice().sort((p, q) => p[0] - q[0] || p[1] - q[1])) {
    const l = out[out.length - 1]; if (l && st <= l[1]) l[1] = Math.max(l[1], en); else out.push([st, en]); }
  return out;
}
function nbLower(a, v) { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m; } return lo; }
function nbUpper(a, v) { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] <= v) lo = m + 1; else hi = m; } return lo; }
function nbDescribe(nb, s, e) {
  return nb.map(([st, en]) => { const i0 = nbLower(s, st), i1 = nbUpper(s, en), u = new Set(); for (let i = i0; i < i1; i++) u.add(e[i]);
    return { start: st, end: en, duration: en - st, spikes: i1 - i0, electrodes: u.size }; });
}
function nbSummary(rows, M, T) {
  const n = rows.length, med = a => a.length ? quantileSorted(Float64Array.from(a).sort(), .5) : NaN, mean = a => a.reduce((x, y) => x + y, 0) / a.length;
  const ibi = []; for (let i = 1; i < n; i++) ibi.push(rows[i].start - rows[i - 1].start);
  const mi = ibi.length ? mean(ibi) : NaN; let sp = 0; for (const r of rows) sp += r.spikes;
  return { n, perMin: 60 * n / T, durMedian: med(rows.map(r => r.duration)), durMean: n ? mean(rows.map(r => r.duration)) : NaN, ibiMean: mi,
    ibiCv: ibi.length > 1 ? sampleSD(ibi, mi) / mi : NaN, spikesMedian: med(rows.map(r => r.spikes)), electrodesMedian: med(rows.map(r => r.electrodes)), pctSpikes: M ? 100 * sp / M : NaN };
}
function nbAgreement(x, y) { // % of the network bursts in x that overlap at least one in y
  if (!x.length) return NaN; let hit = 0;
  for (const [st, en] of x) for (const [ys, ye] of y) if (ys <= en && ye >= st) { hit++; break; }
  return 100 * hit / x.length;
}
function nbIsinHist(d, perDecade, lo, hi) { // counts per log-spaced bin, edges lo·10^(k/perDecade)
  const nb = Math.round(Math.log10(hi / lo) * perDecade), h = new Float64Array(nb);
  for (const v of d) { const k = Math.floor(Math.log10(v / lo) * perDecade + 1e-9); if (v > 0 && k >= 0 && k < nb) h[k]++; }
  return { h, nb, n: d.length, edge: k => lo * 10 ** (k / perDecade) };
}
function nbValley(h, lo, perDecade) { // deepest bin after the fullest one; threshold = its upper edge
  const nb = h.length; let p = 0; for (let k = 1; k < nb; k++) if (h[k] > h[p]) p = k;
  if (!nb || !h[p]) return NaN;
  const suf = new Float64Array(nb + 1); for (let k = nb - 1; k >= 0; k--) suf[k] = Math.max(suf[k + 1], h[k]);
  let best = 0, bv = -1; for (let v = p + 1; v < nb - 1; v++) { const depth = Math.min(h[p], suf[v + 1]) - h[v]; if (depth > best) { best = depth; bv = v; } }
  return bv >= 0 ? lo * 10 ** ((bv + 1) / perDecade) : NaN;
}
function nbWell(trains, a, T, p) {
  const A = trains.length, { s, e } = nbPooled(trains), nb = {
    isin: nbIsin(s, e, p.isinN, p.isinThr, p.isinMinSpikes, p.minEl),
    group: nbGroup(trains, p.burstIsi, p.burstMinSpikes, p.onsetGap, nbRequired(p.groupPct, A, p.minEl)),
    recruit: nbRecruit(s, a, T, p.recBin, trains, nbRequired(p.recPct, A, p.minEl)) };
  const bursts = {}, summary = {}, agreement = {};
  for (const k of NB_METHODS) { nb[k] = nbMerge(nb[k]); bursts[k] = nbDescribe(nb[k], s, e); summary[k] = nbSummary(bursts[k], s.length, T); }
  for (const x of NB_METHODS) { agreement[x] = {}; for (const y of NB_METHODS) agreement[x][y] = nbAgreement(nb[x], nb[y]); }
  return { active: A, spikes: s.length, s, e, bursts, summary, agreement };
}
/* NB-END */

/* ---------------------------------------------------------------------
 * Comparison between recordings
 * --------------------------------------------------------------------- */
function nbCompare(res, T, well) { // {method: {n, rate, dur}} in one well, or over all wells ("*": counts summed, durations pooled)
  const o = {};
  for (const m of NB_METHODS) {
    const rs = well === "*" ? [...res.values()] : res.has(well) ? [res.get(well)] : [];
    if (!rs.length) { o[m] = { n: NaN, rate: NaN, dur: NaN }; continue; }
    const d = rs.flatMap(r => r.bursts[m].map(b => b.duration)), n = d.length;
    o[m] = { n, rate: 60 * n / T, dur: n ? quantileSorted(Float64Array.from(d).sort(), .5) * 1000 : NaN };
  }
  return o;
}

function firingCompare(R, el, T, well) { // firing of one well, or of all recorded wells ("*"), from the per-electrode statistics el
  // net: all spikes of the recorded electrodes / T (spikes/s); mean, sd: firing rate over the active electrodes (Hz); active, recorded: counts
  const cs = R.rows.filter(c => well === "*" || R.channels[c].well === well);
  if (!cs.length) return { net: NaN, mean: NaN, sd: NaN, active: NaN, recorded: 0 };
  let N = 0; for (const c of cs) N += el[c].N;
  const r = cs.filter(c => el[c].active).map(c => el[c].rate), m = r.length ? r.reduce((s, x) => s + x, 0) / r.length : NaN;
  return { net: N / T, mean: m, sd: r.length >= 2 ? sampleSD(r, m) : NaN, active: r.length, recorded: cs.length };
}
