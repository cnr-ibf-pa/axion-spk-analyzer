/* =====================================================================
 * REFERENCE ONLY: not loaded by any page.
 *
 * The PCA + spike-sorting code of the former spk_viewer.html, copied
 * verbatim (lines 817–1365 of the last committed version, commit 6fc708d),
 * kept as the starting point for the PCA step of the MEA pipeline.
 * See docs/pca_spike_sorting.md for what it does, how it was checked, and
 * what it depends on from the old page (F, S, $, C, fmt, setup, frame, ...).
 * The Python twin is python/axion_sort.py.
 * ===================================================================== */
/* =====================================================================
 * Spike sorting: a port of the SAMS automatic pipeline
 * (Ren et al. 2025, bioRxiv doi:10.1101/2025.02.08.637245; src/process_electrode.m)
 * Each electrode is sorted on its own:
 *  1. PCA of the waveforms
 *  2. initial clusters on PC1-PC2: k-means (k-means++ seeding, best of 3 runs) for
 *     K = 2...maxK, K picked by the Davies-Bouldin index
 *     [SAMS uses spectral clustering here; k-means keeps it fast in a browser]
 *  3. outliers removed per cluster: waveform RMSE more than sdCut SD from the cluster mean RMSE
 *  4. under-sorting fix: Hartigan's dip test on PC1 and PC2 of every cluster; when either
 *     is multimodal (p < 0.05) the cluster is clustered again (K = 2...maxK, Davies-Bouldin).
 *     Clusters of 10 spikes or fewer are dropped.
 *  5. template merge: clusters whose min-max normalised mean waveforms are closer than
 *     mergeThr by dynamic time warping, and whose peaks and troughs lie within
 *     (1 - overlap) x height of each other, are merged (repeated until nothing changes)
 *  6. final clean-up per unit: iterative RMSE outlier removal, the second spike of every
 *     pair closer than the refractory period removed, units slower than minRate dropped
 * Spikes removed along the way are "unsorted" (unit 0).
 * ===================================================================== */
const SORT_DEFAULTS = { mergeThr: 1.5, overlap: 0.7, refractoryMs: 1.5, sdCut: 3, minRate: 0.1, maxK: 5 };

function mulberry(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

// symmetric eigen-decomposition (cyclic Jacobi); eigenvalues descending, vecs[k] = k-th eigenvector
function eigSym(A, n) {
  const a = Float64Array.from(A), v = new Float64Array(n * n);
  for (let i = 0; i < n; i++) v[i * n + i] = 1;
  let scale = 0; for (const x of a) scale += x * x;
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p * n + q] ** 2;
    if (off <= 1e-26 * scale) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      const apq = a[p * n + q]; if (Math.abs(apq) < 1e-300) continue;
      const theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
      const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const kp = a[k * n + p], kq = a[k * n + q]; a[k * n + p] = c * kp - s * kq; a[k * n + q] = s * kp + c * kq; }
      for (let k = 0; k < n; k++) { const pk = a[p * n + k], qk = a[q * n + k]; a[p * n + k] = c * pk - s * qk; a[q * n + k] = s * pk + c * qk; }
      for (let k = 0; k < n; k++) { const kp = v[k * n + p], kq = v[k * n + q]; v[k * n + p] = c * kp - s * kq; v[k * n + q] = s * kp + c * kq; }
    }
  }
  const order = [...Array(n).keys()].sort((i, j) => a[j * n + j] - a[i * n + i]);
  return { vals: order.map(i => a[i * n + i]), vecs: order.map(k => Float64Array.from({ length: n }, (_, i) => v[i * n + k])) };
}

// PCA scores (n x nc, row-major) of the rows of X (n x d)
function pcaScores(X, n, d, nc) {
  const mu = new Float64Array(d);
  for (let i = 0; i < n; i++) for (let k = 0; k < d; k++) mu[k] += X[i * d + k];
  for (let k = 0; k < d; k++) mu[k] /= n;
  const cov = new Float64Array(d * d), r = new Float64Array(d);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < d; k++) r[k] = X[i * d + k] - mu[k];
    for (let a = 0; a < d; a++) { const ra = r[a]; for (let b = a; b < d; b++) cov[a * d + b] += ra * r[b]; }
  }
  for (let a = 0; a < d; a++) for (let b = a; b < d; b++) { cov[a * d + b] /= Math.max(1, n - 1); cov[b * d + a] = cov[a * d + b]; }
  const { vals, vecs } = eigSym(cov, d), tot = vals.reduce((s, v) => s + Math.max(v, 0), 0) || 1;
  for (const v of vecs.slice(0, nc)) { let m = 0; for (let k = 1; k < d; k++) if (Math.abs(v[k]) > Math.abs(v[m])) m = k; if (v[m] < 0) for (let k = 0; k < d; k++) v[k] = -v[k]; }
  const S = new Float64Array(n * nc);
  for (let i = 0; i < n; i++) for (let c = 0; c < nc; c++) { const v = vecs[c]; let s = 0; for (let k = 0; k < d; k++) s += (X[i * d + k] - mu[k]) * v[k]; S[i * nc + c] = s; }
  return { scores: S, explained: vals.slice(0, nc).map(v => Math.max(v, 0) / tot) };
}

// k-means (k-means++ seeding, best of `reps` runs by within-cluster sum of squares)
function kmeans(P, n, dim, K, rnd, reps = 3, maxIt = 100) {
  const d2 = (i, C, k) => { let s = 0; for (let e = 0; e < dim; e++) { const z = P[i * dim + e] - C[k * dim + e]; s += z * z; } return s; };
  let best = null;
  for (let rep = 0; rep < reps; rep++) {
    const C = new Float64Array(K * dim), D = new Float64Array(n).fill(Infinity);
    let f = Math.floor(rnd() * n); for (let e = 0; e < dim; e++) C[e] = P[f * dim + e];
    for (let k = 1; k < K; k++) {
      let tot = 0; for (let i = 0; i < n; i++) { const z = d2(i, C, k - 1); if (z < D[i]) D[i] = z; tot += D[i]; }
      let r = rnd() * tot, i = 0; for (; i < n - 1; i++) { r -= D[i]; if (r <= 0) break; }
      for (let e = 0; e < dim; e++) C[k * dim + e] = P[i * dim + e];
    }
    const lab = new Int32Array(n), cnt = new Int32Array(K); let sse = 0;
    for (let it = 0; it < maxIt; it++) {
      let changed = 0; sse = 0;
      for (let i = 0; i < n; i++) {
        let bk = 0, bd = Infinity; for (let k = 0; k < K; k++) { const z = d2(i, C, k); if (z < bd) { bd = z; bk = k; } }
        if (lab[i] !== bk || it === 0) changed++; lab[i] = bk; sse += bd;
      }
      if (!changed) break;
      C.fill(0); cnt.fill(0);
      for (let i = 0; i < n; i++) { cnt[lab[i]]++; for (let e = 0; e < dim; e++) C[lab[i] * dim + e] += P[i * dim + e]; }
      for (let k = 0; k < K; k++) {
        if (cnt[k]) for (let e = 0; e < dim; e++) C[k * dim + e] /= cnt[k];
        else { const j = Math.floor(rnd() * n); for (let e = 0; e < dim; e++) C[k * dim + e] = P[j * dim + e]; } // re-seed an empty cluster
      }
    }
    if (!best || sse < best.sse) best = { lab, C, sse };
  }
  return best;
}

// Davies-Bouldin index (lower = better separated)
function daviesBouldin(P, n, dim, lab, C, K) {
  const s = new Float64Array(K), cnt = new Int32Array(K);
  for (let i = 0; i < n; i++) { let z = 0; for (let e = 0; e < dim; e++) z += (P[i * dim + e] - C[lab[i] * dim + e]) ** 2; s[lab[i]] += Math.sqrt(z); cnt[lab[i]]++; }
  for (let k = 0; k < K; k++) { if (!cnt[k]) return Infinity; s[k] /= cnt[k]; }
  let db = 0;
  for (let i = 0; i < K; i++) {
    let m = 0;
    for (let j = 0; j < K; j++) { if (j === i) continue; let z = 0; for (let e = 0; e < dim; e++) z += (C[i * dim + e] - C[j * dim + e]) ** 2; m = Math.max(m, (s[i] + s[j]) / (Math.sqrt(z) || 1e-12)); }
    db += m;
  }
  return db / K;
}

function clusterBest(P, n, dim, kmin, kmax, rnd) {
  let best = null;
  for (let K = kmin; K <= Math.min(kmax, n - 1); K++) {
    const km = kmeans(P, n, dim, K, rnd), db = daviesBouldin(P, n, dim, km.lab, km.C, K);
    if (!best || db < best.db) best = { lab: km.lab, K, db };
  }
  return best || { lab: new Int32Array(n), K: 1, db: NaN };
}

// Hartigan's dip statistic of a sorted sample (port of AS 217, as in SAMS' HartigansDipTest.m)
function dipStat(xs) {
  const N = xs.length; if (N < 4 || !(xs[N - 1] > xs[0])) return 0;
  const x = new Float64Array(N + 1); for (let i = 0; i < N; i++) x[i + 1] = xs[i];
  // unimodal-by-construction shortcut (ifault 5 in the original)
  let minNeg = Infinity; const sg = new Int8Array(N);
  for (let i = 1; i <= N - 2; i++) { const dd = x[i + 2] - 2 * x[i + 1] + x[i]; sg[i] = dd > 0 ? -1 : dd < 0 ? 1 : 0; if (sg[i] < 0 && i < minNeg) minNeg = i; }
  let posAfter = false, anyPos = false; for (let i = 1; i <= N - 2; i++) if (sg[i] > 0) { anyPos = true; if (i > minNeg) posAfter = true; }
  if (!anyPos || minNeg === Infinity || !posAfter) return 0;
  const mn = new Int32Array(N + 1), mj = new Int32Array(N + 1), gcm = new Int32Array(N + 2), lcm = new Int32Array(N + 2), fn = N;
  let low = 1, high = N, dip = 1 / fn;
  mn[1] = 1;
  for (let j = 2; j <= N; j++) {
    mn[j] = j - 1; let mnj = mn[j], mnmnj = mn[mnj], a = mnj - mnmnj, b = j - mnj;
    while (!(mnj === 1 || (x[j] - x[mnj]) * a < (x[mnj] - x[mnmnj]) * b)) { mn[j] = mnmnj; mnj = mn[j]; mnmnj = mn[mnj]; a = mnj - mnmnj; b = j - mnj; }
  }
  mj[N] = N;
  for (let jk = 1; jk <= N - 1; jk++) {
    const k = N - jk; mj[k] = k + 1; let mjk = mj[k], mjmjk = mj[mjk], a = mjk - mjmjk, b = k - mjk;
    while (!(mjk === N || (x[k] - x[mjk]) * a < (x[mjk] - x[mjmjk]) * b)) { mj[k] = mjmjk; mjk = mj[k]; mjmjk = mj[mjk]; a = mjk - mjmjk; b = k - mjk; }
  }
  for (let guard = 0; guard < 4 * N; guard++) {
    let ic = 1; gcm[1] = high; let igcm1 = gcm[ic]; ic++; gcm[ic] = mn[igcm1];
    while (gcm[ic] > low) { igcm1 = gcm[ic]; ic++; gcm[ic] = mn[igcm1]; }
    const icx = ic;
    ic = 1; lcm[1] = low; let lcm1 = lcm[ic]; ic++; lcm[ic] = mj[lcm1];
    while (lcm[ic] < high) { lcm1 = lcm[ic]; ic++; lcm[ic] = mj[lcm1]; }
    const icv = ic;
    let ig = icx, ih = icv, ix = icx - 1, iv = 2, d = 0;
    if (icx === 2 && icv === 2) d = 1 / fn;
    else {
      for (let g2 = 0; g2 < 4 * N; g2++) {
        const igcmx = gcm[ix], lcmiv = lcm[iv];
        if (!(igcmx > lcmiv)) {
          const lcmiv1 = lcm[iv - 1], a = lcmiv - lcmiv1, b = igcmx - lcmiv1 - 1;
          const dx = (x[igcmx] - x[lcmiv1]) * a / (fn * (x[lcmiv] - x[lcmiv1])) - b / fn;
          ix--;
          if (!(dx < d)) { d = dx; ig = ix + 1; ih = iv; }
        } else {
          const igcm = gcm[ix], ig1 = gcm[ix + 1], a = lcmiv - ig1 + 1, b = igcm - ig1;
          const dx = a / fn - ((x[lcmiv] - x[ig1]) * b) / (fn * (x[igcm] - x[ig1]));
          iv++;
          if (!(dx < d)) { d = dx; ig = ix + 1; ih = iv - 1; }
        }
        if (ix < 1) ix = 1; if (iv > icv) iv = icv;
        if (gcm[ix] === lcm[iv]) break;
      }
    }
    if (d < dip) break;
    let dl = 0;
    if (ig !== icx) for (let j = ig; j <= icx - 1; j++) {
      let temp = 1 / fn; const jb = gcm[j + 1], je = gcm[j];
      if (je - jb > 1 && x[je] !== x[jb]) { const cst = (je - jb) / (fn * (x[je] - x[jb])); for (let jr = jb; jr <= je; jr++) { const t = (jr - jb + 1) / fn - (x[jr] - x[jb]) * cst; if (t > temp) temp = t; } }
      if (dl < temp) dl = temp;
    }
    let du = 0;
    if (ih !== icv) for (let k = ih; k <= icv - 1; k++) {
      let temp = 1 / fn; const kb = lcm[k], ke = lcm[k + 1];
      if (ke - kb > 1 && x[ke] !== x[kb]) { const cst = (ke - kb) / (fn * (x[ke] - x[kb])); for (let kr = kb; kr <= ke; kr++) { const t = (x[kr] - x[kb]) * cst - (kr - kb - 1) / fn; if (t > temp) temp = t; } }
      if (du < temp) du = temp;
    }
    dip = Math.max(dip, dl, du);
    low = gcm[ig]; high = lcm[ih];
  }
  return 0.5 * dip;
}

// dip test p-value against a uniform null (bootstrap, cached per sample size).
// Large clusters are tested on a random subsample of `cap` spikes to keep it fast.
const dipNull = new Map();
function dipTest(vals, rnd, nboot = 200, cap = 2000) {
  let v = Float64Array.from(vals);
  if (v.length > cap) { for (let i = 0; i < cap; i++) { const j = i + Math.floor(rnd() * (v.length - i)); const t = v[i]; v[i] = v[j]; v[j] = t; } v = v.slice(0, cap); }
  v.sort(); const n = v.length, dip = dipStat(v);
  let nul = dipNull.get(n);
  if (!nul) { const r = mulberry(n * 7919), u = new Float64Array(n); nul = new Float64Array(nboot); for (let b = 0; b < nboot; b++) { for (let i = 0; i < n; i++) u[i] = r(); u.sort(); nul[b] = dipStat(u); } dipNull.set(n, nul); }
  let c = 0; for (const z of nul) if (z >= dip) c++;
  return { dip, p: c / nboot };
}

// spikes (local indices g) whose waveform RMSE to the cluster mean is more than sdCut SD from the mean RMSE
function rmseOutliers(W, ns, g, sdCut) {
  const m = g.length; if (m < 3) return { keep: g, out: 0 };
  const mean = new Float64Array(ns); for (const j of g) for (let k = 0; k < ns; k++) mean[k] += W[j * ns + k] / m;
  const r = new Float64Array(m); let mu = 0;
  g.forEach((j, q) => { let s = 0; for (let k = 0; k < ns; k++) s += (W[j * ns + k] - mean[k]) ** 2; r[q] = Math.sqrt(s / ns); mu += r[q] / m; });
  let sd = 0; for (const x of r) sd += (x - mu) ** 2; sd = Math.sqrt(sd / (m - 1));
  const keep = g.filter((_, q) => Math.abs(r[q] - mu) <= sdCut * sd);
  return { keep, out: m - keep.length };
}
function meanWave(W, ns, g) { const m = new Float64Array(ns); for (const j of g) for (let k = 0; k < ns; k++) m[k] += W[j * ns + k] / g.length; return m; }

// dynamic time warping distance (sum of |a-b| along the optimal path, as MATLAB dtw)
function dtw(a, b) {
  const m = b.length; let prev = new Float64Array(m + 1).fill(Infinity), cur = new Float64Array(m + 1); prev[0] = 0;
  for (let i = 1; i <= a.length; i++) {
    cur[0] = Infinity;
    for (let j = 1; j <= m; j++) cur[j] = Math.abs(a[i - 1] - b[j - 1]) + Math.min(prev[j], cur[j - 1], prev[j - 1]);
    [prev, cur] = [cur, prev];
  }
  return prev[m];
}

// SAMS template comparison: which templates should merge (connected components)
function mergePlan(T, prm) {
  let lo = Infinity, hi = -Infinity; for (const t of T) for (const v of t) { if (v < lo) lo = v; if (v > hi) hi = v; }
  const N = T.map(t => Array.from(t, v => (v - lo) / ((hi - lo) || 1))), tol = 1 - prm.overlap;
  const parent = T.map((_, i) => i), root = i => parent[i] === i ? i : (parent[i] = root(parent[i]));
  let closeNotMerged = false; const dist = [];
  for (let i = 0; i < T.length; i++) for (let j = i + 1; j < T.length; j++) {
    const d = dtw(N[i], N[j]); dist.push([i, j, d]);
    const pk1 = Math.max(...N[i]), tr1 = Math.min(...N[i]), pk2 = Math.max(...N[j]), tr2 = Math.min(...N[j]), h1 = pk1 - tr1;
    const same = Math.abs(pk1 - pk2) < h1 * tol && Math.abs(tr1 - tr2) < h1 * tol;
    if (d < prm.mergeThr && same) parent[root(j)] = root(i);
    else if (d < prm.mergeThr) closeNotMerged = true;
  }
  const comps = new Map(); T.forEach((_, i) => { const r = root(i); if (!comps.has(r)) comps.set(r, []); comps.get(r).push(i); });
  return { groups: [...comps.values()], closeNotMerged, dist };
}

function readParams() {
  const num = (id, def, lo, hi) => { const v = parseFloat($(id).value); return Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : def; };
  return { mergeThr: num("pMerge", 1.5, 0, 100), overlap: num("pOverlap", 0.7, 0, 1), refractoryMs: num("pRefr", 1.5, 0, 100),
    sdCut: num("pSd", 3, 1, 20), minRate: num("pRate", 0.1, 0, 1000), maxK: Math.round(num("pK", 5, 2, 8)) };
}

function sortElectrode(c, prm) {
  const t0 = performance.now(), idx = F.perCh[c], n = idx.length, ns = F.header.NumSamplesPerBlock, sc = uV();
  const res = { c, n, prm: { ...prm }, labels: new Int16Array(n), units: [], flags: { hdt: 0, merged: 0, inactive: 0, mua: false }, ms: 0 };
  const minN = Math.max(20, Math.ceil(F.duration * prm.minRate));
  if (n < minN) { res.skipped = `${fmt(n, 0)} spike${n === 1 ? "" : "s"}; at least ${minN} needed`; return res; }
  const W = new Float32Array(n * ns);
  for (let j = 0; j < n; j++) { const o = F.waveOffset(idx[j]); for (let k = 0; k < ns; k++) W[j * ns + k] = F.dv.getInt16(o + 2 * k, true) * sc; }
  const { scores, explained } = pcaScores(W, n, ns, 3);
  res.pcs = Float32Array.from(scores); res.explained = explained;
  const rnd = mulberry(0x5A15 + c);
  const split = g => {
    const P = new Float64Array(g.length * 2); g.forEach((j, q) => { P[2 * q] = scores[3 * j]; P[2 * q + 1] = scores[3 * j + 1]; });
    const { lab, K } = clusterBest(P, g.length, 2, 2, prm.maxK, rnd), out = Array.from({ length: K }, () => []);
    g.forEach((j, q) => out[lab[q]].push(j)); return out.filter(a => a.length);
  };
  // 1-2. initial clusters
  let groups = split(Array.from({ length: n }, (_, j) => j));
  res.initialK = groups.length;
  // 3. outliers
  groups = groups.map(g => rmseOutliers(W, ns, g, prm.sdCut).keep);
  // 4. under-sorting fix (Hartigan's dip test)
  const next = [];
  for (const g of groups) {
    if (g.length <= 10) continue;
    const p1 = dipTest(g.map(j => scores[3 * j]), rnd).p, p2 = dipTest(g.map(j => scores[3 * j + 1]), rnd).p;
    if (Math.min(p1, p2) < 0.05) { const sub = split(g); if (sub.length > 1) { res.flags.hdt++; next.push(...sub); continue; } }
    next.push(g);
  }
  groups = next;
  // 5. template merge (DTW), until stable
  for (let pass = 0; pass < 10 && groups.length > 1; pass++) {
    const plan = mergePlan(groups.map(g => meanWave(W, ns, g)), prm);
    if (plan.groups.length === groups.length) break;
    res.flags.merged += groups.length - plan.groups.length;
    groups = plan.groups.map(comp => comp.flatMap(i => groups[i]));
  }
  // 6. final clean-up
  const refr = prm.refractoryMs / 1000, units = [];
  for (let g of groups) {
    g = g.slice().sort((a, b) => a - b);
    for (let it = 0; it < ns && g.length; it++) { const r = rmseOutliers(W, ns, g, prm.sdCut); if (!r.out) break; g = r.keep; }
    let viol = 0;
    if (g.length > 1) { const drop = new Uint8Array(g.length); for (let q = 1; q < g.length; q++) if (F.times[idx[g[q]]] - F.times[idx[g[q - 1]]] < refr) { drop[q] = 1; viol++; } g = g.filter((_, q) => !drop[q]); }
    if (!g.length || g.length / F.duration <= prm.minRate) { res.flags.inactive++; continue; }
    units.push({ g, viol });
  }
  units.sort((a, b) => b.g.length - a.g.length);
  res.units = units.map(u => unitStats(c, u.g, u.viol));
  relabel(res);
  if (res.units.length > 1) res.flags.mua = mergePlan(res.units.map(u => u.mean), prm).closeNotMerged;
  res.ms = performance.now() - t0;
  return res;
}

// statistics of one unit (g = indices into F.perCh[c], time-ordered)
function unitStats(c, g, viol = 0) {
  const idx = F.perCh[c], ns = F.header.NumSamplesPerBlock, sc = uV(), m = g.length;
  const mean = new Float64Array(ns), m2 = new Float64Array(ns); let noise = 0;
  for (const j of g) { const o = F.waveOffset(idx[j]); for (let k = 0; k < ns; k++) { const v = F.dv.getInt16(o + 2 * k, true) * sc; mean[k] += v; m2[k] += v * v; } noise += F.sd[idx[j]] * 1e6; }
  for (let k = 0; k < ns; k++) { mean[k] /= m; m2[k] = Math.sqrt(Math.max(0, m2[k] / m - mean[k] ** 2)); }
  let mx = -Infinity, mn = Infinity; for (const v of mean) { mx = Math.max(mx, v); mn = Math.min(mn, v); }
  let s1 = 0, s2 = 0, minIsi = Infinity;
  for (let q = 1; q < m; q++) { const d = F.times[idx[g[q]]] - F.times[idx[g[q - 1]]]; s1 += d; s2 += d * d; if (d < minIsi) minIsi = d; }
  const k = m - 1, mi = k ? s1 / k : NaN, cv = k > 1 ? Math.sqrt(Math.max(0, (s2 - k * mi * mi) / (k - 1))) / mi : NaN;
  noise /= m;
  return { g: Int32Array.from(g), n: m, rate: m / F.duration, mean, sd: m2, ptp: mx - mn, peak: Math.abs(mx) > Math.abs(mn) ? mx : mn, noise, snr: (mx - mn) / noise, cv, minIsi, viol };
}
function relabel(res) {
  res.labels.fill(0); res.units.forEach((u, k) => { for (const j of u.g) res.labels[j] = k + 1; });
  const idx = F.perCh[res.c]; for (let j = 0; j < res.n; j++) F.unitOf[idx[j]] = res.labels[j];
}
function storeSort(res) { F.sorting.set(res.c, res); if (!res.units.length) relabel(res); }

/* ---------- curation ---------- */
function mergeUnits(c, list) {
  const res = F.sorting.get(c); if (!res || list.length < 2) return;
  const pick = res.units.filter((_, k) => list.includes(k + 1)), rest = res.units.filter((_, k) => !list.includes(k + 1));
  const g = Int32Array.from(pick.flatMap(u => Array.from(u.g))).sort();
  rest.push(unitStats(c, g, pick.reduce((s, u) => s + u.viol, 0)));
  res.units = rest.sort((a, b) => b.n - a.n); res.curated = true; relabel(res);
}
function discardUnits(c, list) {
  const res = F.sorting.get(c); if (!res) return;
  res.units = res.units.filter((_, k) => !list.includes(k + 1)); res.curated = true; relabel(res);
}

/* ---------- sorting views ---------- */
const SS = { pair: [0, 1], checked: new Set(), hl: 0, running: false, cancel: false, rot: { yaw: -0.6, pitch: 0.35 }, spin: false, spinRaf: 0 };
function unitColor(u) { return u > 0 && u <= 8 ? C.cat[u - 1] : u > 8 ? C["ink-2"] : C.muted; }
function totalUnits() { let t = 0; for (const r of F.sorting.values()) t += r.units.length; return t; }

function drawSort() {
  const res = S.sel >= 0 ? F.sorting.get(S.sel) : null, sorted = F.sorting.size;
  $("sortStatus").textContent = SS.running ? $("sortStatus").textContent : sorted ? `${sorted} of ${F.rows.filter(c => F.counts[c] > 0).length} electrodes with spikes sorted · ${totalUnits()} units` : "";
  $("sortClearBtn").disabled = !sorted || SS.running; $("unitCsvBtn").disabled = !sorted || SS.running; $("sortSelBtn").disabled = SS.running || S.sel < 0;
  const lab = S.sel >= 0 ? F.channels[S.sel].label : "";
  let sub;
  if (!res) sub = S.sel >= 0 ? `${lab} is not sorted yet. Use “Sort selected electrode” or “Sort all electrodes”.` : "";
  else if (res.skipped) sub = `${lab}: too few spikes to sort (${res.skipped}).`;
  else {
    const kept = res.units.reduce((s, u) => s + u.n, 0), un = res.n - kept, f = res.flags;
    const steps = [`${res.initialK} initial clusters`];
    if (f.hdt) steps.push(`${f.hdt} split by the dip test`);
    if (f.merged) steps.push(`${f.merged} merge${f.merged > 1 ? "s" : ""} by template similarity`);
    if (f.inactive) steps.push(`${f.inactive} dropped below ${fmt(res.prm.minRate, 2)} Hz`);
    sub = `${lab}: ${fmt(res.n, 0)} spikes → <b>${res.units.length} unit${res.units.length === 1 ? "" : "s"}</b>, ${fmt(un, 0)} unsorted (${fmt(un / res.n * 100, 1)} %) · ${steps.join(", ")}`
      + (res.curated ? " · edited by hand" : "") + ` · ${fmt(res.ms, 0)} ms`
      + (f.mua && !res.curated ? ` · <span style="color:var(--critical)">⚠ possible multi-unit: two unit templates are similar but their peaks or troughs differ</span>` : "");
  }
  $("sortSub").innerHTML = sub;
  drawSortPca(res); drawSortPca3(res); drawSortTpl(res); drawUnitTable(res);
}

function drawSortPca(res) {
  const cv = $("sortPca"), { ctx, W, H } = setup(cv); hit.pca = null;
  if (!res || !res.pcs) return emptyMsg(ctx, W, H, res ? "Too few spikes to sort" : "Not sorted yet");
  const [a, b] = SS.pair, n = res.n, P = res.pcs, step = Math.max(1, n / 6000);
  const xs = [], ys = []; for (let s = 0; s < n; s += step) { const j = Math.floor(s); xs.push(P[3 * j + a]); ys.push(P[3 * j + b]); }
  const sx = Float64Array.from(xs).sort(), sy = Float64Array.from(ys).sort();
  let x0 = pctl(sx, .002), x1 = pctl(sx, .998), y0 = pctl(sy, .002), y1 = pctl(sy, .998); const px = (x1 - x0) * .06 || 1, py = (y1 - y0) * .06 || 1; x0 -= px; x1 += px; y0 -= py; y1 += py;
  const ex = v => fmt(v * 100, 0);
  const fr = frame(ctx, W, H, { x0, x1, y0, y1, xTicks: niceTicks(x0, x1, 5), yTicks: niceTicks(y0, y1, 5), xLabel: `PC${a + 1} (${ex(res.explained[a])} % of variance)`, yLabel: `PC${b + 1} (${ex(res.explained[b])} %)` });
  ctx.save(); ctx.beginPath(); ctx.rect(fr.L, fr.T, fr.pw, fr.ph); ctx.clip();
  const order = [0, ...res.units.map((_, k) => k + 1)], pts = [], base = Math.max(.25, Math.min(.8, 1500 / Math.min(n, 6000)));
  for (const u of SS.hl ? [...order.filter(u => u !== SS.hl), SS.hl] : order) {
    ctx.fillStyle = unitColor(u); ctx.globalAlpha = SS.hl && u !== SS.hl ? .07 : u === 0 ? base * .6 : base;
    for (let s = 0; s < n; s += step) {
      const j = Math.floor(s); if (res.labels[j] !== u) continue;
      const X = fr.fx(P[3 * j + a]), Y = fr.fy(P[3 * j + b]); ctx.fillRect(X - 1.2, Y - 1.2, 2.4, 2.4); pts.push(X, Y, j);
    }
  }
  ctx.restore(); ctx.globalAlpha = 1;
  hit.pca = { pts, res };
}

/* ---------- 3D feature space: PC1 x PC2 x PC3, orthographic projection, drag to rotate ---------- */
const SS_ROT0 = { yaw: -0.6, pitch: 0.35 };
function pca3Cache(res) { // subsample (<= 6000 spikes) centred and scaled once per sorting result
  if (res._p3) return res._p3;
  const n = res.n, P = res.pcs, step = Math.max(1, n / 6000), js = [];
  for (let s = 0; s < n; s += step) js.push(Math.floor(s));
  const ctr = [0, 0, 0], half = [0, 0, 0];
  for (let a = 0; a < 3; a++) { const v = Float64Array.from(js, j => P[3 * j + a]).sort(), lo = pctl(v, .005), hi = pctl(v, .995); ctr[a] = (lo + hi) / 2; half[a] = (hi - lo) / 2 || 1; }
  const sc = Math.max(...half), xyz = new Float32Array(js.length * 3); // one scale for all axes keeps the true geometry
  js.forEach((j, q) => { for (let a = 0; a < 3; a++) xyz[3 * q + a] = (P[3 * j + a] - ctr[a]) / sc; });
  return (res._p3 = { js: Int32Array.from(js), xyz, ext: half.map(h => Math.min(1.15, h / sc * 1.08)) });
}
function drawSortPca3(res) {
  const cv = $("sortPca3"), { ctx, W, H } = setup(cv); hit.pca3 = null;
  if (!res || !res.pcs) { stopSpin(); return emptyMsg(ctx, W, H, res ? "Too few spikes to sort" : "Not sorted yet"); }
  const { js, xyz, ext } = pca3Cache(res), { yaw, pitch } = SS.rot;
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const s = Math.min(W * .36, H * .42), ox = W / 2, oy = H / 2;
  // PC1 -> x (right), PC2 -> y (up), PC3 -> z (towards the viewer); yaw about y, then pitch about x
  const proj = (x, y, z) => { const x1 = x * cy + z * sy, z1 = -x * sy + z * cy; return [ox + s * x1, oy - s * (y * cp - z1 * sp), y * sp + z1 * cp]; };
  const [ex, ey, ez] = ext, corners = [];
  for (const a of [-1, 1]) for (const b of [-1, 1]) for (const c of [-1, 1]) corners.push([a, b, c]);
  ctx.lineWidth = 1; ctx.strokeStyle = C.grid; ctx.beginPath();
  for (let i = 0; i < 8; i++) for (let j = i + 1; j < 8; j++) {
    const d = corners[i].reduce((t, v, k) => t + (v !== corners[j][k]), 0); if (d !== 1) continue;
    const p = proj(corners[i][0] * ex, corners[i][1] * ey, corners[i][2] * ez), q = proj(corners[j][0] * ex, corners[j][1] * ey, corners[j][2] * ez);
    ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]);
  }
  ctx.stroke();
  // axis triad along the box edges from the (-,-,-) corner, labelled at the far end
  const o = proj(-ex, -ey, -ez), ends = [[ex, -ey, -ez], [-ex, ey, -ez], [-ex, -ey, ez]], pct = v => fmt(v * 100, 0);
  ctx.strokeStyle = C.axis; ctx.lineWidth = 1.5; ctx.fillStyle = C["ink-2"]; ctx.font = `11px ${css("--font")}`; ctx.textBaseline = "middle";
  ends.forEach((e, a) => {
    const p = proj(...e); ctx.beginPath(); ctx.moveTo(o[0], o[1]); ctx.lineTo(p[0], p[1]); ctx.stroke();
    const dx = p[0] - o[0], dy = p[1] - o[1], L = Math.hypot(dx, dy) || 1;
    ctx.textAlign = dx >= 0 ? "left" : "right"; ctx.fillText(`PC${a + 1} (${pct(res.explained[a])} %)`, p[0] + dx / L * 6, p[1] + dy / L * 8);
  });
  // points: unsorted first, then units (largest first); highlighted unit drawn last; nearer points slightly larger
  const m = js.length, sx = new Float32Array(m), syy = new Float32Array(m), sz = new Float32Array(m);
  for (let q = 0; q < m; q++) { const p = proj(xyz[3 * q], xyz[3 * q + 1], xyz[3 * q + 2]); sx[q] = p[0]; syy[q] = p[1]; sz[q] = p[2]; }
  ctx.save(); ctx.beginPath(); ctx.rect(0, 0, W, H); ctx.clip();
  const order = [0, ...res.units.map((_, k) => k + 1)], base = Math.max(.25, Math.min(.8, 1500 / m)), pts = [];
  for (const u of SS.hl ? [...order.filter(u => u !== SS.hl), SS.hl] : order) {
    ctx.fillStyle = unitColor(u); ctx.globalAlpha = SS.hl && u !== SS.hl ? .07 : u === 0 ? base * .6 : base;
    for (let q = 0; q < m; q++) {
      if (res.labels[js[q]] !== u) continue;
      const r = Math.max(.8, Math.min(1.8, 1.2 + sz[q] * .35));
      ctx.fillRect(sx[q] - r, syy[q] - r, 2 * r, 2 * r); pts.push(sx[q], syy[q], js[q]);
    }
  }
  ctx.restore(); ctx.globalAlpha = 1;
  hit.pca3 = { pts, res };
}
function curSort() { return F && S.sel >= 0 ? F.sorting.get(S.sel) : null; }
function stopSpin() { SS.spin = false; cancelAnimationFrame(SS.spinRaf); const b = $("spinBtn"); if (b) b.textContent = "Spin"; }
function startSpin() {
  if (!curSort() || !curSort().pcs) return;
  SS.spin = true; $("spinBtn").textContent = "Stop"; let last = performance.now();
  const step = now => { if (!SS.spin) return; SS.rot.yaw += (now - last) / 1000 * .6; last = now; drawSortPca3(curSort()); if (SS.spin) SS.spinRaf = requestAnimationFrame(step); };
  SS.spinRaf = requestAnimationFrame(step);
}
{
  const cv = $("sortPca3"); let rot = null;
  cv.addEventListener("pointerdown", e => { if (!hit.pca3 || e.button !== 0) return; stopSpin(); hideTip(); rot = { x: e.clientX, y: e.clientY, yaw: SS.rot.yaw, pitch: SS.rot.pitch }; cv.setPointerCapture(e.pointerId); cv.classList.add("grabbing"); e.preventDefault(); });
  cv.addEventListener("pointermove", e => { if (!rot) return; SS.rot.yaw = rot.yaw + (e.clientX - rot.x) * .01; SS.rot.pitch = Math.max(-1.55, Math.min(1.55, rot.pitch + (e.clientY - rot.y) * .01)); drawSortPca3(curSort()); });
  const end = () => { rot = null; cv.classList.remove("grabbing"); };
  cv.addEventListener("pointerup", end); cv.addEventListener("pointercancel", end);
  cv.addEventListener("dblclick", () => { stopSpin(); SS.rot = { ...SS_ROT0 }; drawSortPca3(curSort()); });
  $("view3Btn").addEventListener("click", () => { stopSpin(); SS.rot = { ...SS_ROT0 }; drawSortPca3(curSort()); });
  $("spinBtn").addEventListener("click", e => { SS.spin ? stopSpin() : startSpin(); e.currentTarget.blur(); });
  chartHover("sortPca3", (x, y) => {
    const h = hit.pca3; if (!h || rot || SS.spin) return null; let best = -1, bd = 49;
    for (let q = 0; q < h.pts.length; q += 3) { const d = (h.pts[q] - x) ** 2 + (h.pts[q + 1] - y) ** 2; if (d < bd) { bd = d; best = q; } }
    if (best < 0) return null;
    const j = h.pts[best + 2], u = h.res.labels[j], i = F.perCh[h.res.c][j], P = h.res.pcs;
    return `<b>${u ? "Unit " + u : "Unsorted"}</b><br>spike at ${fmtTime(F.times[i])} (${fmt(F.times[i], 4)} s)<br>PC1 ${fmt(P[3 * j], 1)} · PC2 ${fmt(P[3 * j + 1], 1)} · PC3 ${fmt(P[3 * j + 2], 1)}`;
  });
}

function drawSortTpl(res) {
  const cv = $("sortTpl"), { ctx, W, H } = setup(cv); hit.tpl = null;
  if (!res || !res.units.length) { $("tplLegend").innerHTML = ""; return emptyMsg(ctx, W, H, !res ? "Not sorted yet" : res.skipped ? "Too few spikes to sort" : "No unit passed the quality checks"); }
  const ns = F.header.NumSamplesPerBlock, idx = F.perCh[res.c], tr = Array.from(idx.subarray(0, 200), i => F.trig[i]).sort((a, b) => a - b)[Math.min(idx.length, 200) >> 1];
  const ms = k => (k - tr) / F.fs * 1000;
  let y0 = Infinity, y1 = -Infinity; for (const u of res.units) for (let k = 0; k < ns; k++) { y0 = Math.min(y0, u.mean[k] - u.sd[k]); y1 = Math.max(y1, u.mean[k] + u.sd[k]); }
  const pad = (y1 - y0) * .06; y0 -= pad; y1 += pad;
  const fr = frame(ctx, W, H, { x0: ms(0), x1: ms(ns - 1), y0, y1, xTicks: niceTicks(ms(0), ms(ns - 1), 6), yTicks: niceTicks(y0, y1, 5), xLabel: "Time from trigger (ms)", yLabel: "µV" });
  ctx.save(); ctx.beginPath(); ctx.rect(fr.L, fr.T, fr.pw, fr.ph); ctx.clip();
  const list = res.units.map((u, k) => [u, k + 1]).reverse(); if (SS.hl) list.sort((p, q) => (p[1] === SS.hl) - (q[1] === SS.hl));
  for (const [u, k] of list) {
    const col = unitColor(k), dim = SS.hl && SS.hl !== k;
    ctx.fillStyle = col; ctx.globalAlpha = dim ? .04 : .14; ctx.beginPath();
    for (let s = 0; s < ns; s++) { const X = fr.fx(ms(s)), Y = fr.fy(u.mean[s] + u.sd[s]); s ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); }
    for (let s = ns - 1; s >= 0; s--) ctx.lineTo(fr.fx(ms(s)), fr.fy(u.mean[s] - u.sd[s]));
    ctx.fill();
    ctx.globalAlpha = dim ? .25 : 1; ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath();
    for (let s = 0; s < ns; s++) { const X = fr.fx(ms(s)), Y = fr.fy(u.mean[s]); s ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); }
    ctx.stroke();
  }
  ctx.globalAlpha = 1; ctx.strokeStyle = C.axis; ctx.setLineDash([3, 3]); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(fr.fx(0), fr.T); ctx.lineTo(fr.fx(0), fr.T + fr.ph); ctx.stroke(); ctx.setLineDash([]);
  ctx.restore();
  hit.tpl = { fr, ms, ns, res };
  $("tplLegend").innerHTML = res.units.map((u, k) => `<span><span class="sw" style="background:${unitColor(k + 1)}"></span> Unit ${k + 1}</span>`).join("") + `<span style="color:var(--muted)">line = mean · band = ±1 SD</span>`;
}

function drawUnitTable(res) {
  const t = $("unitTable");
  if (!res || res.skipped) { t.tHead.innerHTML = ""; t.tBodies[0].innerHTML = ""; $("unitBtns").hidden = true; return; }
  const kept = res.units.reduce((s, u) => s + u.n, 0);
  t.tHead.innerHTML = `<tr><th></th><th>Unit</th><th class="n">Spikes</th><th class="n">Rate (Hz)</th><th class="n" title="Peak-to-peak amplitude of the mean waveform">p-p (µV)</th><th class="n" title="Peak-to-peak amplitude / mean noise SD reported by the spike detector">SNR</th><th class="n" title="Coefficient of variation of the inter-spike intervals">ISI CV</th><th class="n" title="Spikes removed for falling inside the refractory period">Refr. removed</th></tr>`;
  t.tBodies[0].innerHTML = res.units.map((u, k) => `<tr data-u="${k + 1}"><td><input type="checkbox" data-u="${k + 1}" ${SS.checked.has(k + 1) ? "checked" : ""} aria-label="Select unit ${k + 1}"></td><td><span class="sw" style="background:${unitColor(k + 1)}"></span> Unit ${k + 1}</td><td class="n">${fmt(u.n, 0)}</td><td class="n">${fmt(u.rate, 3)}</td><td class="n">${fmt(u.ptp, 1)}</td><td class="n">${fmt(u.snr, 1)}</td><td class="n">${fmt(u.cv, 2)}</td><td class="n">${fmt(u.viol, 0)}</td></tr>`).join("")
    + `<tr data-u="0" class="unsorted"><td></td><td><span class="sw" style="background:${C.muted}"></span> Unsorted</td><td class="n">${fmt(res.n - kept, 0)}</td><td class="n">${fmt((res.n - kept) / F.duration, 3)}</td><td colspan="4"></td></tr>`;
  $("unitBtns").hidden = !res.units.length;
  $("mergeBtn").disabled = SS.checked.size < 2; $("noiseBtn").disabled = SS.checked.size < 1;
}

function afterSortChange() { SB.hidden.clear(); drawTiles(); drawTable(); drawRaster(); drawWaveGrid(); drawSort(); drawSpikeBrowser(); }

async function sortAll() {
  const FF = F, prm = readParams(), list = F.rows.filter(c => F.counts[c] > 0);
  SS.running = true; SS.cancel = false; $("sortAllBtn").textContent = "Cancel"; drawSort();
  const t0 = performance.now();
  for (let k = 0; k < list.length; k++) {
    $("sortStatus").textContent = `Sorting ${F.channels[list[k]].label} (${k + 1} of ${list.length})…`;
    await new Promise(r => setTimeout(r, 0));
    if (SS.cancel || F !== FF) break;
    storeSort(sortElectrode(list[k], prm));
    if (list[k] === S.sel) drawSort();
  }
  SS.running = false; $("sortAllBtn").textContent = "Sort all electrodes";
  if (F !== FF) return;
  afterSortChange();
  $("sortStatus").textContent += ` · ${SS.cancel ? "cancelled" : "done"} in ${fmt((performance.now() - t0) / 1000, 1)} s`;
}

$("sortSelBtn").addEventListener("click", () => { if (!F || S.sel < 0) return; SS.checked.clear(); storeSort(sortElectrode(S.sel, readParams())); afterSortChange(); });
$("sortAllBtn").addEventListener("click", () => { if (!F) return; if (SS.running) SS.cancel = true; else sortAll(); });
$("sortClearBtn").addEventListener("click", () => { if (!F) return; F.sorting.clear(); F.unitOf.fill(-1); SS.checked.clear(); SS.hl = 0; afterSortChange(); });
$("pReset").addEventListener("click", () => { const d = SORT_DEFAULTS; $("pMerge").value = d.mergeThr; $("pOverlap").value = d.overlap; $("pRefr").value = d.refractoryMs; $("pSd").value = d.sdCut; $("pRate").value = d.minRate; $("pK").value = d.maxK; });
$("pcPair").addEventListener("change", e => { SS.pair = e.target.value.split(",").map(Number); if (F) drawSortPca(F.sorting.get(S.sel)); });
$("unitTable").addEventListener("change", e => {
  const u = +e.target.dataset.u; if (!u) return; e.target.checked ? SS.checked.add(u) : SS.checked.delete(u);
  $("mergeBtn").disabled = SS.checked.size < 2; $("noiseBtn").disabled = SS.checked.size < 1;
});
$("unitTable").addEventListener("mouseover", e => { const tr = e.target.closest("tr[data-u]"); const u = tr ? +tr.dataset.u : 0; if (u !== SS.hl) { SS.hl = u; const r = F.sorting.get(S.sel); drawSortPca(r); drawSortPca3(r); drawSortTpl(r); } });
$("unitTable").addEventListener("mouseleave", () => { if (SS.hl) { SS.hl = 0; const r = F.sorting.get(S.sel); drawSortPca(r); drawSortPca3(r); drawSortTpl(r); } });
$("mergeBtn").addEventListener("click", () => { mergeUnits(S.sel, [...SS.checked]); SS.checked.clear(); afterSortChange(); });
$("noiseBtn").addEventListener("click", () => { discardUnits(S.sel, [...SS.checked]); SS.checked.clear(); afterSortChange(); });
chartHover("sortPca", (x, y) => {
  const h = hit.pca; if (!h) return null; let best = -1, bd = 64;
  for (let q = 0; q < h.pts.length; q += 3) { const d = (h.pts[q] - x) ** 2 + (h.pts[q + 1] - y) ** 2; if (d < bd) { bd = d; best = q; } }
  if (best < 0) return null;
  const j = h.pts[best + 2], u = h.res.labels[j], i = F.perCh[h.res.c][j];
  return `<b>${u ? "Unit " + u : "Unsorted"}</b><br>spike at ${fmtTime(F.times[i])} (${fmt(F.times[i], 4)} s)`;
});
chartHover("sortTpl", x => {
  const h = hit.tpl; if (!h || x < h.fr.L || x > h.fr.L + h.fr.pw) return null;
  let s = 0, bd = 1e9; for (let k = 0; k < h.ns; k++) { const d = Math.abs(h.fr.fx(h.ms(k)) - x); if (d < bd) { bd = d; s = k; } }
  return `<b>${fmt(h.ms(s), 2)} ms</b>` + h.res.units.map((u, k) => `<br><span class="sw" style="background:${unitColor(k + 1)}"></span> Unit ${k + 1}: ${fmt(u.mean[s], 1)} ± ${fmt(u.sd[s], 1)} µV`).join("");
});
$("unitCsvBtn").addEventListener("click", () => {
  if (!F || !F.sorting.size) return;
  const out = ["electrode,well,treatment,unit,spikes,rate_hz,peak_to_peak_uV,peak_uV,snr,isi_cv,min_isi_ms,refractory_removed,unsorted_on_electrode"];
  for (const c of F.rows) {
    const r = F.sorting.get(c); if (!r) continue;
    const ch = F.channels[c], w = F.wells[ch.well], tr = w && w.what ? `"${String(w.what).replace(/"/g, '""')}"` : "", un = r.n - r.units.reduce((s, u) => s + u.n, 0);
    r.units.forEach((u, k) => out.push([ch.label, ch.well, tr, k + 1, u.n, u.rate.toFixed(4), u.ptp.toFixed(2), u.peak.toFixed(2), u.snr.toFixed(2), Number.isFinite(u.cv) ? u.cv.toFixed(3) : "", Number.isFinite(u.minIsi) ? (u.minIsi * 1000).toFixed(3) : "", u.viol, un].join(",")));
  }
  download(F.name.replace(/\.spk$/i, "") + "_units.csv", out.join("\n"));
});

