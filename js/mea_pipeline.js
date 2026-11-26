"use strict";
/* MEA pipeline: page state, drawing, controls and tabs. The analysis itself is in js/analysis.js,
 * the .spk reader in js/spk_reader.js (both loaded before this file). */
/* =====================================================================
 * State and helpers
 * ===================================================================== */
let F = null;                         // the loaded recording
const S = { sel: -1, t0: 0, t1: 1 };  // selected electrode (index into F.channels), shared time window (s)
const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt = (v, d = 2) => Number.isFinite(v) ? v.toLocaleString(undefined, { maximumFractionDigits: d, minimumFractionDigits: 0 }) : "–";
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
function fmtTime(s) { if (Math.abs(s) >= 120) { const t = Math.round(s), m = Math.floor(t / 60); return `${m}:${String(t - m * 60).padStart(2, "0")}`; } return fmt(s, s < 1 ? 3 : 1) + " s"; }

const C = {};                         // theme colours, refreshed on every full redraw
function readColors() {
  C.stops = [0, 1, 2, 3, 4, 5, 6].map(i => css("--seq-" + i)).map(h => [1, 3, 5].map(k => parseInt(h.slice(k, k + 2), 16)));
  for (const k of ["muted", "grid", "empty", "ink", "ink-2", "accent", "accent-wash", "axis", "surface", "font", "critical", "mean", "nb1", "nb2", "nb3", "meas-n", "meas-r", "meas-d"]) C[k] = css("--" + k);
}
function seqColor(f) { // f in [0,1] -> sequential ramp (light -> dark in light mode, the reverse in dark mode)
  const s = C.stops, x = Math.min(Math.max(f, 0), 1) * 6, i = Math.min(Math.floor(x), 5), r = x - i;
  return `rgb(${s[i].map((v, k) => Math.round(v + (s[i + 1][k] - v) * r)).join(",")})`;
}
const folded = id => $(id).classList.contains("collapsed");
function setup(cv) { // device-pixel-ratio aware canvas
  const dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = cv.clientHeight;
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  const ctx = cv.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  ctx.font = `14px ${C.font}`; return { ctx, W, H };
}
const tip = $("tip");
function showTip(e, html) {
  tip.innerHTML = html; tip.style.display = "block";
  const r = tip.getBoundingClientRect(); let x = e.clientX + 14, y = e.clientY + 14;
  if (x + r.width > innerWidth - 8) x = e.clientX - r.width - 14; if (y + r.height > innerHeight - 8) y = e.clientY - r.height - 14;
  tip.style.left = Math.max(8, x) + "px"; tip.style.top = Math.max(8, y) + "px";
}
function hideTip() { tip.style.display = "none"; }

/* =====================================================================
 * Summary tiles
 * ===================================================================== */
function drawTiles() {
  const h = F.header, active = F.counts.filter(c => c > 0).length, st = h.ExperimentStartTime;
  const treat = [...new Set([...F.recordedWells].map(w => F.wells[w] && F.wells[w].what).filter(Boolean))];
  const tiles = [
    ["Spikes", fmt(F.n, 0), `${fmt(F.n / F.duration, 1)} spikes/s across the plate`],
    ["Active electrodes", `${active} / ${F.recorded.size}`, `with ≥1 spike / recorded (${F.channels.length} on plate)`],
    ["Duration", fmtTime(F.duration), st ? st.toLocaleString() : ""],
    ["Sampling rate", fmt(h.SamplingFrequency / 1000, 3) + " kHz", `${h.NumSamplesPerBlock} samples per spike (${fmt(h.NumSamplesPerBlock / F.fs * 1000, 2)} ms)`],
    ["Plate", F.plateName, `${F.geo.WR}×${F.geo.WC} wells · ${F.geo.EC}×${F.geo.ER} electrodes per well`],
    ["Wells recorded", `${F.recordedWells.size} / ${F.geo.WR * F.geo.WC}`, treat.length ? "Treatments: " + treat.join(", ") : "no treatment labels"],
  ];
  $("recName").textContent = F.name; $("recName").title = F.name;
  $("tiles").innerHTML = tiles.map(([k, v, s]) => `<div class="tile"><div class="k">${esc(k)}</div><div class="v">${esc(v)}</div><div class="s">${esc(s)}</div></div>`).join("");
}

/* =====================================================================
 * Plate map
 * ===================================================================== */
const METRICS = {
  rate:    { name: "Mean firing rate", unit: "Hz", val: c => F.rate[c], d: 2 },
  lograte: { name: "Mean firing rate", unit: "Hz", val: c => F.rate[c], d: 2, log: true },
  count:   { name: "Spike count", unit: "spikes", val: c => F.counts[c], d: 0 },
  amp:     { name: "Median peak-to-peak amplitude", unit: "µV", val: c => F.medPtp[c], d: 1 },
  fano:    { name: "Fano factor", unit: "", val: c => FS.el[c] && FS.el[c].fano, d: 2 },
  cv:      { name: "ISI coefficient of variation (CV)", unit: "", val: c => FS.el[c] && FS.el[c].cv, d: 2 },
  cv2:     { name: "CV2", unit: "", val: c => FS.el[c] && FS.el[c].cv2, d: 2 },
  lv:      { name: "Local variation (LV)", unit: "", val: c => FS.el[c] && FS.el[c].lv, d: 2 },
  medIsi:  { name: "Median ISI", unit: "ms", val: c => FS.el[c] && FS.el[c].medIsi * 1000, d: 1 },
  shortPct:{ name: "ISIs below the short-ISI limit", unit: "%", val: c => FS.el[c] && FS.el[c].shortPct, d: 2 },
};
const hit = { plate: [], wells: [], raster: null, rate: null, wf: null, wg: [], tl: null, plateKey: null, fsRate: null, fsIsi: null, fsAcg: null };
function plateScale() {
  if (P.active) return { m: { name: `Spikes in the last ${fmtWin(P.win)}`, unit: "spikes", val: c => P.counts[c], d: 0 }, lo: 0, hi: P.N, f: v => v / P.N, replay: true };
  const m = METRICS[$("plateMetric").value], vals = [...F.recorded].map(m.val).filter(v => Number.isFinite(v) && v > 0);
  if (!vals.length) return { m, f: () => 0, lo: 0, hi: 1 };
  if (m.log) {
    const lo = Math.min(...vals), hi = Math.max(...vals), a = Math.log10(lo), b = Math.log10(hi);
    return { m, lo, hi, f: v => b > a ? (Math.log10(v) - a) / (b - a) : 1 };
  }
  const hi = Math.max(...vals); return { m, lo: 0, hi, f: v => v / hi };
}
function drawPlate() {
  if (folded("mapCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const cv = $("plate"), { WR, WC, ER, EC } = F.geo;
  // square electrodes: two thirds of the largest size that fits half the panel's width and ~65 % of the window height.
  // The plate's column is then made just wide enough for the plate, and the raster takes the rest of the row.
  const pair = $("mapPair"), stacked = getComputedStyle(pair).gridTemplateColumns.split(" ").length < 2;
  const avail = stacked ? pair.clientWidth : pair.clientWidth / 2, T = 20, gap = Math.max(6, Math.min(14, avail / WC / 14));
  const maxH = Math.max(260, Math.min(640, innerHeight * .65)), PLATE_SCALE = 2 / 3;
  const cell = Math.max(3, PLATE_SCALE * Math.min((avail - 28 - WC * (gap + 8)) / (WC * EC), (maxH - T - 4 - WR * (gap + 8)) / (WR * ER)));
  const ww = cell * EC + gap + 8, wh = cell * ER + gap + 8;
  pair.style.setProperty("--plate-w", Math.ceil(Math.max(340, 28 + WC * ww + 8)) + "px");
  cv.style.height = Math.round(T + WR * wh + 4) + "px";
  const { ctx, W, H } = setup(cv);
  const L = Math.max(24, (W - WC * ww) / 2), sc = plateScale();
  hit.plate = []; hit.wells = [];
  ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.textBaseline = "middle";
  for (let c = 1; c <= WC; c++) ctx.fillText(c, L + (c - .5) * ww, T / 2);
  for (let r = 1; r <= WR; r++) ctx.fillText(String.fromCharCode(64 + r), L / 2, T + (r - .5) * wh);
  for (let r = 1; r <= WR; r++) for (let c = 1; c <= WC; c++) {
    const name = String.fromCharCode(64 + r) + c, x = L + (c - 1) * ww + gap / 2, y = T + (r - 1) * wh + gap / 2, w = ww - gap, h = wh - gap;
    const recd = F.recordedWells.has(name), t = F.wells[name];
    ctx.strokeStyle = recd ? C.axis : C.grid; ctx.lineWidth = 1; ctx.setLineDash(recd ? [] : [3, 3]);
    ctx.beginPath(); ctx.roundRect(x + .5, y + .5, w - 1, h - 1, 6); ctx.stroke(); ctx.setLineDash([]);
    if (t && t.color && recd) { ctx.fillStyle = t.color; ctx.fillRect(x + 6, y - 1, Math.min(28, w - 12), 3); } // treatment colour tab (AxIS colour)
    hit.wells.push({ x, y, w, h, name });
  }
  for (const ch of F.channels) {
    const wx = L + (ch.wc - 1) * ww + gap / 2 + 4, wy = T + (ch.wr - 1) * wh + gap / 2 + 4;
    const cw = (ww - gap - 8) / EC, chh = (wh - gap - 8) / ER, pad = Math.min(1.5, cw * .08);
    const x = wx + (ch.ec - 1) * cw + pad / 2, y = wy + (ch.er - 1) * chh + pad / 2, w = cw - pad, h = chh - pad;
    const rec = F.recorded.has(ch.i), v = sc.m.val(ch.i);
    if (!rec) { ctx.fillStyle = C.empty; ctx.globalAlpha = .45; }
    else if (!(v > 0)) ctx.fillStyle = C.empty;
    else ctx.fillStyle = seqColor(0.08 + 0.92 * Math.min(1, Math.max(0, sc.f(v))));
    ctx.beginPath(); ctx.roundRect(x, y, w, h, Math.min(3, w / 5)); ctx.fill(); ctx.globalAlpha = 1;
    if (ch.i === S.sel) { ctx.strokeStyle = C.ink; ctx.lineWidth = 2; ctx.beginPath(); ctx.roundRect(x - 1.5, y - 1.5, w + 3, h + 3, 4); ctx.stroke(); }
    hit.plate.push({ x, y, w, h, c: ch.i });
  }
  // subtitle + legend (rebuilt only when they change: the map is redrawn on every replay frame)
  const m = sc.m, key = `${sc.replay}|${m.name}|${sc.lo}|${sc.hi}|${C.ink}`;
  if (key === hit.plateKey) return; hit.plateKey = key;
  $("plateTitle").textContent = sc.replay ? "Plate map · replay" : "Plate map";
  $("plateMetric").disabled = !!sc.replay;
  if (sc.replay) {
    $("plateSub").textContent = `Spikes per electrode in the ${fmtWin(P.win)} before the playhead (count, no smoothing) · colours saturate at ${sc.hi} spike${sc.hi > 1 ? "s" : ""} (99th percentile over the recording) · Stop returns to the whole-recording map`;
    $("plateLegend").innerHTML = `<span>1 spike</span><span class="ramp" style="background:linear-gradient(90deg,${[0, .25, .5, .75, 1].map(f => seqColor(.08 + .92 * f)).join(",")})"></span><span>≥ ${sc.hi} spikes</span><span><span class="sw" style="background:var(--empty)"></span> silent in the window</span><span><span class="sw" style="background:var(--empty);opacity:.45"></span> not recorded</span>`;
    return;
  }
  $("plateSub").textContent = `${m.name}${m.log ? " (log scale)" : ""} per electrode over the whole recording · click an electrode to select it (arrow keys move the selection)`;
  const ramp = `linear-gradient(90deg,${[0, .25, .5, .75, 1].map(f => seqColor(.08 + .92 * f)).join(",")})`;
  const treat = [...F.recordedWells].map(w => F.wells[w]).filter(t => t && t.what);
  const uniq = [...new Map(treat.map(t => [t.what + t.color, t])).values()];
  const sig = v => v > 0 && v < 0.1 ? String(+v.toPrecision(2)) : fmt(v, m.d); // keep small rates readable (e.g. 0.0033 Hz)
  $("plateLegend").innerHTML = `<span>${sig(sc.lo)} ${esc(m.unit)}</span><span class="ramp" style="background:${ramp}"></span><span>${sig(sc.hi)} ${esc(m.unit)}</span>`
    + `<span><span class="sw" style="background:var(--empty)"></span> no spikes</span><span><span class="sw" style="background:var(--empty);opacity:.45"></span> not recorded</span>`
    + uniq.map(t => `<span><span class="sw" style="background:${esc(t.color || "transparent")};height:4px;border:0"></span> ${esc(t.what)}</span>`).join("");
}
function plateHit(e) {
  const r = $("plate").getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
  const el = hit.plate.find(b => x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h);
  return el ? { el } : { well: hit.wells.find(b => x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h) };
}
function electrodeTip(c) {
  const ch = F.channels[c], t = F.wells[ch.well];
  if (!F.recorded.has(c)) return `<b>${ch.label}</b><br>not recorded`;
  const now = P.active ? `<br>${fmt(P.counts[c], 0)} spike${P.counts[c] === 1 ? "" : "s"} in ${fmtClock(Math.max(0, P.t - P.win))} – ${fmtClock(P.t)}` : "";
  return `<b>${ch.label}</b>${t && t.what ? ` · ${esc(t.what)}` : ""}<br>${fmt(F.counts[c], 0)} spikes · ${fmt(F.rate[c], 3)} Hz`
    + (F.counts[c] ? `<br>median peak-to-peak ${fmt(F.medPtp[c], 1)} µV · noise SD ${fmt(F.medNoise[c], 2)} µV` : "") + now;
}
function wellTip(name) {
  const t = F.wells[name], chans = F.channels.filter(c => c.well === name && F.recorded.has(c.i));
  if (!chans.length) return `<b>Well ${name}</b><br>not recorded`;
  const n = chans.reduce((s, c) => s + F.counts[c.i], 0), act = chans.filter(c => F.counts[c.i] > 0).length;
  return `<b>Well ${name}</b>${t && t.what ? ` · ${esc(t.what)}` : ""}<br>${act} of ${chans.length} electrodes active · ${fmt(n, 0)} spikes (${fmt(n / F.duration, 1)} /s)`;
}

/* ---------- selection ---------- */
function select(c) {
  if (c === undefined || c < 0 || !F.recorded.has(c)) return;
  S.sel = c; drawPlate(); drawSelection(); drawRaster(); drawWaveforms(); drawPlayheads(); if (FS.prm) { drawStatsTable(); drawStatCharts(); }
}
function drawSelection() {
  const c = S.sel; if (c < 0) { $("selbar").innerHTML = ""; return; }
  const ch = F.channels[c], t = F.wells[ch.well];
  $("selbar").innerHTML = `<b>${ch.label}</b><span>well ${ch.well}${t && t.what ? ` · ${esc(t.what)}` : ""}</span>`
    + `<span>${fmt(F.counts[c], 0)} spikes · ${fmt(F.rate[c], 3)} Hz</span>`
    + (F.counts[c] ? `<span>median peak-to-peak ${fmt(F.medPtp[c], 1)} µV</span><span>detector noise SD ${fmt(F.medNoise[c], 2)} µV · threshold ${fmt(F.thrMult[c], 1)} × SD</span>` : "<span>no spikes</span>");
}
function moveSel(dx, dy) { // arrow keys: nearest recorded electrode in that direction on the plate grid
  const g = F.geo, pos = c => { const ch = F.channels[c]; return [(ch.wc - 1) * g.EC + ch.ec, (ch.wr - 1) * g.ER + ch.er]; };
  const [x0, y0] = pos(S.sel); let best = -1, bd = Infinity;
  for (const c of F.recorded) {
    const [x, y] = pos(c), ax = (x - x0) * dx + (y - y0) * dy; if (ax <= 0) continue;
    const d = ax + 3 * Math.abs(dx ? y - y0 : x - x0); if (d < bd) { bd = d; best = c; }
  }
  if (best >= 0) select(best);
}

/* =====================================================================
 * Spike raster + population rate (shared time window S.t0 – S.t1)
 * ===================================================================== */
function parseTime(str) { // "385.2", "6:25.2", "1:02:03.5" -> seconds
  const t = String(str).trim().replace(",", "."); if (!t) return NaN;
  if (!t.includes(":")) return /^\d*\.?\d+(e[+-]?\d+)?$/i.test(t) ? +t : NaN;
  const parts = t.split(":"); if (parts.length > 3 || parts.some(p => !/^\d*\.?\d+$/.test(p))) return NaN;
  return parts.reduce((a, p) => a * 60 + +p, 0);
}
function niceTicks(a, b, n = 6) {
  const span = b - a; if (!(span > 0)) return [a];
  const raw = span / n, mag = 10 ** Math.floor(Math.log10(raw)), e = raw / mag, step = (e >= 5 ? 10 : e >= 2 ? 5 : e >= 1 ? 2 : 1) * mag;
  const out = []; for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-9; v += step) out.push(+v.toPrecision(12)); return out;
}
function timeLabel(v, span) { return span >= 120 ? fmtTime(v) : fmt(v, span < 0.05 ? 4 : span < 2 ? 3 : span < 20 ? 2 : 1) + " s"; }
function rasterRows() { // rows shown: recorded electrodes in plate order, filtered by well and activity
  const well = $("rWell").value, all = $("rRows").value === "all";
  return F.rows.filter(c => (well === "*" || F.channels[c].well === well) && (all || F.counts[c] > 0 || c === S.sel));
}
function setWindow(a, b) {
  const minSpan = 2 / F.fs;                                   // at least two samples
  a = Math.max(0, Math.min(a, F.duration - minSpan)); b = Math.min(F.duration, Math.max(b, a + minSpan));
  S.t0 = a; S.t1 = b; drawRaster(); drawWaveforms(); drawPlayheads(); if (FS.prm && $("fsPeriod").value === "window") drawStats();
}
const RL = 64, RR = 12;                                       // left label column / right margin (px), shared by both canvases
function drawRaster() {
  if (folded("mapCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const rows = rasterRows(), rowOf = new Map(rows.map((c, r) => [c, r])), span = S.t1 - S.t0;
  // population rate strip
  const cvR = $("rate"), { ctx: cr, W } = setup(cvR), pw = W - RL - RR, HR = cvR.clientHeight, TR = 6, BR = 4, phR = HR - TR - BR;
  const fx = t => RL + (t - S.t0) / span * pw;
  const nb = Math.max(10, Math.min(400, Math.floor(pw / 3))), bw = span / nb, bins = new Float64Array(nb);
  const lo = lowerBound(F.times, S.t0), hi = lowerBound(F.times, S.t1 + 1e-12);
  for (let i = lo; i < hi; i++) { if (!rowOf.has(F.ch[i])) continue; const b = Math.min(nb - 1, Math.floor((F.times[i] - S.t0) / bw)); bins[b]++; }
  for (let b = 0; b < nb; b++) bins[b] /= bw;
  const ymax = Math.max(1e-9, ...bins), yt = niceTicks(0, ymax, 2), ytop = yt[yt.length - 1] >= ymax ? yt[yt.length - 1] : ymax;
  cr.textAlign = "right"; cr.textBaseline = "middle"; cr.fillStyle = C.muted; cr.font = `13px ${C.font}`;
  for (const v of yt) { const y = Math.round(TR + phR - v / ytop * phR) + .5; cr.strokeStyle = C.grid; cr.lineWidth = 1; cr.beginPath(); cr.moveTo(RL, y); cr.lineTo(RL + pw, y); cr.stroke(); cr.fillText(fmt(v, v < 10 ? 1 : 0), RL - 6, y); }
  cr.save(); cr.translate(10, TR + phR / 2); cr.rotate(-Math.PI / 2); cr.textAlign = "center"; cr.fillText("spikes/s", 0, 0); cr.restore();
  cr.fillStyle = C.accent; const bpx = pw / nb;
  for (let b = 0; b < nb; b++) if (bins[b] > 0) { const h = bins[b] / ytop * phR; cr.fillRect(RL + b * bpx + (bpx > 3 ? .5 : 0), TR + phR - h, Math.max(1, bpx - (bpx > 3 ? 1 : 0)), h); }
  hit.rate = { nb, bw, bins, pw };
  // raster
  const cv = $("raster");
  const wantH = Math.round(Math.min(1400, Math.max(innerWidth <= 600 ? 300 : 420, rows.length * 8 + 40))) + "px";
  if (cv.style.height !== wantH) cv.style.height = wantH;
  const { ctx, H } = setup(cv), T = 4, B = 30, ph = H - T - B, rh = ph / Math.max(1, rows.length);
  let lastWell = null, wellY = -1e9;
  ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.font = `13px ${C.font}`;
  rows.forEach((c, r) => {
    const ch = F.channels[c], y = T + r * rh, first = ch.well !== lastWell;
    if (c === S.sel) { ctx.fillStyle = C["accent-wash"]; ctx.fillRect(RL, y, pw, rh); }
    if (first) {
      if (lastWell !== null) { ctx.strokeStyle = C.axis; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(0, Math.round(y) + .5); ctx.lineTo(RL + pw, Math.round(y) + .5); ctx.stroke(); }
      lastWell = ch.well;
    }
    if (rh >= 10) { ctx.fillStyle = c === S.sel ? C.ink : C.muted; ctx.fillText(ch.label, RL - 6, y + rh / 2); }
    else if (first) { ctx.fillStyle = C["ink-2"]; ctx.textBaseline = "top"; ctx.fillText(ch.well, RL - 6, y + 2); ctx.textBaseline = "middle"; wellY = y; }
    else if (c === S.sel && y - wellY > 14) { ctx.fillStyle = C.ink; ctx.fillText(ch.label, RL - 6, y + rh / 2); }
  });
  const tickH = Math.max(1, rh - (rh > 5 ? 2 : 0)), tw = Math.max(1, Math.min(3, pw / (span * F.fs)));
  ctx.save(); ctx.beginPath(); ctx.rect(RL, T, pw, ph); ctx.clip();
  ctx.fillStyle = C["ink-2"]; ctx.globalAlpha = .85;
  for (let i = lo; i < hi; i++) { const c = F.ch[i], r = rowOf.get(c); if (r === undefined || c === S.sel) continue; ctx.fillRect(Math.floor(fx(F.times[i])), T + r * rh + (rh - tickH) / 2, tw, tickH); }
  ctx.globalAlpha = 1;
  if (rowOf.has(S.sel)) {
    ctx.fillStyle = C.accent; const r = rowOf.get(S.sel), idx = F.perCh[S.sel];
    for (let k = lowerBound(Array.from(idx, i => F.times[i]), S.t0); k < idx.length && F.times[idx[k]] <= S.t1; k++) ctx.fillRect(Math.floor(fx(F.times[idx[k]])), T + r * rh + (rh - tickH) / 2, tw, tickH);
  }
  ctx.restore();
  ctx.strokeStyle = C.axis; ctx.beginPath(); ctx.moveTo(RL, T + ph + .5); ctx.lineTo(RL + pw, T + ph + .5); ctx.stroke();
  ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const v of niceTicks(S.t0, S.t1, Math.max(2, Math.floor(pw / 110)))) { const x = fx(v); if (x >= RL - 1 && x <= RL + pw + 1) ctx.fillText(timeLabel(v, span), x, T + ph + 7); }
  hit.raster = { T, ph, rh, rows, pw, fx };
  // toolbar state
  const nSpk = rows.reduce((s, c) => s + F.counts[c], 0), inWin = bins.reduce((s, v) => s + v * bw, 0);
  $("rInfo").textContent = `${rows.length} electrodes · ${fmt(Math.round(inWin), 0)} spikes in window`;
  $("rReset").disabled = S.t0 <= 0 && S.t1 >= F.duration;
  for (const [id, v] of [["rFrom", S.t0], ["rTo", S.t1]]) { const el = $(id); if (document.activeElement !== el) { el.value = +v.toFixed(4); el.removeAttribute("aria-invalid"); } }
  $("rasterSub").textContent = `Every spike of every ${$("rRows").value === "all" ? "recorded electrode" : "electrode with spikes"}${$("rWell").value === "*" ? "" : " of well " + $("rWell").value}, grouped by well · top: spikes/s summed over the rows shown (bins of ${bw < 1 ? fmt(bw * 1000, bw < 0.01 ? 2 : 0) + " ms" : fmt(bw, 2) + " s"}) · drag to zoom · click a row to select its electrode`;
  void nSpk;
}
function fillWellSelect() {
  const sel = $("rWell"), wells = [...F.recordedWells].sort((a, b) => a.charCodeAt(0) - b.charCodeAt(0) || parseInt(a.slice(1)) - parseInt(b.slice(1)));
  sel.innerHTML = `<option value="*">All wells</option>` + wells.map(w => `<option value="${w}">${w}${F.wells[w] && F.wells[w].what ? " · " + esc(F.wells[w].what) : ""}</option>`).join("");
}
{
  const wrap = $("rWrap"), br = $("rBrush"); let d = null;
  const px = e => e.clientX - wrap.getBoundingClientRect().left;
  const tOf = x => S.t0 + Math.min(Math.max((x - RL) / hit.raster.pw, 0), 1) * (S.t1 - S.t0);
  wrap.addEventListener("pointerdown", e => { if (!F || e.button !== 0) return; d = { x: px(e) }; wrap.setPointerCapture(e.pointerId); e.preventDefault(); });
  wrap.addEventListener("pointermove", e => {
    if (!F) return;
    if (d) { const x = px(e), a = Math.max(RL, Math.min(d.x, x)), b = Math.min(RL + hit.raster.pw, Math.max(d.x, x)); if (Math.abs(x - d.x) > 4) { hideTip(); Object.assign(br.style, { display: "block", left: a + "px", width: (b - a) + "px" }); } return; }
    const rc = $("raster").getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top, g = hit.raster;
    if (x < RL || x > RL + g.pw) return hideTip();
    if (y < 0) { // over the rate strip
      const b = Math.floor((x - RL) / (g.pw / hit.rate.nb)); if (b < 0 || b >= hit.rate.nb) return hideTip();
      const t = S.t0 + b * hit.rate.bw; return showTip(e, `<b>${timeLabel(t, S.t1 - S.t0)} – ${timeLabel(t + hit.rate.bw, S.t1 - S.t0)}</b><br>${fmt(hit.rate.bins[b], 1)} spikes/s over ${g.rows.length} electrodes`);
    }
    const r = Math.floor((y - g.T) / g.rh); if (r < 0 || r >= g.rows.length) return hideTip();
    const c = g.rows[r], t = tOf(x - 0 + (wrap.getBoundingClientRect().left - rc.left)), idx = F.perCh[c];
    const a = lowerBound(Array.from(idx, i => F.times[i]), S.t0), b = lowerBound(Array.from(idx, i => F.times[i]), S.t1 + 1e-12);
    showTip(e, `<b>${F.channels[c].label}</b>${F.wells[F.channels[c].well] && F.wells[F.channels[c].well].what ? " · " + esc(F.wells[F.channels[c].well].what) : ""} · ${timeLabel(t, S.t1 - S.t0)}<br>${fmt(b - a, 0)} spikes in window (${fmt((b - a) / (S.t1 - S.t0), 2)} Hz)`);
  });
  const end = (e, cancel) => {
    if (!d) return; const x = px(e), x0 = d.x; d = null; br.style.display = "none"; if (cancel || !F) return;
    if (Math.abs(x - x0) > 4) return setWindow(tOf(Math.min(x, x0)), tOf(Math.max(x, x0)));
    const rc = $("raster").getBoundingClientRect(), y = e.clientY - rc.top, g = hit.raster, r = Math.floor((y - g.T) / g.rh);
    if (y >= 0 && r >= 0 && r < g.rows.length) select(g.rows[r]);
  };
  wrap.addEventListener("pointerup", e => end(e, false)); wrap.addEventListener("pointercancel", e => end(e, true));
  wrap.addEventListener("mouseleave", hideTip);
  const shift = f => { const s = S.t1 - S.t0; let a = S.t0 + f * s; a = Math.max(0, Math.min(a, F.duration - s)); setWindow(a, a + s); };
  $("rLeft").addEventListener("click", () => F && shift(-.5));
  $("rRight").addEventListener("click", () => F && shift(.5));
  $("rOut").addEventListener("click", () => { if (!F) return; const m = (S.t0 + S.t1) / 2, s = (S.t1 - S.t0) * 2; setWindow(m - s / 2, m + s / 2); });
  $("rReset").addEventListener("click", () => F && setWindow(0, F.duration));
  const apply = () => {
    if (!F) return; const a = parseTime($("rFrom").value), b = parseTime($("rTo").value);
    $("rFrom").setAttribute("aria-invalid", String(!Number.isFinite(a))); $("rTo").setAttribute("aria-invalid", String(!Number.isFinite(b)));
    if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) return;
    $("rFrom").blur(); $("rTo").blur(); setWindow(Math.min(a, b), Math.max(a, b));
  };
  for (const id of ["rFrom", "rTo"]) {
    $(id).addEventListener("keydown", e => { if (e.key === "Enter") { e.preventDefault(); apply(); } else if (e.key === "Escape") { e.target.blur(); drawRaster(); } });
    $(id).addEventListener("change", apply);
  }
  $("rWell").addEventListener("change", () => F && drawRaster());
  $("rRows").addEventListener("change", () => F && drawRaster());
}

/* =====================================================================
 * Waveforms per electrode (spikes in the raster's time window)
 * Samples are drawn exactly as stored in the file: straight lines between
 * recorded samples, no interpolation, filtering or smoothing. Alignment only
 * shifts whole samples.
 * ===================================================================== */
function windowSpikes(c) { // spikes of electrode c inside [S.t0, S.t1]
  const idx = F.perCh[c], t = F.times; let lo = 0, hi = idx.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (t[idx[m]] < S.t0) lo = m + 1; else hi = m; }
  let e = lo; hi = idx.length; while (e < hi) { const m = (e + hi) >> 1; if (t[idx[m]] <= S.t1) e = m + 1; else hi = m; }
  return idx.subarray(lo, e);
}
function pickSpikes(idx, max) { if (idx.length <= max) return Array.from(idx); const out = [], step = idx.length / max; for (let k = 0; k < max; k++) out.push(idx[Math.floor(k * step)]); return out; }
const WF = { cache: null };
function wfSettings() { return { n: $("wfN").value === "all" ? Infinity : +$("wfN").value, mode: $("wfAlign").value, view: $("wfView").value, shared: $("wfScale").value === "shared" }; }
function wfWell() { const v = $("wfWell").value; return v === "sel" ? (S.sel >= 0 ? F.channels[S.sel].well : null) : v; }

function drawWaveforms() {
  if (folded("wfCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const st0 = wfSettings(), well = wfWell();
  const chans = F.channels.filter(ch => ch.well === well);
  // per-electrode stats (grid: up to 60 spikes each; detail: the chosen number)
  const per = new Map();
  for (const ch of chans) if (F.recorded.has(ch.i)) { const idx = windowSpikes(ch.i); per.set(ch.i, { all: idx.length, st: waveStats(F, pickSpikes(idx, 60), st0.mode) }); }
  const selIdx = S.sel >= 0 ? windowSpikes(S.sel) : new Int32Array(0);
  const det = S.sel >= 0 ? waveStats(F, pickSpikes(selIdx, st0.view === "density" ? Infinity : st0.n), st0.mode) : null;
  // shared amplitude range: robust (0.2–99.8 %) over the shown spikes of the well
  let pooled = [];
  for (const { st } of per.values()) if (st) for (const w of st.waves) for (let k = 0; k < w.length; k += 2) pooled.push(w[k]);
  pooled.sort((a, b) => a - b);
  const q = p => pooled.length ? pooled[Math.min(pooled.length - 1, Math.floor(p * pooled.length))] : 0;
  const gm = Math.max(10, Math.abs(q(.002)), Math.abs(q(.998))) * 1.1;
  WF.cache = { per, det, gm, well, chans };
  drawWaveDetail(det, selIdx.length, st0);
  drawWaveGrid(st0);
  const span = `${timeLabel(S.t0, S.t1 - S.t0)} – ${timeLabel(S.t1, S.t1 - S.t0)}`;
  $("wfSub").textContent = `Spikes in the raster window (${span}) · samples drawn as recorded (straight lines between samples, no smoothing)${st0.view === "density" ? " · density: every spike, normalised per time column (darker = more spikes at that voltage at that instant)" : ""} · time 0 = ${{ trig: "threshold crossing (as stored by the detector)", neg: "most negative sample", pos: "most positive sample", abs: "largest absolute sample" }[st0.mode]}`;
}

function drawWaveDetail(st, nAll, set) {
  const cv = $("wfMain"), { ctx, W, H } = setup(cv); hit.wf = null;
  if (S.sel < 0 || !st) { $("wfStats").innerHTML = ""; ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.fillText(S.sel < 0 ? "Select an electrode" : "No spikes in this window", W / 2, H / 2); return; }
  const c = S.sel, ms = j => (j + st.lagMin) / F.fs * 1000;
  const sh = set.shared ? WF.cache.gm : null;
  let lo, hi;
  if (sh) { lo = -sh; hi = sh; }
  else { const v = []; for (const w of st.waves) for (const x of w) v.push(x); v.sort((a, b) => a - b); lo = v[Math.floor(v.length * .002)]; hi = v[Math.min(v.length - 1, Math.floor(v.length * .998))]; const p = (hi - lo) * .1 || 5; lo -= p; hi += p; }
  const thr = F.thrMult[c] * F.medNoise[c];
  lo = Math.min(lo, -thr * 1.1); hi = Math.max(hi, thr * 1.1);
  const Lm = 52, Rm = 12, Tm = 10, Bm = 34, pw = W - Lm - Rm, ph = H - Tm - Bm;
  const x0 = ms(0), x1 = ms(st.mean.length - 1), fx = t => Lm + (t - x0) / (x1 - x0) * pw, fy = v => Tm + ph - (v - lo) / (hi - lo) * ph;
  ctx.font = `13px ${C.font}`; ctx.fillStyle = C.muted; ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const v of niceTicks(lo, hi, 6)) { const y = Math.round(fy(v)) + .5; ctx.strokeStyle = v === 0 ? C.axis : C.grid; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(Lm, y); ctx.lineTo(Lm + pw, y); ctx.stroke(); ctx.fillText(fmt(v, 0), Lm - 6, y); }
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const v of niceTicks(x0, x1, 7)) ctx.fillText(fmt(v, 2), fx(v), Tm + ph + 6);
  ctx.fillText("ms from time 0", Lm + pw / 2, H - 13);
  ctx.save(); ctx.translate(12, Tm + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("µV", 0, 0); ctx.restore();
  ctx.save(); ctx.beginPath(); ctx.rect(Lm, Tm, pw, ph); ctx.clip();
  if (set.view === "density") { // 2-D histogram of every sample of every spike
    const gx = Math.max(1, Math.floor(pw / 3)), gy = Math.max(1, Math.floor(ph / 3)), G = new Float64Array(gx * gy);
    st.waves.forEach((w, qq) => { for (let k = 0; k < w.length - 1; k++) { // rasterise each straight segment between consecutive samples
      const xa = fx(ms(k - st.al[qq] - st.lagMin)), xb = fx(ms(k + 1 - st.al[qq] - st.lagMin)), ya = fy(w[k]), yb = fy(w[k + 1]);
      const steps = Math.max(1, Math.ceil(Math.max(Math.abs(xb - xa), Math.abs(yb - ya)) / 2));
      for (let s = 0; s < steps; s++) { const f = s / steps, X = Math.floor((xa + (xb - xa) * f - Lm) / 3), Y = Math.floor((ya + (yb - ya) * f - Tm) / 3); if (X >= 0 && X < gx && Y >= 0 && Y < gy) G[Y * gx + X]++; }
    } });
    // normalise each time column (share of the spikes at each voltage at that instant), colour up to the 98th percentile
    for (let X = 0; X < gx; X++) { let sum = 0; for (let Y = 0; Y < gy; Y++) sum += G[Y * gx + X]; if (sum) for (let Y = 0; Y < gy; Y++) G[Y * gx + X] /= sum; }
    const nz = Array.from(G).filter(v => v > 0).sort((a, b) => a - b), cap = nz.length ? nz[Math.floor(nz.length * .98)] : 1;
    for (let Y = 0; Y < gy; Y++) for (let X = 0; X < gx; X++) { const v = G[Y * gx + X]; if (!v) continue; ctx.fillStyle = seqColor(.04 + .96 * Math.sqrt(Math.min(1, v / cap))); ctx.fillRect(Lm + X * 3, Tm + Y * 3, 3, 3); }
  } else {
    ctx.strokeStyle = C["ink-2"]; ctx.lineWidth = 1; ctx.globalAlpha = Math.max(.05, Math.min(.5, 8 / st.n));
    st.waves.forEach((w, qq) => { ctx.beginPath(); for (let k = 0; k < w.length; k++) { const X = fx(ms(k - st.al[qq] - st.lagMin)), Y = fy(w[k]); k ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); } ctx.stroke(); });
    ctx.globalAlpha = 1;
  }
  // detection threshold (± threshold multiplier × noise SD, median over the electrode's spikes) and time 0, on top of the spikes
  ctx.strokeStyle = C.critical; ctx.globalAlpha = .85; ctx.setLineDash([6, 4]); ctx.lineWidth = 2;
  for (const v of [-thr, thr]) { ctx.beginPath(); ctx.moveTo(Lm, Math.round(fy(v))); ctx.lineTo(Lm + pw, Math.round(fy(v))); ctx.stroke(); }
  ctx.setLineDash([]); ctx.globalAlpha = 1;
  ctx.strokeStyle = C.axis; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(Math.round(fx(0)) + .5, Tm); ctx.lineTo(Math.round(fx(0)) + .5, Tm + ph); ctx.stroke(); ctx.setLineDash([]);
  // mean ± SD over the spikes (lags covered by at least half of them)
  ctx.fillStyle = C.mean; ctx.globalAlpha = .16; ctx.beginPath(); let started = false;
  const okJ = []; for (let j = 0; j < st.mean.length; j++) if (st.ok[j]) okJ.push(j);
  okJ.forEach(j => { const X = fx(ms(j)), Y = fy(st.mean[j] + st.sd[j]); started ? ctx.lineTo(X, Y) : (ctx.moveTo(X, Y), started = true); });
  [...okJ].reverse().forEach(j => ctx.lineTo(fx(ms(j)), fy(st.mean[j] - st.sd[j]))); ctx.fill(); ctx.globalAlpha = 1;
  ctx.strokeStyle = C.mean; ctx.lineWidth = 3; ctx.lineJoin = "round"; ctx.beginPath(); okJ.forEach((j, q) => { const X = fx(ms(j)), Y = fy(st.mean[j]); q ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); }); ctx.stroke();
  ctx.fillStyle = C.mean; okJ.forEach(j => { ctx.beginPath(); ctx.arc(fx(ms(j)), fy(st.mean[j]), 3, 0, 2 * Math.PI); ctx.fill(); });     // one dot per recorded sample
  ctx.restore();
  hit.wf = { fx, ms, okJ, st, Lm, pw };
  const m = waveMeasures(F, st, c), ch = F.channels[c];
  $("wfStats").innerHTML = [
    ["Electrode", `${ch.label}${F.wells[ch.well] && F.wells[ch.well].what ? " · " + esc(F.wells[ch.well].what) : ""}`],
    ["Spikes", `${fmt(nAll, 0)} in window${st.n < nAll ? ` (${fmt(st.n, 0)} drawn and averaged)` : ""}`],
    ["Mean waveform", `trough ${fmt(m.trough, 1)} µV · peak ${fmt(m.peak, 1)} µV · peak-to-peak ${fmt(m.ptp, 1)} µV`],
    ["Trough to next peak", `${fmt(m.troughToPeakMs, 2)} ms`, "wfTtp"],
    ["Noise SD (detector)", `${fmt(m.noise, 2)} µV · threshold ±${fmt(thr, 1)} µV (${fmt(F.thrMult[c], 1)} × SD)`, "wfNoise"],
    ["Signal-to-noise", `${fmt(m.snr, 1)} (peak-to-peak / noise SD)`, "wfSnr"],
  ].map(([k, v, i]) => `<div><span>${esc(k)}${i ? ` <button class="ib" data-info="${i}" aria-label="About ${esc(k)}">i</button>` : ""}</span><b>${v}</b></div>`).join("");
}

function drawWaveGrid(set) {
  const cv = $("wfGrid"), { per, gm, well, chans } = WF.cache, { ER, EC } = F.geo;
  const Wc = cv.clientWidth, g = 4, cw = (Wc - g * (EC - 1)) / EC, chh = Math.max(34, Math.min(110, cw * .75));
  cv.style.height = Math.round(ER * chh + g * (ER - 1)) + "px";
  const { ctx } = setup(cv); hit.wg = [];
  $("wfGridTitle").textContent = well ? `Well ${well}${F.wells[well] && F.wells[well].what ? " · " + F.wells[well].what : ""} · ${EC}×${ER} electrodes` : "";
  ctx.font = `12px ${C.font}`;
  for (const ch of chans) {
    const x = (ch.ec - 1) * (cw + g), y = (ch.er - 1) * (chh + g), d = per.get(ch.i), sel = ch.i === S.sel;
    ctx.lineWidth = sel ? 2 : 1; ctx.strokeStyle = sel ? C.ink : C.grid; ctx.fillStyle = C.surface;
    ctx.beginPath(); ctx.roundRect(x + .5, y + .5, cw - 1, chh - 1, 5); ctx.fill(); ctx.stroke();
    hit.wg.push({ x, y, w: cw, h: chh, c: ch.i, d });
    ctx.textBaseline = "top"; ctx.textAlign = "left"; ctx.fillStyle = sel ? C.ink : C.muted; if (cw > 28) ctx.fillText(`${ch.ec}${ch.er}`, x + 4, y + 3);
    if (!d) { if (cw > 60) { ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText("not recorded", x + cw / 2, y + chh / 2); } continue; }
    if (!d.st) { if (cw > 50) { ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText("no spikes", x + cw / 2, y + chh / 2); } continue; }
    { const lab = `n=${fmt(d.all, 0)}`; if (cw > 28 && ctx.measureText(`${ch.ec}${ch.er}`).width + ctx.measureText(lab).width + 14 <= cw) { ctx.textAlign = "right"; ctx.fillText(lab, x + cw - 4, y + 3); } } // spike count only where it fits beside the name
    const st = d.st, T = 16, B = 3, px = 3, ph = chh - T - B, pw = cw - 2 * px;
    const m = set.shared ? gm : Math.max(5, ...st.waves.map(w => Math.max(...w.map(Math.abs)))) * 1.1;
    const L = st.mean.length, fx = j => x + px + j / (L - 1) * pw, fy = v => y + T + ph / 2 - v / m * ph / 2;
    ctx.save(); ctx.beginPath(); ctx.rect(x + 1, y + T, cw - 2, ph); ctx.clip();
    ctx.strokeStyle = C["ink-2"]; ctx.lineWidth = 1; ctx.globalAlpha = Math.max(.06, Math.min(.4, 5 / st.n));
    st.waves.forEach((w, q) => { ctx.beginPath(); for (let k = 0; k < w.length; k++) { const X = fx(k - st.al[q] - st.lagMin), Y = fy(w[k]); k ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); } ctx.stroke(); });
    ctx.globalAlpha = 1; ctx.strokeStyle = C.mean; ctx.lineWidth = 2; ctx.beginPath(); let s0 = false;
    for (let j = 0; j < L; j++) { if (!st.ok[j]) continue; const X = fx(j), Y = fy(st.mean[j]); s0 ? ctx.lineTo(X, Y) : (ctx.moveTo(X, Y), s0 = true); }
    ctx.stroke(); ctx.restore();
  }
}
function fillWaveWellSelect() {
  const wells = [...F.recordedWells].sort((a, b) => a.charCodeAt(0) - b.charCodeAt(0) || parseInt(a.slice(1)) - parseInt(b.slice(1)));
  $("wfWell").innerHTML = `<option value="sel">Well of the selected electrode</option>` + wells.map(w => `<option value="${w}">${w}${F.wells[w] && F.wells[w].what ? " · " + esc(F.wells[w].what) : ""}</option>`).join("");
}
{
  for (const id of ["wfN", "wfAlign", "wfView", "wfScale", "wfWell"]) $(id).addEventListener("change", () => F && drawWaveforms());
  const grid = $("wfGrid");
  const gh = e => { const r = grid.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top; return (hit.wg || []).find(b => x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h); };
  grid.addEventListener("mousemove", e => {
    const h = gh(e); if (!h) return hideTip(); const ch = F.channels[h.c];
    if (!h.d) return showTip(e, `<b>${ch.label}</b><br>not recorded`);
    if (!h.d.st) return showTip(e, `<b>${ch.label}</b><br>no spikes in this window`);
    const m = waveMeasures(F, h.d.st, h.c);
    showTip(e, `<b>${ch.label}</b> · ${fmt(h.d.all, 0)} spikes in window${h.d.st.n < h.d.all ? ` (${h.d.st.n} drawn)` : ""}<br>mean peak-to-peak ${fmt(m.ptp, 1)} µV · SNR ${fmt(m.snr, 1)}<br>click to select`);
  });
  grid.addEventListener("mouseleave", hideTip);
  grid.addEventListener("click", e => { const h = gh(e); if (h && h.d) select(h.c); });
  $("wfMain").addEventListener("mousemove", e => {
    const h = hit.wf; if (!h) return; const x = e.clientX - $("wfMain").getBoundingClientRect().left; if (x < h.Lm || x > h.Lm + h.pw) return hideTip();
    let best = -1, bd = 1e9; for (const j of h.okJ) { const d = Math.abs(h.fx(h.ms(j)) - x); if (d < bd) { bd = d; best = j; } }
    if (best < 0) return hideTip();
    showTip(e, `<b>${fmt(h.ms(best), 3)} ms</b> (sample ${best + h.st.lagMin >= 0 ? "+" : ""}${best + h.st.lagMin})<br>mean ${fmt(h.st.mean[best], 2)} µV · SD ${fmt(h.st.sd[best], 2)} µV`);
  });
  $("wfMain").addEventListener("mouseleave", hideTip);
}

/* =====================================================================
 * Activity replay
 * Activity of an electrode at time t = number of its spikes with
 * t - Δ < spike time <= t (causal window of width Δ). A plain count:
 * no kernel, no smoothing. The plate map shows it during replay.
 * ===================================================================== */
const P = { active: false, playing: false, t: 0, speed: 1, win: 0.1, raf: 0, last: 0, counts: null, N: 1 };
function replayScale() { // fixed colour scale for the whole replay: 99th percentile of non-zero window counts over the recording
  const tmp = new Int32Array(F.channels.length), vals = [], K = 2000;
  for (let k = 1; k <= K; k++) { windowCounts(F, k / K * F.duration, P.win, tmp); for (const c of F.recorded) if (tmp[c]) vals.push(tmp[c]); }
  vals.sort((a, b) => a - b);
  P.N = Math.max(1, vals.length ? vals[Math.min(vals.length - 1, Math.floor(vals.length * .99))] : 1);
}
const fmtClock = t => { const m = Math.floor(t / 60), s = t - m * 60; return `${m}:${s.toFixed(3).padStart(6, "0")}`; };
const fmtWin = s => s >= 1 ? fmt(s, 2) + " s" : fmt(s * 1000, 0) + " ms";

function drawTimeline() { // population spike rate over the whole recording (all recorded electrodes)
  const cv = $("tl"), { ctx, W, H } = setup(cv), pw = W - RL - RR, T = 4, B = 18, ph = H - T - B;
  const nb = Math.max(20, Math.min(1200, Math.floor(pw / 2))), bw = F.duration / nb, bins = new Float64Array(nb);
  for (let i = 0; i < F.n; i++) bins[Math.min(nb - 1, Math.floor(F.times[i] / bw))]++;
  for (let b = 0; b < nb; b++) bins[b] /= bw;
  const mx = Math.max(1e-9, ...bins);
  ctx.fillStyle = C.accent; ctx.globalAlpha = .8; const bpx = pw / nb;
  for (let b = 0; b < nb; b++) if (bins[b]) { const h = bins[b] / mx * ph; ctx.fillRect(RL + b * bpx, T + ph - h, Math.max(1, bpx - (bpx > 3 ? 1 : 0)), h); }
  ctx.globalAlpha = 1; ctx.strokeStyle = C.axis; ctx.beginPath(); ctx.moveTo(RL, T + ph + .5); ctx.lineTo(RL + pw, T + ph + .5); ctx.stroke();
  ctx.fillStyle = C.muted; ctx.font = `13px ${C.font}`; ctx.textAlign = "right"; ctx.textBaseline = "top"; ctx.fillText(`${fmt(mx, 0)}/s`, RL - 6, T);
  ctx.textAlign = "center";
  for (const v of niceTicks(0, F.duration, Math.max(3, Math.floor(pw / 90)))) ctx.fillText(fmtTime(v), RL + v / F.duration * pw, T + ph + 4);
  hit.tl = { pw, T, ph, nb, bw, bins };
  drawPlayheads();
}
function drawPlayheads() {
  for (const [id, g, map] of [["tlPH", hit.tl, t => RL + t / F.duration * hit.tl.pw], ["rPH", hit.raster, t => RL + (t - S.t0) / (S.t1 - S.t0) * hit.raster.pw]]) {
    const cv = $(id), { ctx, H } = setup(cv); if (!g || !P.active) continue;
    const x = map(P.t); if (x < RL - 1 || x > RL + g.pw + 1) continue;
    const bottom = id === "rPH" ? H - 30 : g.T + g.ph;
    ctx.fillStyle = C.ink; ctx.fillRect(Math.round(x) - 1, 0, 2, bottom);
    ctx.beginPath(); ctx.moveTo(x - 5, 0); ctx.lineTo(x + 5, 0); ctx.lineTo(x, 7); ctx.fill();
  }
}
function renderReplay() {
  if (!F) return;
  if (!P.counts || P.counts.length !== F.channels.length) P.counts = new Int32Array(F.channels.length);
  const n = windowCounts(F, P.t, P.win, P.counts);
  $("clock").textContent = `${fmtClock(P.t)} / ${fmtClock(F.duration)}`;
  $("replayInfo").textContent = P.active ? `${fmt(n, 0)} spikes on the plate in the last ${fmtWin(P.win)}` : "";
  $("playBtn").textContent = P.playing ? "❚❚ Pause" : "▶ Play"; $("stopBtn").disabled = !P.active;
  // keep the playhead inside the raster window (playing or seeking): move the window, keeping its length
  if (P.active && (P.t > S.t1 || P.t < S.t0) && S.t1 - S.t0 < F.duration) { const s = S.t1 - S.t0, a = Math.max(0, Math.min(P.t, F.duration - s)); setWindow(a, a + s); }
  drawPlate(); drawPlayheads();
}
function tick(now) {
  const dt = Math.min(.25, (now - P.last) / 1000); P.last = now;
  P.t = Math.min(F.duration, P.t + dt * P.speed);
  if (P.t >= F.duration) P.playing = false;
  renderReplay();
  if (P.playing) P.raf = requestAnimationFrame(tick);
}
function play() { if (!F) return; if (P.t >= F.duration) P.t = 0; P.active = true; P.playing = true; P.last = performance.now(); cancelAnimationFrame(P.raf); P.raf = requestAnimationFrame(tick); renderReplay(); }
function pause() { P.playing = false; cancelAnimationFrame(P.raf); renderReplay(); }
function stopReplay() { P.playing = false; P.active = false; P.t = 0; cancelAnimationFrame(P.raf); if (F) renderReplay(); }
function seek(t) { P.t = Math.min(Math.max(0, t), F.duration); P.active = true; renderReplay(); }
{
  $("playBtn").addEventListener("click", e => { P.playing ? pause() : play(); e.currentTarget.blur(); });
  $("stopBtn").addEventListener("click", e => { stopReplay(); e.currentTarget.blur(); });
  $("speedSel").addEventListener("change", e => { P.speed = +e.target.value; e.target.blur(); });
  $("winSel").addEventListener("change", e => { P.win = +e.target.value; e.target.blur(); if (F) { replayScale(); renderReplay(); } });
  const tl = $("tl"); let drag = false;
  const tOf = e => Math.min(Math.max((e.clientX - tl.getBoundingClientRect().left - RL) / hit.tl.pw, 0), 1) * F.duration;
  tl.addEventListener("pointerdown", e => { if (!F || e.button !== 0) return; drag = true; tl.setPointerCapture(e.pointerId); seek(tOf(e)); });
  tl.addEventListener("pointermove", e => {
    if (!F || !hit.tl) return; const t = tOf(e);
    if (drag) { hideTip(); return seek(t); }
    const b = Math.min(hit.tl.nb - 1, Math.floor(t / hit.tl.bw));
    showTip(e, `<b>${fmtClock(t)}</b><br>${fmt(hit.tl.bins[b], 0)} spikes/s on the plate · click to move the playhead`);
  });
  tl.addEventListener("pointerup", () => { drag = false; }); tl.addEventListener("pointercancel", () => { drag = false; });
  tl.addEventListener("mouseleave", hideTip);
  document.addEventListener("keydown", e => {
    if (!F || /^(INPUT|SELECT|TEXTAREA|BUTTON)$/.test(e.target.tagName) || e.target === $("plate")) return;
    if (e.code === "Space") { e.preventDefault(); P.playing ? pause() : play(); }
    else if (e.key === "ArrowRight" || e.key === "ArrowLeft") { e.preventDefault(); seek(P.t + (e.shiftKey ? 10 : 1) * P.win * (e.key === "ArrowRight" ? 1 : -1)); }
  });
}

/* =====================================================================
 * Firing statistics panel
 * ===================================================================== */
const FS = { prm: null, a: 0, T: 1, el: [], wells: [], sort: { key: "well", dir: 1 } };
function fsParams() {
  const num = (id, def, lo) => { const v = parseFloat($(id).value); return Number.isFinite(v) && v >= lo ? v : def; };
  return { activePerMin: num("fsActive", 5, 0), fanoBin: num("fsBin", 1, 0.001), shortIsi: num("fsShort", 1.5, 0) / 1000,
    minIsi: Math.round(num("fsMinIsi", 20, 2)), acgW: num("fsAcgW", 100, 1) / 1000, acgB: num("fsAcgB", 1, 0.08) / 1000 };
}
function periodTimes(c) { // spike times of electrode c inside the analysis period, as Float64Array
  const idx = F.perCh[c], t = F.times, a = FS.a, b = FS.a + FS.T; let lo = 0, hi = idx.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (t[idx[m]] < a) lo = m + 1; else hi = m; }
  let e = lo; hi = idx.length; while (e < hi) { const m = (e + hi) >> 1; if (t[idx[m]] <= b) e = m + 1; else hi = m; }
  return Float64Array.from(idx.subarray(lo, e), i => t[i]);
}
function computeStats() {
  FS.prm = fsParams();
  if ($("fsPeriod").value === "window") { FS.a = S.t0; FS.T = S.t1 - S.t0; } else { FS.a = 0; FS.T = F.duration; }
  FS.el = [];
  for (const c of F.rows) FS.el[c] = electrodeStats(periodTimes(c), FS.a, FS.T, FS.prm);
  // per well: summaries over active electrodes
  const byWell = new Map();
  for (const c of F.rows) { const w = F.channels[c].well; if (!byWell.has(w)) byWell.set(w, []); byWell.get(w).push(c); }
  FS.wells = [...byWell.entries()].map(([w, cs]) => {
    const act = cs.filter(c => FS.el[c].active), tot = cs.reduce((s, c) => s + FS.el[c].N, 0);
    const sum = key => { const v = act.map(c => FS.el[c][key]).filter(Number.isFinite); if (!v.length) return { n: 0 };
      const mu = v.reduce((s, x) => s + x, 0) / v.length; return { n: v.length, med: quantileSorted(Float64Array.from(v).sort(), .5), mean: mu, sd: sampleSD(v, mu) }; };
    return { well: w, treat: F.wells[w] && F.wells[w].what || "", recorded: cs.length, active: act.length, wellRate: tot / FS.T,
      rate: sum("rate"), cv: sum("cv"), cv2: sum("cv2"), lv: sum("lv"), medIsi: sum("medIsi"), fano: sum("fano") };
  }).sort((p, q) => (p.treat || "~").localeCompare(q.treat || "~") || p.well.localeCompare(q.well, undefined, { numeric: true }));
}
const EL_COLS = [
  ["label", "Electrode", null], ["well", "Well", null], ["treat", "Treatment", null], ["active", "Active", "active"], ["N", "Spikes", "rate"],
  ["rate", "Rate (Hz)", "rate"], ["fano", "Fano", "fano"], ["medIsi", "Median ISI (ms)", "isi"], ["meanIsi", "Mean ISI (ms)", "isi"],
  ["cv", "CV", "cv"], ["cv2", "CV2", "cv2"], ["lv", "LV", "lv"], ["shortPct", "ISI < limit (%)", "short"]];
function elValue(c, k) {
  const ch = F.channels[c], s = FS.el[c];
  if (k === "label") return ch.label; if (k === "well") return ch.well; if (k === "treat") return F.wells[ch.well] && F.wells[ch.well].what || "";
  if (k === "active") return s.active ? 1 : 0; if (k === "medIsi" || k === "meanIsi") return s[k] * 1000; return s[k];
}
function drawStatsTable() {
  const only = $("fsOnlyActive").checked, { key, dir } = FS.sort;
  const rows = F.rows.filter(c => !only || FS.el[c].active || c === S.sel);
  const ord = c => { const ch = F.channels[c]; return ch.wr * 1e6 + ch.wc * 1e4 + ch.er * 100 + ch.ec; };
  rows.sort((p, q) => {
    if (key === "label" || key === "well") return (ord(p) - ord(q)) * dir;
    const A = elValue(p, key), B = elValue(q, key);
    if (typeof A === "string") return A.localeCompare(B) * dir || ord(p) - ord(q);
    const an = Number.isFinite(A), bn = Number.isFinite(B); if (!an || !bn) return an === bn ? ord(p) - ord(q) : an ? -1 : 1; // missing values last
    return (A - B) * dir || ord(p) - ord(q);
  });
  const d = { rate: 3, fano: 2, medIsi: 1, meanIsi: 1, cv: 2, cv2: 2, lv: 2, shortPct: 2, N: 0 };
  $("fsTable").tHead.innerHTML = "<tr>" + EL_COLS.map(([k, t, info]) => `<th data-k="${k}" class="${["label", "well", "treat"].includes(k) ? "" : "n"}">${t}${key === k ? (dir > 0 ? " ▲" : " ▼") : ""}${info ? ` <button class="ib" data-info="${info}" aria-label="About ${t}">i</button>` : ""}</th>`).join("") + "</tr>";
  $("fsTable").tBodies[0].innerHTML = rows.map(c => {
    const s = FS.el[c];
    return `<tr data-c="${c}" class="${c === S.sel ? "sel" : ""}${s.active ? "" : " inactive"}">` + EL_COLS.map(([k]) => {
      const v = elValue(c, k);
      if (k === "active") return `<td class="n">${v ? "yes" : "no"}</td>`;
      if (typeof v === "string") return `<td>${esc(v)}</td>`;
      return `<td class="n">${fmt(v, d[k])}</td>`;
    }).join("") + "</tr>";
  }).join("");
  $("fsCount").textContent = `${rows.length} electrodes · ${F.rows.filter(c => FS.el[c].active).length} active of ${F.rows.length} recorded`;
}
function drawWellTable() {
  const mode = $("fsSummary").value, cell = (o, dd) => !o.n ? "–" : mode === "median" ? fmt(o.med, dd) : `${fmt(o.mean, dd)} ± ${fmt(o.sd, dd)}`;
  const ms = o => o.n ? { ...o, med: o.med * 1000, mean: o.mean * 1000, sd: o.sd * 1000 } : o;
  $("fsWells").tHead.innerHTML = `<tr><th>Treatment</th><th>Well</th><th class="n">Active / recorded <button class="ib" data-info="active" aria-label="About active">i</button></th><th class="n">Well rate (spikes/s) <button class="ib" data-info="wellrate" aria-label="About well rate">i</button></th>`
    + [["Rate (Hz)", "rate"], ["CV", "cv"], ["CV2", "cv2"], ["LV", "lv"], ["Median ISI (ms)", "isi"], ["Fano", "fano"]].map(([t, i]) => `<th class="n">${t} <button class="ib" data-info="${i}" aria-label="About ${t}">i</button></th>`).join("") + "</tr>";
  let last = null;
  $("fsWells").tBodies[0].innerHTML = FS.wells.map(w => {
    const first = w.treat !== last; last = w.treat;
    return `<tr class="${first ? "grp" : ""}"><td>${first ? esc(w.treat || "no treatment") : ""}</td><td>${w.well}</td><td class="n">${w.active} / ${w.recorded}</td><td class="n">${fmt(w.wellRate, 2)}</td>`
      + `<td class="n">${cell(w.rate, 3)}</td><td class="n">${cell(w.cv, 2)}</td><td class="n">${cell(w.cv2, 2)}</td><td class="n">${cell(w.lv, 2)}</td><td class="n">${cell(ms(w.medIsi), 1)}</td><td class="n">${cell(w.fano, 2)}</td></tr>`;
  }).join("");
  $("fsWellNote").textContent = `${mode === "median" ? "Median" : "Mean ± SD"} over the active electrodes of each well. Electrodes in a well are not independent samples: compare treatments across wells.`;
}

/* ---------- charts of the selected electrode ---------- */
function axesFrame(ctx, W, H, o) { // o: x0,x1,y0,y1, xTicks, yTicks, xFmt, yFmt, xLabel, yLabel, logX
  const L = 52, R = 10, T = 8, B = 34, pw = W - L - R, ph = H - T - B;
  const lx = v => Math.log10(v), fx = o.logX ? v => L + (lx(v) - lx(o.x0)) / (lx(o.x1) - lx(o.x0)) * pw : v => L + (v - o.x0) / (o.x1 - o.x0) * pw;
  const fy = v => T + ph - (v - o.y0) / (o.y1 - o.y0) * ph;
  ctx.font = `13px ${C.font}`; ctx.fillStyle = C.muted; ctx.lineWidth = 1;
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const v of o.yTicks) { const y = Math.round(fy(v)) + .5; ctx.strokeStyle = C.grid; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke(); ctx.fillText(o.yFmt(v), L - 6, y); }
  ctx.strokeStyle = C.axis; ctx.beginPath(); ctx.moveTo(L, T + ph + .5); ctx.lineTo(L + pw, T + ph + .5); ctx.stroke();
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const v of o.xTicks) { const x = fx(v); if (x >= L - 1 && x <= L + pw + 1) ctx.fillText(o.xFmt(v), x, T + ph + 5); }
  ctx.fillText(o.xLabel, L + pw / 2, H - 13);
  ctx.save(); ctx.translate(11, T + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = "middle"; ctx.fillText(o.yLabel, 0, 0); ctx.restore();
  return { L, T, pw, ph, fx, fy };
}
function drawStatCharts() {
  if (folded("fsCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const c = S.sel, t = c >= 0 ? periodTimes(c) : new Float64Array(0), prm = FS.prm, s = c >= 0 ? FS.el[c] : null;
  $("fsSel").textContent = c >= 0 ? `${F.channels[c].label} · ${fmt(t.length, 0)} spikes in the period` : "";
  // 1. rate over time
  { const cv = $("fsRate"), { ctx, W, H } = setup(cv); hit.fsRate = null;
    const n = binCounts(t, FS.a, FS.T, prm.fanoBin), K = n.length;
    if (!K) { ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.fillText("Period shorter than one bin", W / 2, H / 2); }
    else {
      const r = Array.from(n, v => v / prm.fanoBin), mx = Math.max(1e-9, ...r, s.rate) * 1.05, yt = niceTicks(0, mx, 4);
      const f = axesFrame(ctx, W, H, { x0: FS.a, x1: FS.a + K * prm.fanoBin, y0: 0, y1: Math.max(mx, yt[yt.length - 1]), xTicks: niceTicks(FS.a, FS.a + K * prm.fanoBin, 6), yTicks: yt, xFmt: v => timeLabel(v, FS.T), yFmt: v => fmt(v, 1), xLabel: "time", yLabel: "spikes/s" });
      ctx.fillStyle = C.accent; const bp = f.pw / K;
      for (let k = 0; k < K; k++) if (r[k]) { const y = f.fy(r[k]); ctx.fillRect(f.L + k * bp + (bp > 3 ? .5 : 0), y, Math.max(1, bp - (bp > 3 ? 1 : 0)), f.T + f.ph - y); }
      ctx.strokeStyle = C.ink; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(f.L, f.fy(s.rate)); ctx.lineTo(f.L + f.pw, f.fy(s.rate)); ctx.stroke(); ctx.setLineDash([]);
      hit.fsRate = { f, K, bp, r };
    }
  }
  // 2. log-ISI histogram
  { const cv = $("fsIsi"), { ctx, W, H } = setup(cv); hit.fsIsi = null;
    const ih = isiHistogram(t, 10, 1e-4, 100);
    if (!ih.n) { ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.fillText("Fewer than two spikes", W / 2, H / 2); }
    else {
      const mx = Math.max(...ih.h) * 100 * 1.08 || 1, yt = niceTicks(0, mx, 4);
      const f = axesFrame(ctx, W, H, { x0: 1e-4, x1: 100, y0: 0, y1: Math.max(mx, yt[yt.length - 1]), logX: true, xTicks: [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], yTicks: yt,
        xFmt: v => v < 1e-3 ? "0.1 ms" : v < 1 ? fmt(v * 1000, 0) + " ms" : fmt(v, 0) + " s", yFmt: v => fmt(v, 1), xLabel: "inter-spike interval (log scale)", yLabel: "% of ISIs" });
      if (F.minIsiAll > 1e-4 && Number.isFinite(F.minIsiAll)) { ctx.fillStyle = C.empty; ctx.fillRect(f.L, f.T, f.fx(F.minIsiAll) - f.L, f.ph); }
      ctx.fillStyle = C.accent;
      for (let k = 0; k < ih.nb; k++) if (ih.h[k]) { const x0 = f.fx(ih.edge(k)), x1 = f.fx(ih.edge(k + 1)), y = f.fy(ih.h[k] * 100); ctx.fillRect(x0 + .5, y, Math.max(1, x1 - x0 - 1), f.T + f.ph - y); }
      ctx.strokeStyle = C.critical; ctx.setLineDash([4, 3]); ctx.lineWidth = 1.5; const xs = f.fx(prm.shortIsi);
      if (prm.shortIsi > 1e-4) { ctx.beginPath(); ctx.moveTo(xs, f.T); ctx.lineTo(xs, f.T + f.ph); ctx.stroke(); } ctx.setLineDash([]);
      hit.fsIsi = { f, ih };
    }
  }
  // 3. autocorrelogram
  { const cv = $("fsAcg"), { ctx, W, H } = setup(cv); hit.fsAcg = null;
    if (t.length < 2) { ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.fillText("Fewer than two spikes", W / 2, H / 2); }
    else {
      const a = autocorrelogram(t, prm.acgW, prm.acgB), mx = Math.max(...a.c, s.rate) * 1.08 || 1, yt = niceTicks(0, mx, 4), Wms = prm.acgW * 1000;
      const f = axesFrame(ctx, W, H, { x0: -Wms, x1: Wms, y0: 0, y1: Math.max(mx, yt[yt.length - 1]), xTicks: niceTicks(-Wms, Wms, 8), yTicks: yt, xFmt: v => fmt(v, 0), yFmt: v => fmt(v, 1), xLabel: "lag (ms)", yLabel: "spikes/s" });
      const dm = F.minIsiAll * 1000; if (Number.isFinite(dm) && dm > 0 && dm < Wms) { ctx.fillStyle = C.empty; ctx.fillRect(f.fx(-dm), f.T, f.fx(dm) - f.fx(-dm), f.ph); }
      ctx.fillStyle = C.accent; const bms = prm.acgB * 1000;
      for (let k = 0; k < 2 * a.K; k++) if (a.c[k]) { const x0 = f.fx((k - a.K) * bms), x1 = f.fx((k - a.K + 1) * bms), y = f.fy(a.c[k]); ctx.fillRect(x0 + (x1 - x0 > 3 ? .5 : 0), y, Math.max(1, x1 - x0 - (x1 - x0 > 3 ? 1 : 0)), f.T + f.ph - y); }
      ctx.strokeStyle = C.ink; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(f.L, f.fy(s.rate)); ctx.lineTo(f.L + f.pw, f.fy(s.rate)); ctx.stroke(); ctx.setLineDash([]);
      hit.fsAcg = { f, a, bms };
    }
  }
  // headline numbers of the selected electrode
  $("fsSelStats").innerHTML = !s ? "" : [["Rate", `${fmt(s.rate, 3)} Hz`, "rate"], ["Active", s.active ? "yes" : "no", "active"], ["Fano", fmt(s.fano, 2), "fano"], ["Median ISI", `${fmt(s.medIsi * 1000, 1)} ms`, "isi"],
    ["CV", fmt(s.cv, 2), "cv"], ["CV2", fmt(s.cv2, 2), "cv2"], ["LV", fmt(s.lv, 2), "lv"], [`ISI < ${fmt(prm.shortIsi * 1000, 1)} ms`, `${fmt(s.shortPct, 2)} %`, "short"]]
    .map(([k, v, i]) => `<div><span>${k} <button class="ib" data-info="${i}" aria-label="About ${k}">i</button></span><b>${v}</b></div>`).join("")
    + (s.nIsi < prm.minIsi ? `<div style="grid-column:1/-1"><span>ISI statistics need at least ${prm.minIsi} intervals; this electrode has ${s.nIsi}.</span></div>` : "");
}
function drawStats() {
  computeStats(); drawStatsTable(); drawWellTable(); drawStatCharts();
  const dm = F.minIsiAll * 1000, post = F.header.NumSamplesPerBlock - (F.n ? F.trig[0] : 0);
  $("fsDead").innerHTML = Number.isFinite(dm) ? `<b>Detector limit:</b> the shortest interval between two spikes of the same electrode in this file is ${fmt(dm, 2)} ms`
    + (F.overlap === "false" && dm >= post / F.fs * 1000 * 0.99 ? ` — just over the ${fmt(post / F.fs * 1000, 2)} ms that each cut-out lasts after the threshold crossing, and the detector settings say "Allow Sample Overlap: false". The detector cannot register a new spike before the previous cut-out ends, so intervals below ${fmt(dm, 2)} ms, the short-ISI percentage and the empty centre of the autocorrelogram reflect the detector, not the neurons.` : ".")
    + ` <button class="ib" data-info="dead" aria-label="About the detector limit">i</button>` : "";
  const p = $("fsPeriod").value === "window" ? `raster window ${timeLabel(FS.a, FS.T)} – ${timeLabel(FS.a + FS.T, FS.T)}` : `whole recording (${fmtTime(FS.T)})`;
  $("fsSub").textContent = `Per electrode over the ${p} · active: ≥ ${fmt(FS.prm.activePerMin, 2)} spikes/min · rate bins and Fano factor: ${fmt(FS.prm.fanoBin, 3)} s · ISI statistics from ≥ ${FS.prm.minIsi} intervals · click a row to select the electrode · click i for definitions`;
  if (["fano", "cv", "cv2", "lv", "medIsi", "shortPct"].includes($("plateMetric").value)) { hit.plateKey = null; drawPlate(); }
  computeNB(); drawNB();
}

/* ---------- info panels (definitions and equations) ---------- */
const M = s => `<math display="block">${s}</math>`;
const INFO = {
  rate: ["Mean firing rate", `<p>The number of spikes of the electrode in the analysis period divided by the period's length.</p>` + M(`<mi>r</mi><mo>=</mo><mfrac><mi>N</mi><mi>T</mi></mfrac>`) +
    `<p><i>N</i>: spikes with time in the period [<i>a</i>, <i>a</i> + <i>T</i>]; <i>T</i>: period length in seconds (the recording duration from the file header, or the raster window). Unit: spikes/s (Hz). An unsorted electrode can pick up several neurons, so this is a multi-unit rate.</p>`],
  active: ["Active electrode", `<p>An electrode is active when its mean firing rate reaches the threshold θ, expressed in spikes per minute (default 5 spikes/min, as in Axion's software).</p>` + M(`<mn>60</mn><mo>·</mo><mi>r</mi><mo>≥</mo><mi>θ</mi>`) +
    `<p>Well summaries use active electrodes only: silent electrodes would otherwise pull every average towards zero.</p>`],
  fano: ["Fano factor", `<p>The spike counts in consecutive bins of width Δ (default 1 s) are compared with a Poisson process, whose count variance equals its mean.</p>` + M(`<mi>F</mi><mo>=</mo><mfrac><msup><mi>s</mi><mn>2</mn></msup><mover><mi>n</mi><mo>‾</mo></mover></mfrac><mo>,</mo><mspace width="1em"/><mover><mi>n</mi><mo>‾</mo></mover><mo>=</mo><mfrac><mn>1</mn><mi>K</mi></mfrac><munderover><mo>∑</mo><mrow><mi>k</mi><mo>=</mo><mn>1</mn></mrow><mi>K</mi></munderover><msub><mi>n</mi><mi>k</mi></msub><mo>,</mo><mspace width="1em"/><msup><mi>s</mi><mn>2</mn></msup><mo>=</mo><mfrac><mn>1</mn><mrow><mi>K</mi><mo>−</mo><mn>1</mn></mrow></mfrac><munderover><mo>∑</mo><mrow><mi>k</mi><mo>=</mo><mn>1</mn></mrow><mi>K</mi></munderover><msup><mrow><mo>(</mo><msub><mi>n</mi><mi>k</mi></msub><mo>−</mo><mover><mi>n</mi><mo>‾</mo></mover><mo>)</mo></mrow><mn>2</mn></msup>`) +
    `<p><i>n<sub>k</sub></i>: spikes in bin <i>k</i> = [<i>a</i> + (<i>k</i>−1)Δ, <i>a</i> + <i>k</i>Δ); <i>K</i> = ⌊<i>T</i>/Δ⌋ whole bins (a final partial bin is left out). F = 1 for Poisson firing, F &gt; 1 for bursty or drifting activity, F &lt; 1 for regular firing. F depends on Δ: compare values computed with the same bin width.</p>`],
  isi: ["Inter-spike intervals: mean and median", `<p>The intervals between consecutive spikes of the same electrode inside the period:</p>` + M(`<msub><mi>τ</mi><mi>i</mi></msub><mo>=</mo><msub><mi>t</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>t</mi><mi>i</mi></msub><mo>,</mo><mspace width="1em"/><mi>i</mi><mo>=</mo><mn>1</mn><mo>…</mo><mi>n</mi><mo>,</mo><mspace width="1em"/><mi>n</mi><mo>=</mo><mi>N</mi><mo>−</mo><mn>1</mn>`) + M(`<mover><mi>τ</mi><mo>‾</mo></mover><mo>=</mo><mfrac><mn>1</mn><mi>n</mi></mfrac><munderover><mo>∑</mo><mrow><mi>i</mi><mo>=</mo><mn>1</mn></mrow><mi>n</mi></munderover><msub><mi>τ</mi><mi>i</mi></msub>`) +
    `<p>The median is the middle interval (the mean of the two middle ones when <i>n</i> is even). Spike times have the resolution of one sample (80 µs at 12.5 kHz). ISI statistics are reported only when <i>n</i> ≥ the minimum number of intervals (default 20); otherwise they show "–".</p>`],
  cv: ["Coefficient of variation (CV)", `<p>The spread of the intervals relative to their mean.</p>` + M(`<mi>CV</mi><mo>=</mo><mfrac><msub><mi>s</mi><mi>τ</mi></msub><mover><mi>τ</mi><mo>‾</mo></mover></mfrac><mo>,</mo><mspace width="1em"/><msub><mi>s</mi><mi>τ</mi></msub><mo>=</mo><msqrt><mfrac><mn>1</mn><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></mfrac><munderover><mo>∑</mo><mrow><mi>i</mi><mo>=</mo><mn>1</mn></mrow><mi>n</mi></munderover><msup><mrow><mo>(</mo><msub><mi>τ</mi><mi>i</mi></msub><mo>−</mo><mover><mi>τ</mi><mo>‾</mo></mover><mo>)</mo></mrow><mn>2</mn></msup></msqrt>`) +
    `<p>CV = 1 for a Poisson process, CV &lt; 1 for regular firing, CV &gt; 1 for bursting. Slow changes of the rate over the recording also raise the CV, even without bursts: CV2 and LV are insensitive to them.</p>`],
  cv2: ["CV2 (Holt et al., 1996)", `<p>A local measure: each pair of neighbouring intervals is compared, then the values are averaged.</p>` + M(`<msub><mi>CV</mi><mn>2</mn></msub><mo>=</mo><mfrac><mn>1</mn><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></mfrac><munderover><mo>∑</mo><mrow><mi>i</mi><mo>=</mo><mn>1</mn></mrow><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></munderover><mfrac><mrow><mn>2</mn><mo>|</mo><msub><mi>τ</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>τ</mi><mi>i</mi></msub><mo>|</mo></mrow><mrow><msub><mi>τ</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>+</mo><msub><mi>τ</mi><mi>i</mi></msub></mrow></mfrac>`) +
    `<p>Ranges from 0 to 2; about 1 for Poisson firing, lower for regular, higher for irregular or bursting firing. Because only neighbouring intervals are compared, slow rate changes do not inflate it. Holt GR, Softky WR, Koch C, Douglas RJ (1996) J Neurophysiol 75:1806–1814.</p>`],
  lv: ["Local variation (LV, Shinomoto et al., 2003)", `<p>A local measure of irregularity, standard for classifying firing patterns.</p>` + M(`<mi>LV</mi><mo>=</mo><mfrac><mn>3</mn><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></mfrac><munderover><mo>∑</mo><mrow><mi>i</mi><mo>=</mo><mn>1</mn></mrow><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></munderover><msup><mrow><mo>(</mo><mfrac><mrow><msub><mi>τ</mi><mi>i</mi></msub><mo>−</mo><msub><mi>τ</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub></mrow><mrow><msub><mi>τ</mi><mi>i</mi></msub><mo>+</mo><msub><mi>τ</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub></mrow></mfrac><mo>)</mo></mrow><mn>2</mn></msup>`) +
    `<p>LV = 1 for a Poisson process, LV &lt; 1 for regular firing, LV &gt; 1 for bursting; like CV2, it is insensitive to slow rate changes. Shinomoto S, Shima K, Tanji J (2003) Neural Comput 15:2823–2842.</p>`],
  short: ["Intervals shorter than the limit", `<p>The percentage of intervals shorter than the limit <i>L</i> (default 1.5 ms).</p>` + M(`<mn>100</mn><mo>·</mo><mfrac><mrow><mo>#</mo><mo>{</mo><mi>i</mi><mo>:</mo><msub><mi>τ</mi><mi>i</mi></msub><mo>&lt;</mo><mi>L</mi><mo>}</mo></mrow><mi>n</mi></mfrac>`) +
    `<p>A single neuron cannot fire again within its refractory period (about 1–2 ms). On an unsorted electrode, short intervals are expected when several neurons are recorded. After spike sorting, the same number becomes a quality measure of each unit (refractory-period violations).</p>`],
  ratebin: ["Rate bin Δ", `<p>The analysis period [<i>a</i>, <i>a</i> + <i>T</i>] is cut into <i>K</i> whole bins of width Δ (default 1 s); the remainder shorter than Δ at the end is left out:</p>` + M(`<mi>K</mi><mo>=</mo><mo>⌊</mo><mfrac><mi>T</mi><mi>Δ</mi></mfrac><mo>⌋</mo>`) + M(`<msub><mi>n</mi><mi>k</mi></msub><mo>=</mo><mo>#</mo><mo>{</mo><mi>t</mi><mo>:</mo><mi>a</mi><mo>+</mo><mi>k</mi><mi>Δ</mi><mo>≤</mo><mi>t</mi><mo>&lt;</mo><mi>a</mi><mo>+</mo><mo>(</mo><mi>k</mi><mo>+</mo><mn>1</mn><mo>)</mo><mi>Δ</mi><mo>}</mo><mo>,</mo><mspace width="1em"/><mi>k</mi><mo>=</mo><mn>0</mn><mo>…</mo><mi>K</mi><mo>−</mo><mn>1</mn>`) +
    `<p>The same counts <i>n<sub>k</sub></i> are used twice: the <b>Rate over time</b> chart plots <i>n<sub>k</sub></i>/Δ, and the <b>Fano factor</b> is their variance over their mean. Changing Δ changes both; the mean rate, ISI statistics and autocorrelogram do not depend on it. The Fano factor grows with Δ when the firing is bursty, so values are comparable only at the same Δ. See the <b>i</b> buttons of Rate over time and Fano for their formulas.</p>`],
  minisi: ["Minimum number of intervals", `<p>The interval statistics of an electrode are computed only when it has enough inter-spike intervals in the period:</p>` + M(`<mi>n</mi><mo>=</mo><mi>N</mi><mo>−</mo><mn>1</mn><mo>≥</mo><msub><mi>n</mi><mtext>min</mtext></msub>`) +
    `<p><i>N</i>: spikes of the electrode in the period; <i>n</i><sub>min</sub>: this setting (default 20, at least 2). Below it, mean ISI, median ISI, CV, CV2, LV and the % of short intervals are left empty (–) rather than estimated from a handful of intervals, and the electrode is left out of those columns in the well summaries. Spike count, rate, active flag and Fano factor do not depend on it, and the ISI histogram and autocorrelogram of the selected electrode are drawn whatever the number of intervals.</p>`],
  ratebins: ["Rate over time", `<p>The spike count in consecutive bins of width Δ divided by Δ: a plain histogram, no smoothing.</p>` + M(`<msub><mi>r</mi><mi>k</mi></msub><mo>=</mo><mfrac><msub><mi>n</mi><mi>k</mi></msub><mi>Δ</mi></mfrac>`) +
    `<p>The same bins are used for the Fano factor. The dashed line is the mean rate <i>r</i>.</p>`],
  isihist: ["Inter-spike interval histogram", `<p>The fraction of intervals in each bin of a logarithmic axis, 10 bins per decade from 0.1 ms to 100 s; bin <i>k</i> spans</p>` + M(`<mo>[</mo><msup><mn>10</mn><mrow><mo>−</mo><mn>4</mn><mo>+</mo><mi>k</mi><mo>/</mo><mn>10</mn></mrow></msup><mo>,</mo><mspace width="0.3em"/><msup><mn>10</mn><mrow><mo>−</mo><mn>4</mn><mo>+</mo><mo>(</mo><mi>k</mi><mo>+</mo><mn>1</mn><mo>)</mo><mo>/</mo><mn>10</mn></mrow></msup><mo>)</mo><mspace width="0.3em"/><mtext>s</mtext>`) +
    `<p>A log axis shows intervals over many orders of magnitude. Bursting electrodes give two peaks: short intervals inside bursts and long ones between them. The red dashed line marks the short-interval limit.</p>`],
  acg: ["Autocorrelogram", `<p>For every spike, the other spikes of the same electrode within ±<i>W</i> (default 100 ms) are counted by their lag, in bins of width <i>b</i> (default 1 ms), and normalised to a rate:</p>` + M(`<mi>A</mi><mo>(</mo><mi>k</mi><mo>)</mo><mo>=</mo><mfrac><mrow><mo>#</mo><mo>{</mo><mo>(</mo><mi>i</mi><mo>,</mo><mi>j</mi><mo>)</mo><mo>:</mo><mi>i</mi><mo>≠</mo><mi>j</mi><mo>,</mo><mspace width="0.2em"/><mi>k</mi><mi>b</mi><mo>≤</mo><mo>|</mo><msub><mi>t</mi><mi>j</mi></msub><mo>−</mo><msub><mi>t</mi><mi>i</mi></msub><mo>|</mo><mo>&lt;</mo><mo>(</mo><mi>k</mi><mo>+</mo><mn>1</mn><mo>)</mo><mi>b</mi><mo>}</mo></mrow><mrow><mn>2</mn><mo>·</mo><mi>N</mi><mo>·</mo><mi>b</mi></mrow></mfrac>`) +
    `<p>Each pair is counted once at +|lag| and once at −|lag|, so the plot is symmetric; a spike is never paired with itself. The value is the rate of the electrode at a given lag after (or before) one of its spikes; far from 0 it approaches the mean rate (dashed line). A gap near 0 shows refractoriness, a peak at short lags shows bursting, side peaks show rhythmic firing.</p>`],
  dead: ["Detector limit (shortest measurable interval)", `<p>The page finds the shortest interval between two consecutive spikes of the same electrode over the whole file:</p>` + M(`<msub><mi>τ</mi><mtext>min</mtext></msub><mo>=</mo><munder><mo>min</mo><mrow><mi>e</mi><mo>,</mo><mi>i</mi></mrow></munder><mo>(</mo><msub><mi>t</mi><mrow><mi>e</mi><mo>,</mo><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>t</mi><mrow><mi>e</mi><mo>,</mo><mi>i</mi></mrow></msub><mo>)</mo>`) +
    `<p>When the AxIS spike detector runs with "Allow Sample Overlap: false", it does not detect a new spike on an electrode until the cut-out of the previous one has ended (the post-spike duration after the threshold crossing, 2.16 ms in these files). Intervals shorter than that are never recorded, whatever the neurons do. Below τ<sub>min</sub> the ISI histogram and the autocorrelogram are therefore empty by construction (shaded), and the percentage of short intervals is 0. Refractoriness of single neurons can only be judged above this limit, and only after spike sorting.</p>`],
  wfTtp: ["Trough to next peak", `<p>Measured on the mean waveform of the selected electrode: the average, sample by sample, of the spikes drawn in the panel (those in the raster window, at most the number chosen in the panel, taken evenly through the window), after the alignment chosen in the panel. Only lags covered by at least half of the spikes are used.</p>` + M(`<msub><mi>j</mi><mtext>tr</mtext></msub><mo>=</mo><munder><mrow><mi>arg</mi><mo>⁡</mo><mi>min</mi></mrow><mi>j</mi></munder><mspace width="0.2em"/><mover><mi>w</mi><mo>‾</mo></mover><mo>(</mo><mi>j</mi><mo>)</mo><mo>,</mo><mspace width="1em"/><msub><mi>j</mi><mtext>pk</mtext></msub><mo>=</mo><munder><mrow><mi>arg</mi><mo>⁡</mo><mi>max</mi></mrow><mrow><mi>j</mi><mo>&gt;</mo><msub><mi>j</mi><mtext>tr</mtext></msub></mrow></munder><mspace width="0.2em"/><mover><mi>w</mi><mo>‾</mo></mover><mo>(</mo><mi>j</mi><mo>)</mo><mo>,</mo><mspace width="1em"/><msub><mi>t</mi><mtext>tp</mtext></msub><mo>=</mo><mfrac><mrow><msub><mi>j</mi><mtext>pk</mtext></msub><mo>−</mo><msub><mi>j</mi><mtext>tr</mtext></msub></mrow><msub><mi>f</mi><mi>s</mi></msub></mfrac>`) +
    `<p><i>j</i>: sample index; <i>f<sub>s</sub></i>: sampling frequency. The next peak is the largest value of the mean after the trough, within the stored cut-out. The time is a whole number of samples (1/<i>f<sub>s</sub></i> = 0.08 ms at 12.5 kHz), with no interpolation. No value is shown when the trough is the last usable sample. On an unsorted electrode the mean can mix several neurons, so this is a property of the average waveform, not of one cell.</p>`],
  wfNoise: ["Noise SD (detector)", `<p>The noise level is not computed by this page. The AxIS spike detector stores its own estimate of the electrode's noise standard deviation σ with every spike it detects (the field read from the file, in volts). The page shows the median over all the electrode's spikes in the whole file, in µV:</p>` + M(`<mover><mi>σ</mi><mo>~</mo></mover><mo>=</mo><mtext>median</mtext><mo>(</mo><msub><mi>σ</mi><mn>1</mn></msub><mo>,</mo><mo>…</mo><mo>,</mo><msub><mi>σ</mi><mi>N</mi></msub><mo>)</mo><mo>,</mo><mspace width="1em"/><mtext>threshold</mtext><mo>=</mo><mo>±</mo><mi>k</mi><mo>·</mo><mover><mi>σ</mi><mo>~</mo></mover>`) +
    `<p><i>k</i>: the threshold multiplier stored with the spikes (the detector setting, e.g. 6 × SD). The threshold is drawn as the dashed red lines. The detector's method is listed under recording details (Spike Detector Settings). With an adaptive threshold the estimate changes during the recording (on the test file's electrode B2_15, from 1.8 to 3.5 µV, median 2.1 µV), so the median is shown rather than one spike's value.</p>`],
  wfSnr: ["Signal-to-noise", `<p>The peak-to-peak amplitude of the mean waveform divided by the detector's noise SD:</p>` + M(`<mi>SNR</mi><mo>=</mo><mfrac><mrow><munder><mi>max</mi><mi>j</mi></munder><mover><mi>w</mi><mo>‾</mo></mover><mo>(</mo><mi>j</mi><mo>)</mo><mo>−</mo><munder><mi>min</mi><mi>j</mi></munder><mover><mi>w</mi><mo>‾</mo></mover><mo>(</mo><mi>j</mi><mo>)</mo></mrow><mover><mi>σ</mi><mo>~</mo></mover></mfrac>`) +
    `<p>w̄: the mean waveform described under Trough to next peak (same spikes, alignment and usable lags); σ̃: the median noise SD of the detector (see Noise SD). The maximum and minimum are taken over the whole mean, in either order. Dimensionless. Averaging removes much of the noise from w̄, so a single spike usually stands out less than this ratio suggests. Other software defines SNR differently (for example the trough alone, or a different noise estimate), so values are only comparable within this page.</p>`],
  wellrate: ["Well rate", `<p>All spikes of the well's recorded electrodes in the period divided by its length.</p>` + M(`<msub><mi>R</mi><mtext>well</mtext></msub><mo>=</mo><mfrac><mrow><msub><mo>∑</mo><mi>e</mi></msub><msub><mi>N</mi><mi>e</mi></msub></mrow><mi>T</mi></mfrac>`) +
    `<p>The other columns summarise the per-electrode values over the well's active electrodes, as the median or as mean ± sample SD.</p>`],
  nbCommon: ["Network bursts: common rules", `<p>Each well is analysed on its own. The input is the pooled spike train of the well's <b>active</b> electrodes (the threshold θ of Firing statistics), over the same period as Firing statistics:</p>` + M(`<msub><mi>s</mi><mn>1</mn></msub><mo>≤</mo><msub><mi>s</mi><mn>2</mn></msub><mo>≤</mo><mo>…</mo><mo>≤</mo><msub><mi>s</mi><mi>M</mi></msub><mo>,</mo><mspace width="1em"/><msub><mi>e</mi><mi>j</mi></msub><mo>=</mo><mtext>electrode of </mtext><msub><mi>s</mi><mi>j</mi></msub>`) +
    `<p>Spikes with the same time are ordered by electrode. <i>A</i> is the number of active electrodes of the well, and every method requires at least <i>m</i> electrodes (default 3; <b>Electrodes ≥ <i>m</i></b> in the All methods box). Methods 2 and 3 also require a fraction <i>f</i> of the active electrodes:</p>` + M(`<msub><mi>n</mi><mtext>req</mtext></msub><mo>=</mo><mo>max</mo><mo>(</mo><mi>m</mi><mo>,</mo><mo>⌈</mo><mi>f</mi><mo>·</mo><mi>A</mi><mo>⌉</mo><mo>)</mo>`) +
    `<p>After detection, the network bursts of a method whose intervals overlap or touch (the start of one ≤ the end of the previous) are merged into one, which lasts from the earliest start to the latest end. A merged network burst is not checked again against the method's criteria. Merging matters mostly for burst grouping, whose network bursts can overlap when a long single-electrode burst runs past the next group's onset. Each network burst is the interval [start, end]; its spikes and electrodes are those of the pooled train with start ≤ <i>s<sub>j</sub></i> ≤ end. No smoothing is applied anywhere.</p>`],
  nbIsin: ["Method 1 · ISI<sub>N</sub> (Bakkum et al., 2014)", `<p>For every run of <i>N</i> consecutive spikes of the pooled train (default <i>N</i> = 10), the time it spans:</p>` + M(`<msub><mi>ISI</mi><mi>N</mi></msub><mo>(</mo><mi>i</mi><mo>)</mo><mo>=</mo><msub><mi>s</mi><mrow><mi>i</mi><mo>+</mo><mi>N</mi><mo>−</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>s</mi><mi>i</mi></msub>`) +
    `<p>Window <i>i</i> (spikes <i>i</i> … <i>i</i>+<i>N</i>−1) qualifies when ISI<sub>N</sub>(<i>i</i>) ≤ θ<sub>N</sub>. Qualifying windows that share at least one spike are joined; the spikes they cover form a candidate, kept when it has at least the minimum number of spikes (default 10) on at least <i>m</i> distinct electrodes (<i>m</i>: <b>Electrodes ≥ <i>m</i></b> in the All methods box, shared by the three methods; default 3). Start and end are its first and last spike. Only the minimum <i>m</i> applies here, not the fraction <i>f</i>. Every candidate already contains at least <i>N</i> spikes, so the minimum number of spikes has an effect only when it is larger than <i>N</i> (with the defaults, <i>N</i> = 10 and minimum 10, it has none). With fewer than <i>N</i> spikes in the well there are no windows and no network bursts.</p>` +
    `<p>θ<sub>N</sub> should sit in the valley between the short ISI<sub>N</sub> inside bursts and the long ones between them (histogram below). <b>Valley</b> sets it to the upper edge of the deepest bin after the fullest one, on the plain histogram: for each bin <i>v</i> after the fullest bin <i>p</i>, depth(<i>v</i>) = min(<i>h<sub>p</sub></i>, max<sub><i>k</i>&gt;<i>v</i></sub> <i>h<sub>k</sub></i>) − <i>h<sub>v</sub></i>. Only bins from <i>p</i>+1 to the second-to-last are candidates; when two bins are equally deep, the earlier one is taken; if no bin has a positive depth, no threshold is proposed and θ<sub>N</sub> is left unchanged. The rule of ≥ 10 spikes on ≥ 3 electrodes is the one of the MEA-NAP pipeline. Bakkum DJ et al. (2014) Parameters for burst detection. Front Comput Neurosci 7:193.</p>`],
  nbGroup: ["Method 2 · Burst grouping (after Chiappalone et al., 2006; Pasquale et al., 2010)", `<p>First, bursts of single electrodes: maximal runs of consecutive spikes of one electrode, <i>t<sub>j</sub></i> … <i>t<sub>l</sub></i>, with</p>` + M(`<msub><mi>t</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>t</mi><mi>i</mi></msub><mo>≤</mo><msub><mi>ISI</mi><mtext>max</mtext></msub><mspace width="0.6em"/><mtext>for</mtext><mspace width="0.4em"/><mi>j</mi><mo>≤</mo><mi>i</mi><mo>&lt;</mo><mi>l</mi><mo>,</mo><mspace width="1em"/><mi>l</mi><mo>−</mo><mi>j</mi><mo>+</mo><mn>1</mn><mo>≥</mo><msub><mi>n</mi><mi>b</mi></msub>`) +
    `<p>(defaults ISI<sub>max</sub> = 100 ms, <i>n<sub>b</sub></i> = 5 spikes). Only the spikes in the period are used, so a burst that crosses the edge of the period is cut there. The onsets of all the well's bursts are sorted by time (ties by electrode). Groups are formed by chaining: each onset joins the current group when it is at most τ after the previous onset (default 100 ms), otherwise it starts a new group. A group can therefore last much longer than τ. A group is a network burst when at least <i>n</i><sub>req</sub> distinct electrodes have a burst in it (default <i>f</i> = 20 % of the active electrodes). It starts at the first onset and ends at the latest end of its bursts.</p>` +
    `<p>Because whole single-electrode bursts are joined, a network burst made of several reverberations comes out as one long event. Chiappalone M et al. (2006) Brain Res 1093:41–53; Pasquale V, Martinoia S, Chiappalone M (2010) J Comput Neurosci 29:213–229. Their adaptive (logISI) thresholds are not used here: the thresholds are fixed and shown.</p>`],
  nbRecruit: ["Method 3 · Recruitment in time bins", `<p>The period is cut into <i>K</i> = ⌊<i>T</i>/<i>b</i>⌋ whole bins of width <i>b</i> (default 50 ms). For each bin, the number of active electrodes that fire at least once in it:</p>` + M(`<msub><mi>R</mi><mi>k</mi></msub><mo>=</mo><mo>#</mo><mo>{</mo><mi>e</mi><mo>:</mo><mo>∃</mo><mi>t</mi><mo>∈</mo><msub><mi>T</mi><mi>e</mi></msub><mo>,</mo><mspace width="0.3em"/><mi>a</mi><mo>+</mo><mi>k</mi><mi>b</mi><mo>≤</mo><mi>t</mi><mo>&lt;</mo><mi>a</mi><mo>+</mo><mo>(</mo><mi>k</mi><mo>+</mo><mn>1</mn><mo>)</mo><mi>b</mi><mo>}</mo>`) +
    `<p>Bin <i>k</i> qualifies when <i>R<sub>k</sub></i> ≥ <i>n</i><sub>req</sub> (default <i>f</i> = 25 % of the active electrodes). A run of consecutive qualifying bins is a network burst, from its first to its last spike, so it can be shorter than one bin. Bins start at the beginning <i>a</i> of the period; spikes in the last, incomplete bin (the remainder <i>T</i> − <i>K b</i>) are not counted. Counting electrodes rather than spikes keeps one fast-firing electrode from creating network bursts by itself. No published paper defines this exact rule. It belongs to the family of detectors that call a network burst when many channels are active at the same time (reviewed by Välkki et al., 2017, who cite e.g. van Pelt et al., 2004 and Raichman &amp; Ben-Jacob, 2008); those detectors mostly threshold the total firing of the array, whereas this one counts electrodes per bin. The bin width <i>b</i> = 50 ms and the fraction <i>f</i> = 25 % are this pipeline's own defaults, with no published reference; change them to suit the recordings. The bars above the raster show <i>R<sub>k</sub></i>, the dashed line <i>n</i><sub>req</sub>.</p><p>Välkki IA, Lenk K, Mikkonen JE, Kapucu FE, Hyttinen JAK (2017) Network-wide adaptive burst detection depicts neuronal activity with improved accuracy. Front Comput Neurosci 11:40.</p>`],
  nbHist: ["ISI<sub>N</sub> histogram", `<p>The fraction of windows of <i>N</i> spikes of the pooled train per bin of a log axis, 10 bins per decade from 0.1 ms to 100 s, without smoothing:</p>` + M(`<msub><mi>ISI</mi><mi>N</mi></msub><mo>(</mo><mi>i</mi><mo>)</mo><mo>=</mo><msub><mi>s</mi><mrow><mi>i</mi><mo>+</mo><mi>N</mi><mo>−</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mi>s</mi><mi>i</mi></msub><mo>,</mo><mspace width="1em"/><mi>i</mi><mo>=</mo><mn>1</mn><mo>…</mo><mi>M</mi><mo>−</mo><mi>N</mi><mo>+</mo><mn>1</mn>`) +
    `<p>With network bursts the histogram has two modes: windows inside bursts (short) and windows that include quiet time (long). The dashed line is θ<sub>N</sub>.</p>` +
    `<p>Bin <i>k</i> covers [0.1 ms · 10<sup><i>k</i>/10</sup>, 0.1 ms · 10<sup>(<i>k</i>+1)/10</sup>): the lower edge is included, the upper one is not. The denominator is all <i>M</i> − <i>N</i> + 1 windows; windows that span 0 s (<i>N</i> spikes at the same time) or fall outside 0.1 ms–100 s have no bar but still count in it.</p>`],
  nbAgree: ["Agreement between methods", `<p>For the network bursts <i>X</i> of one method and <i>Y</i> of another, the percentage of <i>X</i> that overlap at least one member of <i>Y</i>:</p>` + M(`<mn>100</mn><mo>·</mo><mfrac><mrow><mo>#</mo><mo>{</mo><mi>x</mi><mo>∈</mo><mi>X</mi><mo>:</mo><mo>∃</mo><mi>y</mi><mo>∈</mo><mi>Y</mi><mo>,</mo><mspace width="0.3em"/><msub><mi>y</mi><mtext>start</mtext></msub><mo>≤</mo><msub><mi>x</mi><mtext>end</mtext></msub><mo>∧</mo><msub><mi>y</mi><mtext>end</mtext></msub><mo>≥</mo><msub><mi>x</mi><mtext>start</mtext></msub><mo>}</mo></mrow><mrow><mo>|</mo><mi>X</mi><mo>|</mo></mrow></mfrac>`) +
    `<p>It is not symmetric. If one method splits a long network burst into several reverberations and another keeps it whole, the long one overlaps the short ones (100 %), but short events detected outside the long ones lower the other direction.</p>`],
  cmp: ["Multi-recording analysis", `<p>Every recording is analysed in its own tab, named after the file. The comparison uses each tab's results as they are: the parameters, period and active threshold set in that tab. A file added here that is not open yet gets its own tab, analysed with the default parameters; a file that is already open (same name, size and date) is not read again. Changed parameters are listed under the recording.</p>` +
    `<p>For each network-burst method (ISI<sub>N</sub>, burst grouping, recruitment) and each recording, over the analysis period of length <i>T</i> of that recording, with the <i>n</i> network bursts [start<sub><i>q</i></sub>, end<sub><i>q</i></sub>] found in one well, or in all wells together (each well analysed on its own, as in Network bursts):</p>` +
    M(`<mi>n</mi><mo>=</mo><munder><mo>∑</mo><mi>w</mi></munder><msub><mi>n</mi><mi>w</mi></msub><mo>,</mo><mspace width="1em"/><mtext>rate</mtext><mo>=</mo><mfrac><mrow><mn>60</mn><mi>n</mi></mrow><mi>T</mi></mfrac><mtext>&#160;per min</mtext><mo>,</mo><mspace width="1em"/><msub><mi>d</mi><mi>q</mi></msub><mo>=</mo><msub><mtext>end</mtext><mi>q</mi></msub><mo>−</mo><msub><mtext>start</mtext><mi>q</mi></msub>`) +
    `<p>Count: <i>n</i>, not divided by the duration. Rate: <i>n</i> per minute of the period (for all wells, the same as adding the wells' rates). Median duration: the median of <i>d<sub>q</sub></i> over the network bursts of the well, or over those of all wells pooled; the median interpolates between the two middle values when <i>n</i> is even, and there is no value when <i>n</i> = 0. A recording without the chosen well has no values.</p>` +
    `<p>The spider plot has one axis per measure and method (up to 9: the ticked measures × 3 methods) and one outline per selected recording. The centre is 0 and each radius is linear. With <i>Each axis to its own maximum</i>, the outer ring of an axis is the smallest round number (1, 2, 2.5 or 5 × 10<sup><i>k</i></sup>) at or above the largest value on that axis among the selected recordings; with <i>One scale per measure</i>, the three methods of a measure share the largest of their three. The range of each axis is written at its end. A missing value has no point, and the outline joins the points that exist. When all recordings have the same period, the count and rate axes have the same shape. Colours follow the recording, in the order it was selected, up to 8 at once.</p>`],
  cmpFiring: ["Firing per recording", `<p>For each selected recording, over the Firing statistics period [<i>a</i>, <i>a</i> + <i>T</i>] of its tab, in the chosen well or over all its recorded wells (electrodes pooled):</p>` + M(`<msub><mi>R</mi><mtext>net</mtext></msub><mo>=</mo><mfrac><mrow><msub><mo>∑</mo><mrow><mi>e</mi><mo>∈</mo><mtext>recorded</mtext></mrow></msub><msub><mi>N</mi><mi>e</mi></msub></mrow><mi>T</mi></mfrac>`) + M(`<mover><mi>r</mi><mo>‾</mo></mover><mo>=</mo><mfrac><mn>1</mn><mi>A</mi></mfrac><munder><mo>∑</mo><mrow><mi>e</mi><mo>∈</mo><mtext>active</mtext></mrow></munder><msub><mi>r</mi><mi>e</mi></msub><mo>,</mo><mspace width="1em"/><msub><mi>r</mi><mi>e</mi></msub><mo>=</mo><mfrac><msub><mi>N</mi><mi>e</mi></msub><mi>T</mi></mfrac>`) +
    `<p><b>Network firing rate</b> <i>R</i><sub>net</sub>: every spike of the recorded electrodes divided by the period, in spikes/s (the "well rate" of Firing statistics; over all wells, the sum). <b>Mean electrode rate</b> <i>r̄</i>: the mean of the firing rates of the <i>A</i> <b>active</b> electrodes (rate ≥ the threshold θ of that tab, default 5 spikes/min), in Hz; the table also gives their sample SD. <b>Active electrodes</b>: <i>A</i>, out of the recorded ones. Silent electrodes count in <i>R</i><sub>net</sub> (they add nothing) but not in <i>r̄</i>. These are spikes of unsorted electrodes, so each rate can mix several neurons.</p>` +
    `<p>The spider plot scales each axis from 0 to the round number at or above its largest value among the selected recordings. Comparisons are meaningful only with the same period type and active threshold in every tab; the threshold of each recording is listed in the table.</p>`],
  nbStats: ["Network burst statistics", `<p>Per well and method, over the <i>n</i> network bursts found in the period of length <i>T</i>:</p>` + M(`<mtext>rate</mtext><mo>=</mo><mfrac><mrow><mn>60</mn><mi>n</mi></mrow><mi>T</mi></mfrac><mspace width="0.4em"/><mtext>per min</mtext><mo>,</mo><mspace width="1em"/><msub><mi>d</mi><mi>q</mi></msub><mo>=</mo><msub><mtext>end</mtext><mi>q</mi></msub><mo>−</mo><msub><mtext>start</mtext><mi>q</mi></msub><mo>,</mo><mspace width="1em"/><msub><mi>I</mi><mi>q</mi></msub><mo>=</mo><msub><mtext>start</mtext><mrow><mi>q</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>−</mo><msub><mtext>start</mtext><mi>q</mi></msub>`) +
    `<p>Duration <i>d</i>: median and mean. Interval <i>I</i> between onsets: mean, and CV (sample SD / mean). Spikes and electrodes per network burst: median. "Spikes in NBs": the percentage of the pooled spikes of the well's active electrodes that fall inside network bursts.</p>` +
    `<p>Medians interpolate linearly between the two middle values when the count is even (numpy's default). The mean interval needs at least 2 network bursts, and its CV at least 3 (2 intervals, SD with <i>n</i> − 1). A dash means the value is not defined; an agreement of a method with no network bursts is shown as a dash too.</p>`],
};
function showInfo(key) {
  const [title, html] = INFO[key]; $("infoTitle").innerHTML = title; $("infoBody").innerHTML = html; hideTip(); $("infoDlg").showModal();
}
{
  document.addEventListener("click", e => { const b = e.target.closest(".ib"); if (b) { e.preventDefault(); e.stopPropagation(); showInfo(b.dataset.info); } }, true);
  $("infoClose").addEventListener("click", () => $("infoDlg").close());
  $("infoDlg").addEventListener("click", e => { if (e.target === $("infoDlg")) $("infoDlg").close(); });
  for (const id of ["fsActive", "fsBin", "fsShort", "fsMinIsi", "fsAcgW", "fsAcgB", "fsPeriod"]) $(id).addEventListener("change", () => F && drawStats());
  $("fsSummary").addEventListener("change", () => F && drawWellTable());
  $("fsOnlyActive").addEventListener("change", () => F && drawStatsTable());
  $("fsTable").addEventListener("click", e => {
    const th = e.target.closest("th[data-k]"); if (th) { const k = th.dataset.k; FS.sort = { key: k, dir: FS.sort.key === k ? -FS.sort.dir : (["label", "well", "treat"].includes(k) ? 1 : -1) }; return drawStatsTable(); }
    const tr = e.target.closest("tr[data-c]"); if (tr) select(+tr.dataset.c);
  });
  const hov = (id, fn) => { $(id).addEventListener("mousemove", e => { const r = e.target.getBoundingClientRect(), html = fn(e.clientX - r.left, e.clientY - r.top); html ? showTip(e, html) : hideTip(); }); $(id).addEventListener("mouseleave", hideTip); };
  hov("fsRate", x => { const h = hit.fsRate; if (!h) return null; const k = Math.floor((x - h.f.L) / h.bp); if (k < 0 || k >= h.K) return null; const a = FS.a + k * FS.prm.fanoBin;
    return `<b>${timeLabel(a, FS.T)} – ${timeLabel(a + FS.prm.fanoBin, FS.T)}</b><br>${fmt(h.r[k] * FS.prm.fanoBin, 0)} spikes · ${fmt(h.r[k], 2)} spikes/s`; });
  hov("fsIsi", x => { const h = hit.fsIsi; if (!h) return null; for (let k = 0; k < h.ih.nb; k++) { const x0 = h.f.fx(h.ih.edge(k)), x1 = h.f.fx(h.ih.edge(k + 1)); if (x >= x0 && x < x1) { const f = v => v < 1 ? fmt(v * 1000, v < 0.01 ? 2 : 1) + " ms" : fmt(v, 2) + " s";
      return `<b>${f(h.ih.edge(k))} – ${f(h.ih.edge(k + 1))}</b><br>${fmt(h.ih.h[k] * 100, 2)} % of ISIs (${fmt(Math.round(h.ih.h[k] * h.ih.n), 0)})`; } } return null; });
  hov("fsAcg", x => { const h = hit.fsAcg; if (!h) return null; const k = Math.floor((x - h.f.L) / h.f.pw * 2 * h.a.K); if (k < 0 || k >= 2 * h.a.K) return null; const l0 = (k - h.a.K) * h.bms;
    return `<b>${fmt(l0, 1)} to ${fmt(l0 + h.bms, 1)} ms</b><br>${fmt(h.a.c[k], 2)} spikes/s`; });
  const csv = (name, header, rows) => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([[header, ...rows].join("\n")], { type: "text/csv" })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); };
  const base = () => F.name.replace(/\.spk$/i, "");
  const num = (v, d) => Number.isFinite(v) ? v.toFixed(d) : "";
  $("fsCsvEl").addEventListener("click", () => {
    if (!F) return; const p = FS.prm;
    csv(`${base()}_electrode_stats.csv`, `electrode,well,treatment,active,spikes,period_start_s,period_s,rate_hz,fano,median_isi_ms,mean_isi_ms,cv,cv2,lv,isi_below_${p.shortIsi * 1000}ms_pct,n_isi`,
      F.rows.map(c => { const s = FS.el[c], ch = F.channels[c], tr = F.wells[ch.well] && F.wells[ch.well].what || "";
        return [ch.label, ch.well, `"${tr.replace(/"/g, '""')}"`, s.active ? 1 : 0, s.N, num(FS.a, 6), num(FS.T, 6), num(s.rate, 6), num(s.fano, 6), num(s.medIsi * 1000, 6), num(s.meanIsi * 1000, 6), num(s.cv, 6), num(s.cv2, 6), num(s.lv, 6), num(s.shortPct, 6), s.nIsi].join(","); }));
  });
  $("fsCsvWell").addEventListener("click", () => {
    if (!F) return; const keys = ["rate", "cv", "cv2", "lv", "medIsi", "fano"];
    csv(`${base()}_well_stats.csv`, "well,treatment,active,recorded,well_rate_spikes_per_s," + keys.flatMap(k => [`${k}_median`, `${k}_mean`, `${k}_sd`, `${k}_n`]).join(","),
      FS.wells.map(w => [w.well, `"${w.treat.replace(/"/g, '""')}"`, w.active, w.recorded, num(w.wellRate, 6), ...keys.flatMap(k => { const o = w[k], f = k === "medIsi" ? 1000 : 1; return [num(o.med * f, 6), num(o.mean * f, 6), num(o.sd * f, 6), o.n]; })].join(",")));
  });
}

/* =====================================================================
 * Network bursts panel
 * ===================================================================== */

const NB = { prm: null, res: new Map(), well: null, v0: 0, v1: 1, key: "", hit: null, ov: null, hist: null };
const NB_NAME = { isin: "ISI<sub>N</sub>", group: "Burst grouping", recruit: "Recruitment" };
const NB_SHORT = { isin: "ISI_N", group: "Grouping", recruit: "Recruit." };
const NB_COL = { isin: "nb1", group: "nb2", recruit: "nb3" };
function nbParams() {
  const num = (id, def, lo) => { const v = parseFloat($(id).value); return Number.isFinite(v) && v >= lo ? v : def; };
  return { minEl: Math.round(num("nbMinEl", 3, 1)), isinN: Math.round(num("nbN", 10, 2)), isinThr: num("nbThr", 25, 0.1) / 1000, isinMinSpikes: Math.round(num("nbMinSp", 10, 1)),
    burstIsi: num("nbBIsi", 100, 1) / 1000, burstMinSpikes: Math.round(num("nbBMin", 5, 2)), onsetGap: num("nbGap", 100, 1) / 1000, groupPct: Math.min(100, num("nbGPct", 20, 0)),
    recBin: num("nbBin", 50, 1) / 1000, recPct: Math.min(100, num("nbRPct", 25, 0)) };
}
const wellOrder = (a, b) => a.charCodeAt(0) - b.charCodeAt(0) || parseInt(a.slice(1)) - parseInt(b.slice(1));
function computeNB() {
  NB.prm = nbParams(); NB.res = new Map();
  for (const w of [...F.recordedWells].sort(wellOrder)) {
    const trains = F.rows.filter(c => F.channels[c].well === w && FS.el[c].active).sort((p, q) => p - q).map(c => ({ c, t: periodTimes(c) }));
    const r = nbWell(trains, FS.a, FS.T, NB.prm);
    r.trains = trains; r.rec = nbRecruitment(trains, FS.a, FS.T, NB.prm.recBin); r.recReq = nbRequired(NB.prm.recPct, trains.length, NB.prm.minEl);
    r.groupReq = nbRequired(NB.prm.groupPct, trains.length, NB.prm.minEl); r.isinD = nbIsinValues(r.s, NB.prm.isinN);
    NB.res.set(w, r);
  }
  const sel = $("nbWell"), wells = [...NB.res.keys()];
  if (!NB.res.has(NB.well)) { NB.well = wells.reduce((b, w) => !b || NB.res.get(w).active > NB.res.get(b).active ? w : b, null); }
  sel.innerHTML = wells.map(w => `<option value="${w}"${w === NB.well ? " selected" : ""}>${w}${F.wells[w] && F.wells[w].what ? " · " + esc(F.wells[w].what) : ""} · ${NB.res.get(w).active} active</option>`).join("");
  const key = `${FS.a}|${FS.T}`; if (key !== NB.key) { NB.key = key; NB.v0 = FS.a; NB.v1 = FS.a + Math.min(FS.T, 60); }
}
function drawNB() {
  if (!NB.prm) return;
  const r = NB.res.get(NB.well), p = NB.prm;
  const per = $("fsPeriod").value === "window" ? `raster window ${timeLabel(FS.a, FS.T)} – ${timeLabel(FS.a + FS.T, FS.T)}` : `whole recording (${fmtTime(FS.T)})`;
  $("nbSub").innerHTML = !r ? "" : `Well ${NB.well}: ${r.active} active electrodes (≥ ${fmt(FS.prm.activePerMin, 2)} spikes/min), ${fmt(r.spikes, 0)} spikes in the ${per}. Period and active threshold are those of Firing statistics. `
    + `Required electrodes: grouping ≥ ${r.groupReq}, recruitment ≥ ${r.recReq}, ISI<sub>N</sub> ≥ ${p.minEl}.`;
  drawNbTimeline(); drawNbOverview(); drawNbHist(); drawNbAgree(); drawNbTable();
}
function nbSetView(a, b) {
  const lo = FS.a, hi = FS.a + FS.T, span = Math.min(Math.max(b - a, 2 / F.fs), hi - lo);
  NB.v0 = Math.max(lo, Math.min(a, hi - span)); NB.v1 = NB.v0 + span; drawNbTimeline(); drawNbOverview();
}
function drawNbTimeline() {
  if (folded("nbCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const r = NB.res.get(NB.well), cv = $("nbTl"); NB.hit = null; if (!r) return;
  const A = r.active, rowH = A ? Math.max(3, Math.min(8, Math.floor(260 / A))) : 0, TOP = 64, G = 10, LANE = 18;
  const yR = TOP + G, yL = yR + A * rowH + G, H = yL + 3 * LANE + 26;
  cv.style.height = H + "px"; const { ctx, W } = setup(cv);
  const L = RL, pw = W - RL - RR, v0 = NB.v0, v1 = NB.v1, span = v1 - v0, fx = t => L + (t - v0) / span * pw;
  for (const [id, v] of [["nbFrom", v0], ["nbTo", v1]]) { const el = $(id); if (document.activeElement !== el) { el.value = +v.toFixed(4); el.removeAttribute("aria-invalid"); } }
  // grid and time axis
  const ticks = niceTicks(v0, v1, 8); ctx.font = `13px ${C.font}`; ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const v of ticks) { const x = Math.round(fx(v)) + .5; ctx.strokeStyle = C.grid; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, yL + 3 * LANE); ctx.stroke(); ctx.fillStyle = C.muted; ctx.fillText(timeLabel(v, span), x, yL + 3 * LANE + 7); }
  // 1. electrodes firing per recruitment bin (one bar per bin, or the largest bin under each pixel column)
  const b = NB.prm.recBin, K = r.rec.length, k0 = Math.max(0, Math.floor((v0 - FS.a) / b)), k1 = Math.min(K - 1, Math.floor((v1 - FS.a) / b));
  const fy = n => TOP - (A ? n / A : 0) * (TOP - 4);
  ctx.fillStyle = C.accent; const bpx = pw * b / span;
  if (bpx >= 2) { for (let k = k0; k <= k1; k++) if (r.rec[k]) { const x0 = Math.max(L, fx(FS.a + k * b)), x1 = Math.min(L + pw, fx(FS.a + (k + 1) * b)); if (x1 > x0) ctx.fillRect(x0 + (bpx > 4 ? .5 : 0), fy(r.rec[k]), Math.max(1, x1 - x0 - (bpx > 4 ? 1 : 0)), TOP - fy(r.rec[k])); } }
  else { const col = new Int32Array(Math.ceil(pw) + 1); for (let k = k0; k <= k1; k++) { const x = Math.floor(fx(FS.a + (k + .5) * b) - L); if (x >= 0 && x < col.length && r.rec[k] > col[x]) col[x] = r.rec[k]; }
    for (let x = 0; x < col.length; x++) if (col[x]) ctx.fillRect(L + x, fy(col[x]), 1, TOP - fy(col[x])); }
  ctx.strokeStyle = C.axis; ctx.beginPath(); ctx.moveTo(L, TOP + .5); ctx.lineTo(L + pw, TOP + .5); ctx.stroke();
  ctx.strokeStyle = C.nb3; ctx.lineWidth = 2; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(L, Math.round(fy(r.recReq)) + .5); ctx.lineTo(L + pw, Math.round(fy(r.recReq)) + .5); ctx.stroke(); ctx.setLineDash([]); ctx.lineWidth = 1;
  ctx.fillStyle = C.muted; ctx.textAlign = "right"; ctx.textBaseline = "middle";
  ctx.fillText(String(A), L - 6, fy(A) + 4); ctx.fillText("0", L - 6, TOP - 4); ctx.fillStyle = C.nb3; ctx.fillText(String(r.recReq), L - 6, fy(r.recReq));
  // 2. raster of the well's active electrodes
  ctx.fillStyle = C["ink-2"];
  r.trains.forEach((x, row) => { const y = yR + row * rowH; for (let i = lowerBound(x.t, v0); i < x.t.length && x.t[i] <= v1; i++) ctx.fillRect(Math.floor(fx(x.t[i])), y, 1, Math.max(1, rowH - 1)); });
  ctx.fillStyle = C.muted; ctx.textAlign = "right"; ctx.fillText(`${A} el.`, L - 6, yR + A * rowH / 2);
  // 3. one lane of network bursts per method
  const boxes = [];
  NB_METHODS.forEach((m, li) => {
    const y = yL + li * LANE; ctx.fillStyle = C.muted; ctx.textAlign = "right"; ctx.fillText(NB_SHORT[m], L - 6, y + LANE / 2);
    ctx.fillStyle = C.empty; ctx.fillRect(L, y + 3, pw, LANE - 6);
    ctx.fillStyle = C[NB_COL[m]];
    r.bursts[m].forEach((nb, i) => { if (nb.end < v0 || nb.start > v1) return; const x0 = Math.max(L, fx(nb.start)), x1 = Math.min(L + pw, Math.max(fx(nb.end), fx(nb.start) + 2));
      ctx.fillRect(x0, y + 3, Math.max(2, x1 - x0), LANE - 6); boxes.push({ m, i, x0: x0 - 2, x1: Math.max(x1, x0 + 2) + 2, y0: y, y1: y + LANE }); });
  });
  NB.hit = { L, pw, fx, TOP, yR, yL, rowH, LANE, A, boxes, v0, span };
  const vis = NB_METHODS.map(m => r.bursts[m].filter(nb => nb.end >= v0 && nb.start <= v1).length);
  $("nbInfo").textContent = `${timeLabel(v0, span)} – ${timeLabel(v1, span)} · in view: ${vis.map((n, i) => `${NB_SHORT[NB_METHODS[i]]} ${n}`).join(", ")}`;
}
function drawNbOverview() {
  if (folded("nbCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const r = NB.res.get(NB.well), cv = $("nbOv"), { ctx, W, H } = setup(cv); NB.ov = null; if (!r) return;
  const L = RL, pw = W - RL - RR, col = new Float64Array(Math.max(1, Math.floor(pw))), dt = FS.T / col.length;
  for (const v of r.s) { const x = Math.floor((v - FS.a) / dt); if (x >= 0 && x < col.length) col[x]++; }
  let mx = 0; for (const v of col) if (v > mx) mx = v; mx = mx || 1;
  ctx.fillStyle = C.accent; for (let x = 0; x < col.length; x++) if (col[x]) { const h = Math.max(1, Math.sqrt(col[x] / mx) * (H - 4)); ctx.fillRect(L + x, H - h, 1, h); }
  ctx.fillStyle = C.muted; ctx.font = `13px ${C.font}`; ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText("overview", L - 6, H / 2);
  const x0 = L + (NB.v0 - FS.a) / FS.T * pw, x1 = L + (NB.v1 - FS.a) / FS.T * pw;
  ctx.fillStyle = C["accent-wash"]; ctx.fillRect(x0, 0, Math.max(2, x1 - x0), H); ctx.strokeStyle = C.accent; ctx.lineWidth = 1.5; ctx.strokeRect(x0 + .5, .5, Math.max(2, x1 - x0) - 1, H - 1);
  NB.ov = { L, pw, dt, col };
}
function drawNbHist() {
  if (folded("nbCard")) return; // collapsed panel: its canvases have no size; drawn when it is opened again
  const r = NB.res.get(NB.well), cv = $("nbHist"), { ctx, W, H } = setup(cv); NB.hist = null; if (!r) return;
  const ih = nbIsinHist(r.isinD, 10, 1e-4, 100);
  if (!ih.n) { ctx.fillStyle = C.muted; ctx.textAlign = "center"; ctx.fillText(`Fewer than N = ${NB.prm.isinN} spikes`, W / 2, H / 2); return; }
  const mx = Math.max(...ih.h) / ih.n * 100 * 1.08 || 1, yt = niceTicks(0, mx, 4);
  const f = axesFrame(ctx, W, H, { x0: 1e-4, x1: 100, y0: 0, y1: Math.max(mx, yt[yt.length - 1]), logX: true, xTicks: [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], yTicks: yt,
    xFmt: v => v < 1e-3 ? "0.1 ms" : v < 1 ? fmt(v * 1000, 0) + " ms" : fmt(v, 0) + " s", yFmt: v => fmt(v, 1), xLabel: `ISI_N = s(i+${NB.prm.isinN - 1}) − s(i), log scale`, yLabel: "% of windows" });
  ctx.fillStyle = C.accent;
  for (let k = 0; k < ih.nb; k++) if (ih.h[k]) { const x0 = f.fx(ih.edge(k)), x1 = f.fx(ih.edge(k + 1)), y = f.fy(ih.h[k] / ih.n * 100); ctx.fillRect(x0 + .5, y, Math.max(1, x1 - x0 - 1), f.T + f.ph - y); }
  const xs = f.fx(NB.prm.isinThr); if (NB.prm.isinThr > 1e-4 && NB.prm.isinThr < 100) { ctx.strokeStyle = C.nb1; ctx.lineWidth = 2; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(xs, f.T); ctx.lineTo(xs, f.T + f.ph); ctx.stroke(); ctx.setLineDash([]); }
  NB.hist = { f, ih };
}
function drawNbAgree() {
  const r = NB.res.get(NB.well); if (!r) { $("nbAgree").innerHTML = ""; return; }
  $("nbAgree").innerHTML = `<thead><tr><th>Network bursts of ↓ overlapping →</th><th class="n">n</th>${NB_METHODS.map(m => `<th class="n">${NB_NAME[m]}</th>`).join("")}</tr></thead><tbody>`
    + NB_METHODS.map(x => `<tr><td><span class="sw" style="background:var(--${NB_COL[x]});display:inline-block;vertical-align:-1px;margin-right:6px"></span>${NB_NAME[x]}</td><td class="n">${r.summary[x].n}</td>`
      + NB_METHODS.map(y => `<td class="n">${x === y ? "—" : Number.isFinite(r.agreement[x][y]) ? fmt(r.agreement[x][y], 1) + " %" : "–"}</td>`).join("") + "</tr>").join("") + "</tbody>";
}
function drawNbTable() {
  const cols = [["n", "Network bursts", 0, 1], ["perMin", "per min", 2, 1], ["durMedian", "Median duration (ms)", 0, 1000], ["durMean", "Mean duration (ms)", 0, 1000], ["ibiMean", "Mean interval (s)", 2, 1],
    ["ibiCv", "Interval CV", 2, 1], ["spikesMedian", "Median spikes", 0, 1], ["electrodesMedian", "Median electrodes", 1, 1], ["pctSpikes", "Spikes in NBs (%)", 1, 1]];
  $("nbTable").tHead.innerHTML = `<tr><th>Treatment</th><th>Well</th><th class="n">Active</th><th>Method</th>${cols.map(c => `<th class="n">${c[1]}</th>`).join("")}</tr>`;
  const wells = [...NB.res.keys()].sort((a, b) => ((F.wells[a] && F.wells[a].what) || "~").localeCompare((F.wells[b] && F.wells[b].what) || "~") || wellOrder(a, b));
  $("nbTable").tBodies[0].innerHTML = wells.map(w => { const r = NB.res.get(w), tr = F.wells[w] && F.wells[w].what || "no treatment";
    return NB_METHODS.map((m, i) => `<tr data-w="${w}" class="${i === 0 ? "grp" : ""}${w === NB.well ? " sel" : ""}${r.active ? "" : " inactive"}" style="cursor:pointer">`
      + `<td>${i === 0 ? esc(tr) : ""}</td><td>${i === 0 ? w : ""}</td><td class="n">${i === 0 ? r.active : ""}</td><td><span class="sw" style="background:var(--${NB_COL[m]});display:inline-block;vertical-align:-1px;margin-right:6px"></span>${NB_NAME[m]}</td>`
      + cols.map(([k, , d, f]) => `<td class="n">${fmt(r.summary[m][k] * f, d)}</td>`).join("") + "</tr>").join(""); }).join("");
}
function nbTip(m, i) {
  const r = NB.res.get(NB.well), b = r.bursts[m][i];
  return `<b>${NB_NAME[m]} · network burst ${i + 1} of ${r.bursts[m].length}</b><br>${fmt(b.start, 3)} – ${fmt(b.end, 3)} s · ${fmt(b.duration * 1000, 1)} ms<br>${fmt(b.spikes, 0)} spikes on ${b.electrodes} of ${r.active} electrodes`;
}
{
  for (const id of ["nbMinEl", "nbN", "nbThr", "nbMinSp", "nbBIsi", "nbBMin", "nbGap", "nbGPct", "nbBin", "nbRPct"]) $(id).addEventListener("change", () => { if (F) { computeNB(); drawNB(); } });
  $("nbWell").addEventListener("change", e => { NB.well = e.target.value; drawNB(); });
  $("nbValley").addEventListener("click", () => {
    const r = NB.res.get(NB.well); if (!r) return; const ih = nbIsinHist(r.isinD, 10, 1e-4, 100), v = nbValley(ih.h, 1e-4, 10);
    if (!Number.isFinite(v)) { $("nbValley").title = "No valley found in the histogram"; return; }
    $("nbThr").value = +(v * 1000).toPrecision(6); computeNB(); drawNB();
  });
  const apply = () => { if (!F) return; const a = parseTime($("nbFrom").value), b = parseTime($("nbTo").value);
    $("nbFrom").setAttribute("aria-invalid", String(!Number.isFinite(a))); $("nbTo").setAttribute("aria-invalid", String(!Number.isFinite(b)));
    if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) return; $("nbFrom").blur(); $("nbTo").blur(); nbSetView(Math.min(a, b), Math.max(a, b)); };
  for (const id of ["nbFrom", "nbTo"]) { $(id).addEventListener("keydown", e => { if (e.key === "Enter") { e.preventDefault(); apply(); } }); $(id).addEventListener("change", apply); }
  $("nbPrev").addEventListener("click", () => F && nbSetView(NB.v0 - (NB.v1 - NB.v0), NB.v0));
  $("nbNext").addEventListener("click", () => F && nbSetView(NB.v1, NB.v1 + (NB.v1 - NB.v0)));
  $("nbWhole").addEventListener("click", () => F && nbSetView(FS.a, FS.a + FS.T));
  const tl = $("nbTl"), pos = e => { const rc = tl.getBoundingClientRect(); return [e.clientX - rc.left, e.clientY - rc.top]; };
  const boxAt = (x, y) => NB.hit && NB.hit.boxes.find(b => x >= b.x0 && x <= b.x1 && y >= b.y0 && y < b.y1);
  tl.addEventListener("mousemove", e => {
    const h = NB.hit; if (!h) return hideTip(); const [x, y] = pos(e), bx = boxAt(x, y); if (bx) return showTip(e, nbTip(bx.m, bx.i));
    if (x < h.L || x > h.L + h.pw) return hideTip(); const t = h.v0 + (x - h.L) / h.pw * h.span, r = NB.res.get(NB.well);
    if (y <= h.TOP) { const k = Math.floor((t - FS.a) / NB.prm.recBin); if (k < 0 || k >= r.rec.length) return hideTip(); const a = FS.a + k * NB.prm.recBin;
      return showTip(e, `<b>${fmt(a, 3)} – ${fmt(a + NB.prm.recBin, 3)} s</b><br>${r.rec[k]} of ${h.A} electrodes fire (threshold ${r.recReq})`); }
    if (y >= h.yR && y < h.yR + h.A * h.rowH) { const row = Math.floor((y - h.yR) / h.rowH), c = r.trains[row].c; return showTip(e, `<b>${F.channels[c].label}</b> · ${fmt(t, 3)} s<br>${fmt(r.trains[row].t.length, 0)} spikes in the period`); }
    hideTip();
  });
  tl.addEventListener("mouseleave", hideTip);
  tl.addEventListener("click", e => { const [x, y] = pos(e), bx = boxAt(x, y); if (!bx) return; const b = NB.res.get(NB.well).bursts[bx.m][bx.i], w = Math.max(1, 3 * b.duration), c = (b.start + b.end) / 2; nbSetView(c - w / 2, c + w / 2); });
  const ov = $("nbOv"); let drag = null;
  const tOf = x => FS.a + Math.min(Math.max((x - NB.ov.L) / NB.ov.pw, 0), 1) * FS.T, ox = e => e.clientX - ov.getBoundingClientRect().left;
  ov.addEventListener("pointerdown", e => { if (!F || !NB.ov || e.button !== 0) return; drag = { x: ox(e) }; ov.setPointerCapture(e.pointerId); e.preventDefault(); });
  ov.addEventListener("pointermove", e => { if (!drag) { if (NB.ov) { const x = ox(e), k = Math.floor((x - NB.ov.L)); if (k >= 0 && k < NB.ov.col.length) showTip(e, `<b>${fmtTime(tOf(x))}</b><br>${fmt(NB.ov.col[k] / NB.ov.dt, 0)} spikes/s in the well · drag to choose the view`); else hideTip(); } return; }
    const x = ox(e); if (Math.abs(x - drag.x) > 3) { const a = tOf(Math.min(drag.x, x)), b = tOf(Math.max(drag.x, x)); NB.v0 = a; NB.v1 = b; drawNbOverview(); } });
  ov.addEventListener("pointerup", e => { if (!drag) return; const x = ox(e), d = drag; drag = null;
    if (Math.abs(x - d.x) <= 3) { const s = NB.v1 - NB.v0, c = tOf(x); nbSetView(c - s / 2, c + s / 2); } else nbSetView(tOf(Math.min(d.x, x)), tOf(Math.max(d.x, x))); });
  ov.addEventListener("mouseleave", hideTip);
  $("nbHist").addEventListener("mousemove", e => { const h = NB.hist; if (!h) return; const x = e.clientX - e.target.getBoundingClientRect().left;
    for (let k = 0; k < h.ih.nb; k++) { const x0 = h.f.fx(h.ih.edge(k)), x1 = h.f.fx(h.ih.edge(k + 1)); if (x >= x0 && x < x1) { const f = v => v < 1 ? fmt(v * 1000, v < 0.01 ? 2 : 1) + " ms" : fmt(v, 2) + " s";
      return showTip(e, `<b>${f(h.ih.edge(k))} – ${f(h.ih.edge(k + 1))}</b><br>${fmt(h.ih.h[k] / h.ih.n * 100, 2)} % of windows (${fmt(h.ih.h[k], 0)})`); } } hideTip(); });
  $("nbHist").addEventListener("mouseleave", hideTip);
  $("nbTable").addEventListener("click", e => { const tr = e.target.closest("tr[data-w]"); if (tr && tr.dataset.w !== NB.well) { NB.well = tr.dataset.w; $("nbWell").value = NB.well; drawNB(); } });
  const csv = (name, rows) => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([rows.join("\n")], { type: "text/csv" })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); };
  const prmLine = () => { const p = NB.prm; return `# period ${FS.a}-${FS.a + FS.T} s; active >= ${FS.prm.activePerMin} spikes/min; electrodes >= ${p.minEl}; ISI_N: N=${p.isinN}, ISI_N <= ${p.isinThr * 1000} ms, spikes >= ${p.isinMinSpikes}; `
    + `grouping: burst ISI <= ${p.burstIsi * 1000} ms, spikes per burst >= ${p.burstMinSpikes}, onset gap <= ${p.onsetGap * 1000} ms, electrodes >= ${p.groupPct} % of active; recruitment: bin ${p.recBin * 1000} ms, electrodes >= ${p.recPct} % of active`; };
  const q = s => `"${String(s).replace(/"/g, '""')}"`, num = (v, d) => Number.isFinite(v) ? v.toFixed(d) : "";
  $("nbCsv").addEventListener("click", () => { if (!F || !NB.prm) return; const rows = [prmLine(), "well,treatment,method,index,start_s,end_s,duration_ms,spikes,electrodes"];
    for (const [w, r] of NB.res) for (const m of NB_METHODS) r.bursts[m].forEach((b, i) => rows.push([w, q(F.wells[w] && F.wells[w].what || ""), m, i + 1, b.start.toFixed(6), b.end.toFixed(6), (b.duration * 1000).toFixed(3), b.spikes, b.electrodes].join(",")));
    csv(`${F.name.replace(/\.spk$/i, "")}_network_bursts.csv`, rows); });
  $("nbCsvSum").addEventListener("click", () => { if (!F || !NB.prm) return; const rows = [prmLine(), "well,treatment,method,active,n,per_min,dur_median_ms,dur_mean_ms,ibi_mean_s,ibi_cv,spikes_median,electrodes_median,pct_spikes"];
    for (const [w, r] of NB.res) for (const m of NB_METHODS) { const s = r.summary[m]; rows.push([w, q(F.wells[w] && F.wells[w].what || ""), m, r.active, s.n, num(s.perMin, 6), num(s.durMedian * 1000, 6), num(s.durMean * 1000, 6), num(s.ibiMean, 6), num(s.ibiCv, 6), num(s.spikesMedian, 6), num(s.electrodesMedian, 6), num(s.pctSpikes, 6)].join(",")); }
    csv(`${F.name.replace(/\.spk$/i, "")}_network_burst_summary.csv`, rows); });
}

/* =====================================================================
 * Recording details
 * ===================================================================== */
function drawInfo() {
  const h = F.header, d = x => x ? x.toLocaleString() : "–";
  const blocks = [["File", [["Name", F.name], ["Size", fmt(F.size / 1048576, 1) + " MB"], ["Header version", F.version],
    ["Data set", `${h.Name} / ${h.DataSetName || ""} (${DATA_TYPES[h.DataType]})`], ["Experiment start", d(h.ExperimentStartTime)],
    ["Recording start", d(h.BlockVectorStartTime)], ["Duration", fmt(F.duration, 3) + " s"], ["Sampling frequency", fmt(h.SamplingFrequency, 1) + " Hz"],
    ["Voltage scale", h.VoltageScale.toExponential(6) + " V per sample unit"], ["Samples per spike", `${h.NumSamplesPerBlock} (${fmt(F.preTrigMs, 2)} ms before the threshold crossing)`]]],
    ["Metadata", Object.entries(F.meta)],
    ["Wells", F.wellNames.filter(w => F.wells[w]).map(w => { const t = F.wells[w]; return [w, `${t.what || "–"}${F.recordedWells.has(w) ? "" : " (not recorded)"}`]; })]];
  for (const s of F.settings) blocks.push([s.title, s.rows]);
  if (F.warnings.length) blocks.push(["Warnings", F.warnings.map(w => ["", w])]);
  $("info").innerHTML = blocks.filter(b => b[1].length).map(([t, rows]) => `<div><h4>${esc(t)}</h4><table>${rows.map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join("")}</table></div>`).join("");
}

/* =====================================================================
 * Redraw, interaction, loading
 * ===================================================================== */
function drawAll() { if (!F) return; readColors(); hit.plateKey = null; computeStats(); drawTiles(); drawTimeline(); renderReplay(); drawSelection(); drawRaster(); drawWaveforms(); drawStats(); drawInfo(); drawPlayheads(); }

$("plate").tabIndex = 0;
$("plate").addEventListener("mousemove", e => { if (!F) return; const h = plateHit(e); h.el ? showTip(e, electrodeTip(h.el.c)) : h.well ? showTip(e, wellTip(h.well.name)) : hideTip(); });
$("plate").addEventListener("mouseleave", hideTip);
$("plate").addEventListener("click", e => { if (!F) return; const h = plateHit(e); if (h.el) select(h.el.c); });
$("plate").addEventListener("keydown", e => {
  if (!F || S.sel < 0) return;
  const k = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[e.key];
  if (k) { e.preventDefault(); moveSel(...k); }
});
$("plateMetric").addEventListener("change", () => F && drawPlate());
$("themeBtn").addEventListener("click", () => {
  const dark = document.documentElement.dataset.theme ? document.documentElement.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  document.documentElement.dataset.theme = dark ? "light" : "dark"; redraw();
});
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => redraw());
new ResizeObserver(([e]) => document.documentElement.style.setProperty("--topbar-h", Math.ceil(e.target.getBoundingClientRect().height) + "px")).observe(document.querySelector(".topbar"));
function redrawCharts() { if (F) { hit.plateKey = null; drawTimeline(); drawPlate(); drawRaster(); drawWaveforms(); drawPlayheads(); if (FS.prm) drawStatCharts(); if (NB.prm) { drawNbTimeline(); drawNbOverview(); drawNbHist(); } } }
let resizeT; new ResizeObserver(() => { clearTimeout(resizeT); resizeT = setTimeout(redrawCharts, 80); }).observe($("mapCard"));

/* ---------- collapsible panels (plate map and below); the choice applies to every tab and is remembered in this browser ---------- */
const COLL_CARDS = ["mapCard", "wfCard", "fsCard", "nbCard"], COLL_KEY = "meaPipeline.collapsed";
const collSaved = (() => { try { return new Set(JSON.parse(localStorage.getItem(COLL_KEY) || "[]").filter(id => COLL_CARDS.includes(id))); } catch (e) { return new Set(); } })();
function setCollapsed(id, on) {
  const card = $(id), b = card.querySelector(".coll");
  card.classList.toggle("collapsed", on); b.setAttribute("aria-expanded", String(!on)); b.title = on ? "Show this panel" : "Hide this panel";
  on ? collSaved.add(id) : collSaved.delete(id);
  try { localStorage.setItem(COLL_KEY, JSON.stringify([...collSaved])); } catch (e) {}
}
for (const id of COLL_CARDS) {
  const card = $(id), h = card.querySelector(".toolbar h3"), b = document.createElement("button");
  b.className = "coll"; b.type = "button"; b.setAttribute("aria-label", "Show or hide " + h.textContent);
  b.innerHTML = '<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true"><path d="M3.5 6l4.5 4.5L12.5 6" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>'; h.before(b);
  const toggle = () => { const on = !card.classList.contains("collapsed"); setCollapsed(id, on); if (!on) requestAnimationFrame(redrawCharts); };
  b.addEventListener("click", toggle); h.classList.add("colltitle"); h.addEventListener("click", e => { if (!e.target.closest("button")) toggle(); });
  setCollapsed(id, collSaved.has(id));
}
function setAllCollapsed(on) { for (const id of COLL_CARDS) setCollapsed(id, on); collButtons(); if (!on) requestAnimationFrame(redrawCharts); }
function collButtons() { const n = COLL_CARDS.filter(id => $(id).classList.contains("collapsed")).length; $("collAll").disabled = n === COLL_CARDS.length; $("expAll").disabled = n === 0; }
$("collAll").addEventListener("click", () => setAllCollapsed(true));
$("expAll").addEventListener("click", () => setAllCollapsed(false));
for (const id of COLL_CARDS) $(id).querySelector(".coll").addEventListener("click", collButtons);
for (const id of COLL_CARDS) $(id).querySelector(".colltitle").addEventListener("click", collButtons);
collButtons();

/* =====================================================================
 * Tabs: one tab per recording (labelled with its file name), and the Multi-recording analysis tab (id "cmp").
 * The analysis panels are a single set of controls and canvases. Switching
 * tab stores the state and control values of the recording on show and
 * restores those of the other one, so each tab keeps its own parameters.
 * ===================================================================== */
const TABS = [];                 // { id, key, title, F, st, ms }
let CUR = null, tabSeq = 0;      // active tab: a tab id, "cmp", or null (nothing open)
const CTRL_IDS = [...document.querySelectorAll("#app input[id], #app select[id]")].map(x => x.id);
const DYN_SEL = ["rWell", "wfWell", "nbWell"];   // their options depend on the recording
const ctrlSnap = () => Object.fromEntries(CTRL_IDS.map(id => { const x = $(id); return [id, x.type === "checkbox" ? x.checked : x.value]; }));
const CTRL_DEFAULT = ctrlSnap();
// controls that change the network bursts, with their labels (the Multi-recording analysis tab lists the ones that differ from the defaults)
const NB_CTRL = [["fsPeriod", "period"], ["fsActive", "active ≥ (spikes/min)"], ["nbMinEl", "electrodes ≥ m"], ["nbN", "N"], ["nbThr", "ISI_N ≤ (ms)"], ["nbMinSp", "ISI_N spikes ≥"],
  ["nbBIsi", "burst ISI ≤ (ms)"], ["nbBMin", "spikes per burst ≥"], ["nbGap", "onset gap τ ≤ (ms)"], ["nbGPct", "grouping electrodes ≥ (%)"], ["nbBin", "recruitment bin (ms)"], ["nbRPct", "recruitment electrodes ≥ (%)"]];
const fileKey = f => `${f.name}|${f.size}|${f.lastModified}`;

function freshState(R) {
  let best = -1; for (const c of R.recorded) if (best < 0 || R.counts[c] > R.counts[best]) best = c;
  return { S: { sel: best, t0: 0, t1: R.duration }, P: { t: 0, active: false, speed: +CTRL_DEFAULT.speedSel, win: +CTRL_DEFAULT.winSel },
    FS: { prm: null, a: 0, T: R.duration, sort: { key: "well", dir: 1 } }, NB: { prm: null, res: new Map(), well: null, v0: 0, v1: 1, key: "" }, ctrl: { ...CTRL_DEFAULT } };
}
function saveTab() { // the globals and controls -> the active tab
  const t = TABS.find(t => t.id === CUR); if (!t || !F) return;
  if (P.playing) pause();
  t.st = { S: { ...S }, P: { t: P.t, active: P.active, speed: P.speed, win: P.win }, FS: { prm: FS.prm, a: FS.a, T: FS.T, sort: { ...FS.sort }, el: FS.el, wells: FS.wells },
    NB: { prm: NB.prm, res: NB.res, well: NB.well, v0: NB.v0, v1: NB.v1, key: NB.key }, ctrl: ctrlSnap() };
}
function installTab(t) { // a tab -> the globals and controls (no drawing)
  const st = t.st; F = t.F;
  for (const id of CTRL_IDS) if (!DYN_SEL.includes(id) && id in st.ctrl) { const x = $(id); if (x.type === "checkbox") x.checked = st.ctrl[id]; else x.value = st.ctrl[id]; }
  Object.assign(S, st.S);
  fillWellSelect(); fillWaveWellSelect();
  for (const id of ["rWell", "wfWell"]) if ([...$(id).options].some(o => o.value === st.ctrl[id])) $(id).value = st.ctrl[id];
  cancelAnimationFrame(P.raf); Object.assign(P, st.P, { playing: false, counts: null }); replayScale();
  Object.assign(FS, { prm: st.FS.prm, a: st.FS.a, T: st.FS.T, sort: { ...st.FS.sort }, el: st.FS.el || [], wells: st.FS.wells || [] });
  Object.assign(NB, { prm: st.NB.prm, res: st.NB.res, well: st.NB.well, v0: st.NB.v0, v1: st.NB.v1, key: st.NB.key });
  WF.cache = null; hit.plateKey = null;
}
function activate(id) {
  if (id !== CUR) saveTab();
  hideTip(); CUR = id;
  const t = TABS.find(t => t.id === id);
  $("drop").hidden = !!t || id === "cmp"; $("app").hidden = !t; $("cmp").hidden = id !== "cmp";
  if (t) { installTab(t); $("fname").textContent = t.F.name; drawAll(); }
  else {
    if (P.playing) pause(); cancelAnimationFrame(P.raf); F = null;
    $("fname").textContent = id === "cmp" ? `Multi-recording analysis · ${CMP.sel.size} of ${TABS.length} recordings selected` : "No file loaded";
    if (id === "cmp") drawCmp();
  }
  renderTabs();
}
function renderTabs() {
  const tab = (id, label, title, close) => `<div class="tab${id === CUR ? " on" : ""}"><button role="tab" id="tab-${id}" aria-selected="${id === CUR}" data-tab="${id}" title="${esc(title)}">${label}</button>`
    + (close ? `<button class="x" data-close="${id}" aria-label="Close ${esc(title)}" title="Close this tab">×</button>` : "") + `</div>`;
  $("tabs").innerHTML = TABS.map(t => tab(t.id, esc(t.title), t.F.name, true)).join("")
    + tab("cmp", `Multi-recording analysis${CMP.sel.size ? ` (${CMP.sel.size})` : ""}`, "Compare the analyses of several recordings", false);
}
function closeTab(id) {
  const i = TABS.findIndex(t => t.id === id); if (i < 0) return;
  if (CUR === id) { if (P.playing) pause(); cancelAnimationFrame(P.raf); F = null; CUR = null; }
  TABS.splice(i, 1); cmpDeselect(id);
  if (CUR === null) activate(TABS.length ? TABS[Math.min(i, TABS.length - 1)].id : null);
  else { renderTabs(); if (CUR === "cmp") activate("cmp"); }
}
const redraw = () => CUR === "cmp" ? drawCmp() : drawAll();

// Read the files one after the other. A new file gets its own tab, analysed with the default parameters;
// a file already open (same name, size and date) is not read again. show: go to the last tab, else stay.
async function openFiles(files, show) {
  files = [...files].filter(f => f); if (!files.length) return [];
  $("err").textContent = "";
  const back = CUR, out = [], errs = [];
  for (const file of files) {
    const known = TABS.find(t => t.key === fileKey(file)); if (known) { out.push(known); continue; }
    $("fname").textContent = `Reading ${file.name}…`; await new Promise(r => setTimeout(r, 0));
    try {
      const buf = await file.arrayBuffer(), t0 = performance.now();
      const R = parseSpk(buf); R.name = file.name; R.size = file.size; prepare(R);
      const n = TABS.filter(t => t.F.name === R.name).length;
      const t = { id: "t" + ++tabSeq, key: fileKey(file), title: n ? `${R.name} (${n + 1})` : R.name, F: R, st: freshState(R) };
      saveTab(); CUR = t.id; TABS.push(t);
      installTab(t); computeStats(); computeNB(); saveTab();
      t.ms = Math.round(performance.now() - t0); out.push(t);
    } catch (e) { console.error(e); errs.push(`Could not read ${file.name}: ${e.message}`); }
  }
  for (const t of out) cmpSelect(t.id, true);
  if (errs.length) $("err").textContent = errs.join("\n") + (TABS.length ? "\nThe recordings already open are unchanged." : "");
  const target = show && out.length ? out[out.length - 1].id : (back === null && out.length ? out[out.length - 1].id : back);
  activate(target);
  return out;
}

/* ---------- Multi-recording analysis (tab id "cmp") ---------- */
const CMP_MAX = 8;
const CMP_AXIS = { isin: "ISI_N", group: "Burst grouping", recruit: "Recruitment" };                               // categorical colours available
const CMP = { sel: new Map(), hit: [], confirm: false, busy: null };          // selected tab id -> colour index (0..7), kept while selected
function cmpSelect(id, quiet) {
  if (CMP.sel.has(id)) return true;
  const used = new Set(CMP.sel.values()); let k = 0; while (used.has(k)) k++;
  if (k >= CMP_MAX) { if (!quiet) $("cmpNote").textContent = `At most ${CMP_MAX} recordings can be compared at once: unselect one first.`; return false; }
  CMP.sel.set(id, k); return true;
}
function cmpDeselect(id) { CMP.sel.delete(id); }
// What the comparison measures, per method. Each is computed from the tab's network bursts (period T of that tab).
const CMP_MEAS = [
  { k: "n", id: "cmpMeasN", name: "Count", axis: "count", unit: "", d: 0, tip: v => `${fmt(v, 0)} network bursts` },
  { k: "rate", id: "cmpMeasR", name: "Rate (per min)", axis: "per min", unit: "/min", d: 2, tip: v => `${fmt(v, 2)} network bursts per minute` },
  { k: "dur", id: "cmpMeasD", name: "Median duration (ms)", axis: "median ms", unit: "ms", d: 1, tip: v => `median duration ${fmt(v, 1)} ms` }];
function cmpValues(t, well) { return nbCompare(t.st.NB.res, t.st.FS.T, well); }
function cmpChanged(t) { // parameters of this tab that differ from the defaults
  const c = t.st.ctrl; return NB_CTRL.filter(([id]) => String(c[id]) !== String(CTRL_DEFAULT[id])).map(([id, l]) => id === "fsPeriod" ? `period: ${timeLabel(t.st.FS.a, t.st.FS.T)} – ${timeLabel(t.st.FS.a + t.st.FS.T, t.st.FS.T)}` : `${l} ${c[id]}`);
}
// analysis parameters reset by "Re-analyse with the default parameters": Firing statistics and Network bursts settings,
// not the view controls (well shown, time view, table options)
const ANALYSIS_IDS = CTRL_IDS.filter(id => /^(fs|nb)/.test(id) && !["nbWell", "nbFrom", "nbTo", "fsOnlyActive", "fsSummary"].includes(id));
const nbSignature = t => JSON.stringify(NB_CTRL.map(([id]) => String(t.st.ctrl[id])));
function drawCmpWarn(sel) { // warning when the selected recordings were not all analysed with the default parameters
  const box = $("cmpWarn"), odd = sel.filter(t => cmpChanged(t).length);
  if (!odd.length || CMP.busy) { box.hidden = !CMP.busy; if (CMP.busy) box.innerHTML = `<b>Re-analysing…</b> ${esc(CMP.busy)}`; return; }
  const differ = new Set(sel.map(nbSignature)).size > 1, names = odd.map(t => `<b>${esc(t.title)}</b>`).join(", ");
  box.hidden = false;
  box.innerHTML = `<span class="wicon" aria-hidden="true">⚠</span><div><p><b>${differ ? "The selected recordings were not analysed with the same parameters" : "The selected recordings were not analysed with the default parameters"}.</b> `
    + `${differ ? "Network burst counts, rates and durations depend on these parameters, so differences between the recordings may come from the settings rather than from the cultures." : "They share the same settings, so they can be compared with each other, but not with results obtained with the defaults."} `
    + `Changed in: ${names} (details under each recording below).</p>`
    + (CMP.confirm
      ? `<p>This sets the Firing statistics and Network bursts parameters of all ${TABS.length} open recording${TABS.length > 1 ? "s" : ""} back to the defaults and recomputes them. The parameters changed in their tabs are lost; the view and the selected electrode are kept.</p>`
        + `<div class="toolbar"><button id="cmpRedoYes" class="primary">Re-analyse all with the defaults</button><button id="cmpRedoNo">Cancel</button></div>`
      : `<div class="toolbar"><button id="cmpRedo">Re-analyse all recordings with the default parameters…</button><span style="color:var(--muted)">or keep the current parameters</span></div>`)
    + `</div>`;
}
async function reanalyseDefaults() {
  const back = CUR; CMP.confirm = false;
  for (const [k, t] of TABS.entries()) {
    CMP.busy = `${k + 1} of ${TABS.length}: ${t.title}`; drawCmpWarn([]); await new Promise(r => setTimeout(r, 0));
    for (const id of ANALYSIS_IDS) t.st.ctrl[id] = CTRL_DEFAULT[id];
    CUR = t.id; installTab(t); computeStats(); computeNB(); saveTab();
  }
  CMP.busy = null; CUR = back; F = null; activate(back);
}
function niceCeil(v) { if (!(v > 0)) return 1; const e = 10 ** Math.floor(Math.log10(v)); for (const m of [1, 2, 2.5, 5, 10]) if (m * e >= v * (1 - 1e-12)) return m * e; return 10 * e; }
function drawCmp() {
  readColors(); $("cmpNote").textContent = "";
  for (const id of [...CMP.sel.keys()]) if (!TABS.some(t => t.id === id)) CMP.sel.delete(id);
  const sel = TABS.filter(t => CMP.sel.has(t.id)).sort((a, b) => CMP.sel.get(a.id) - CMP.sel.get(b.id));
  // recordings
  $("cmpSub").innerHTML = TABS.length ? `${sel.length} of ${TABS.length} open recording${TABS.length > 1 ? "s" : ""} selected (at most ${CMP_MAX}). Tick a recording to add it to the plot; each uses the parameters of its own tab.`
    : `No recording is open. Add .spk files: each one is analysed with the default parameters and gets its own tab.`;
  $("cmpList").innerHTML = TABS.map(t => {
    const on = CMP.sel.has(t.id), ch = cmpChanged(t), nw = t.F.recordedWells.size;
    return `<div class="cmprow"><label><input type="checkbox" data-cmp="${t.id}"${on ? " checked" : ""}> <span class="sw" style="background:${on ? `var(--cat-${CMP.sel.get(t.id) + 1})` : "transparent"}"></span> <span class="nm">${esc(t.title)}</span></label>`
      + `<span class="meta">${fmtTime(t.F.duration)} · ${nw} well${nw === 1 ? "" : "s"} · ${ch.length ? "changed: " + esc(ch.join(", ")) : "default parameters"}</span><span class="grow"></span>`
      + `<button data-go="${t.id}" title="Show its tab">Open tab</button></div>`;
  }).join("");
  drawCmpWarn(sel);
  // wells: all wells recorded in the selected recordings
  const ws = $("cmpWell"), prev = ws.value || "*", wells = [...new Set(sel.flatMap(t => [...t.F.recordedWells]))].sort(wellOrder);
  ws.innerHTML = `<option value="*">All wells (sum)</option>` + wells.map(w => `<option value="${w}">${w}</option>`).join("");
  ws.value = prev === "*" || wells.includes(prev) ? prev : "*";
  const well = ws.value, rows = sel.map(t => ({ t, k: CMP.sel.get(t.id), v: cmpValues(t, well) }));
  const meas = CMP_MEAS.filter(q => $(q.id).checked), axes = meas.flatMap(q => NB_METHODS.map(m => ({ q, m })));
  const nbAxes = axes.map(a => ({ l1: CMP_AXIS[a.m], l2: a.q.axis, col: C[{ n: "meas-n", rate: "meas-r", dur: "meas-d" }[a.q.k]], grp: a.q.k,
    val: r => r.v[a.m][a.q.k], title: `${NB_NAME[a.m]} · ${a.q.name.toLowerCase()}`, tip: a.q.tip }));
  CMP.hit = drawSpider("cmpSpider", rows, nbAxes, $("cmpScale").value, "Tick at least one measure");
  CMP.well = well;
  drawCmpFiring(sel, well);
  // table: one column per axis, grouped by measure
  const tb = $("cmpTable"), nc = 3 + axes.length;
  tb.innerHTML = `<thead><tr><th rowspan="2">Recording</th><th rowspan="2" class="n">Period</th><th rowspan="2" class="n">Active electrodes</th>${meas.map(q => `<th colspan="${NB_METHODS.length}" class="grp m-${q.k}">${q.name}</th>`).join("")}</tr>`
    + `<tr>${axes.map(a => `<th class="n">${NB_NAME[a.m]}</th>`).join("")}</tr></thead><tbody>`
    + (rows.length ? rows.map(({ t, k, v }) => {
      const res = t.st.NB.res, act = well === "*" ? [...res.values()].reduce((s, r) => s + r.active, 0) : res.has(well) ? res.get(well).active : NaN;
      return `<tr><td><span class="sw" style="background:var(--cat-${k + 1})"></span> ${esc(t.title)}</td><td class="n">${fmtTime(t.st.FS.T)}</td><td class="n">${fmt(act, 0)}</td>`
        + axes.map(a => `<td class="n">${fmt(v[a.m][a.q.k], a.q.d)}</td>`).join("") + `</tr>`; }).join("")
      : `<tr><td colspan="${nc}" style="color:var(--muted)">No recording selected</td></tr>`) + `</tbody>`;
  $("cmpNote").textContent = !rows.length ? "" : `In the analysis period of each recording, ${well === "*" ? "over all its wells (counts summed, durations of all network bursts pooled)" : "in well " + well}. `
    + `Rate = 60 · count / period. Duration = end − start of each network burst; median, with no value when there are none.`;
}
// Spider plot: one axis per entry of axes ({ l1, l2, col, grp, val(row), title, tip(v) }), one outline per row (recording).
// Each axis runs from 0 to the round number at or above its largest value ("measure": shared by the axes of one grp). Returns the hover targets.
// Firing per recording: network firing rate, mean electrode firing rate and active electrodes (period and threshold of each tab)
const FIRING_AXES = [
  { k: "net", l1: "Network firing rate", l2: "spikes/s", d: 1, tip: v => `${fmt(v, 1)} spikes/s` },
  { k: "mean", l1: "Mean electrode rate", l2: "Hz, active electrodes", d: 2, tip: v => `${fmt(v, 2)} Hz` },
  { k: "active", l1: "Active electrodes", l2: "count", d: 0, tip: v => `${fmt(v, 0)} electrodes` }];
function drawCmpFiring(sel, well) {
  const rows = sel.map(t => ({ t, k: CMP.sel.get(t.id), f: firingCompare(t.F, t.st.FS.el, t.st.FS.T, well) }));
  const axes = FIRING_AXES.map(a => ({ ...a, grp: a.k, val: r => r.f[a.k], title: a.l1.toLowerCase() }));
  CMP.hitF = drawSpider("cmpSpiderF", rows, axes, "axis", "");
  const th = (rows.length ? rows : []).map(({ t, k, f }) => `<tr><td><span class="sw" style="background:var(--cat-${k + 1})"></span> ${esc(t.title)}</td><td class="n">${fmtTime(t.st.FS.T)}</td>`
    + `<td class="n">${fmt(f.net, 1)}</td><td class="n">${fmt(f.mean, 2)}${Number.isFinite(f.sd) ? ` ± ${fmt(f.sd, 2)}` : ""}</td><td class="n">${fmt(f.active, 0)} / ${f.recorded}</td><td class="n">${fmt(t.st.FS.prm ? t.st.FS.prm.activePerMin : NaN, 2)}</td></tr>`).join("");
  $("cmpTableF").innerHTML = `<thead><tr><th>Recording</th><th class="n">Period</th><th class="n">Network firing rate (spikes/s)</th><th class="n">Mean electrode rate (Hz, mean ± SD)</th><th class="n">Active / recorded electrodes</th><th class="n">Active ≥ (spikes/min)</th></tr></thead><tbody>`
    + (th || `<tr><td colspan="6" style="color:var(--muted)">No recording selected</td></tr>`) + `</tbody>`;
  $("cmpNoteF").textContent = !rows.length ? "" : `${well === "*" ? "Over all recorded wells of each recording (spikes and electrodes pooled)" : "In well " + well}, over the Firing statistics period of each recording. `
    + `Network firing rate = all spikes of the recorded electrodes / period. Mean electrode rate: mean of the active electrodes' rates (± sample SD); active: rate ≥ the threshold of that recording's tab.`;
}
function drawSpider(id, rows, axes, scale, noAxes) {
  const cv = $(id), { ctx, W, H } = setup(cv), n = axes.length, hits = [];
  ctx.textAlign = "center"; ctx.fillStyle = C.muted;
  if (!rows.length) { ctx.fillText("Tick one or more recordings above", W / 2, H / 2); return hits; }
  if (!n) { ctx.fillText(noAxes, W / 2, H / 2); return hits; }
  const cx = W / 2, cy = H / 2 + 4, R = Math.max(20, Math.min(W / 2 - 100, H / 2 - 62));
  const ang = k => -Math.PI / 2 + k * 2 * Math.PI / n, pt = (k, f) => [cx + Math.cos(ang(k)) * R * f, cy + Math.sin(ang(k)) * R * f];
  const maxOf = a => Math.max(0, ...rows.map(r => a.val(r)).filter(Number.isFinite));
  // outer ring of an axis: the round number at or above its largest value; "measure": shared by the methods of one measure
  const top = axes.map(a => niceCeil(scale === "measure" ? Math.max(...axes.filter(b => b.grp === a.grp).map(maxOf)) : maxOf(a)));
  ctx.lineWidth = 1; ctx.strokeStyle = C.grid;
  for (const f of [0.25, 0.5, 0.75, 1]) { ctx.beginPath(); for (let k = 0; k < n; k++) { const [x, y] = pt(k, f); k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); } if (n > 2) ctx.closePath(); ctx.stroke(); }
  ctx.strokeStyle = C.axis; for (let k = 0; k < n; k++) { const [x, y] = pt(k, 1); ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(x, y); ctx.stroke(); }
  // axis labels: method, measure, range (text ink, not series colours)
  axes.forEach((a, k) => {
    const [x, y] = pt(k, 1), an = ang(k), dx = Math.cos(an), dy = Math.sin(an);
    ctx.textAlign = Math.abs(dx) < 0.2 ? "center" : dx > 0 ? "left" : "right";
    const mc = a.col || C.ink; // label colour: the measure group (never a recording's colour)
    const lines = [[`600 14px ${C.font}`, mc, a.l1], [`13px ${C.font}`, a.col || C["ink-2"], a.l2], [`13px ${C.font}`, C.muted, `0 – ${fmt(top[k], 2)}`]].filter(l => l[2]);
    const lh = 16, h = lines.length * lh, lx = x + dx * 10, y0 = y + dy * 10 + (dy < -0.3 ? -h : dy > 0.3 ? 0 : -h / 2);
    ctx.textBaseline = "top"; lines.forEach(([f, c, s], i) => { ctx.font = f; ctx.fillStyle = c; ctx.fillText(s, lx, y0 + i * lh); });
  });
  // one outline per recording through its points (2 px, light fill), ringed markers; a missing value is left out
  for (const r of rows) {
    const col = css(`--cat-${r.k + 1}`);
    const P2 = axes.map((a, k) => { const v = a.val(r); return Number.isFinite(v) ? pt(k, v / top[k]) : null; }), ok = P2.filter(Boolean);
    if (ok.length > 1) {
      ctx.beginPath(); ok.forEach(([x, y], k) => k ? ctx.lineTo(x, y) : ctx.moveTo(x, y)); if (ok.length > 2) ctx.closePath();
      if (ok.length > 2) { ctx.globalAlpha = 0.10; ctx.fillStyle = col; ctx.fill(); ctx.globalAlpha = 1; }
      ctx.lineWidth = 2; ctx.lineJoin = "round"; ctx.strokeStyle = col; ctx.stroke();
    }
    P2.forEach((p, k) => { if (!p) return; ctx.beginPath(); ctx.arc(p[0], p[1], 4.5, 0, 2 * Math.PI); ctx.fillStyle = col; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = C.surface; ctx.stroke();
      hits.push({ x: p[0], y: p[1], r, a: axes[k], k }); });
  }
  return hits;
}
{
  $("tabs").addEventListener("click", e => {
    const x = e.target.closest("[data-close]"); if (x) return closeTab(x.dataset.close);
    const b = e.target.closest("[data-tab]"); if (b && b.dataset.tab !== CUR) activate(b.dataset.tab);
  });
  $("tabs").addEventListener("keydown", e => { // arrow keys move between tabs
    const b = e.target.closest("[data-tab]"); if (!b || !["ArrowLeft", "ArrowRight"].includes(e.key)) return;
    const ids = [...TABS.map(t => t.id), "cmp"], i = ids.indexOf(b.dataset.tab), j = (i + (e.key === "ArrowRight" ? 1 : -1) + ids.length) % ids.length;
    e.preventDefault(); activate(ids[j]); $("tab-" + ids[j]).focus();
  });
  $("cmpList").addEventListener("change", e => {
    const c = e.target.closest("[data-cmp]"); if (!c) return;
    if (c.checked) { if (!cmpSelect(c.dataset.cmp)) { c.checked = false; return; } } else cmpDeselect(c.dataset.cmp);
    activate("cmp");
  });
  $("cmpList").addEventListener("click", e => { const b = e.target.closest("[data-go]"); if (b) activate(b.dataset.go); });
  $("cmpWarn").addEventListener("click", e => {
    const b = e.target.closest("button"); if (!b || CMP.busy) return;
    if (b.id === "cmpRedo") { CMP.confirm = true; drawCmp(); $("cmpRedoNo").focus(); }
    else if (b.id === "cmpRedoNo") { CMP.confirm = false; drawCmp(); }
    else if (b.id === "cmpRedoYes") reanalyseDefaults();
  });
  $("cmpIn").addEventListener("change", e => { const fs = [...e.target.files]; e.target.value = ""; openFiles(fs, false); });
  for (const id of ["cmpWell", "cmpScale", "cmpMeasN", "cmpMeasR", "cmpMeasD"]) $(id).addEventListener("change", drawCmp);
  for (const [id, key] of [["cmpSpider", "hit"], ["cmpSpiderF", "hitF"]]) {
    const sp = $(id);
    sp.addEventListener("mousemove", e => {
      const rc = sp.getBoundingClientRect(), x = e.clientX - rc.left, y = e.clientY - rc.top, H = CMP[key] || [];
      let best = null, bd = 12; for (const h of H) { const d = Math.hypot(h.x - x, h.y - y); if (d < bd) { bd = d; best = h; } }
      if (!best) return hideTip();
      const same = H.filter(h => h.k === best.k && Math.hypot(h.x - best.x, h.y - best.y) < 1);   // coincident points
      const a = best.a; showTip(e, `<b>${a.title}</b> · ${CMP.well === "*" ? "all wells" : "well " + CMP.well}<br>` + same.map(h => `${esc(h.r.t.title)}: ${a.tip(a.val(h.r))}`).join("<br>"));
    });
    sp.addEventListener("mouseleave", hideTip);
  }
  let rz; new ResizeObserver(() => { clearTimeout(rz); rz = setTimeout(() => { if (CUR === "cmp") drawCmp(); }, 80); }).observe($("cmpSpider"));
}

$("fileIn").addEventListener("change", e => { const fs = [...e.target.files]; e.target.value = ""; openFiles(fs, CUR !== "cmp"); });
["dragover", "dragenter"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); $("drop").classList.add("over"); }));
["dragleave", "drop"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); $("drop").classList.remove("over"); }));
document.addEventListener("drop", e => { const fs = [...e.dataTransfer.files]; if (fs.length) openFiles(fs, CUR !== "cmp"); });
renderTabs();
