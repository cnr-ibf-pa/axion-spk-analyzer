// Detect network bursts with the page's own code (parser + STATS and NB blocks of mea_pipeline.html) in Node,
// for compare_network_bursts.py.   node dump_page_nb.js ../mea_pipeline.html recording.spk out.json [t0 t1]
const fs = require("fs");
const [html, spk, out, t0s, t1s] = process.argv.slice(2);
const page = fs.readFileSync(html, "utf8");
const cut = (a, b) => { const i = page.indexOf(a), j = page.indexOf(b, i); if (i < 0 || j < 0) throw new Error(`block ${a} not found`); return page.slice(i, j); };
const pb = page.indexOf("\nfunction parseSpk("), pe = page.indexOf("\n}\n", pb) + 3;
const lib = new Function(page.slice(page.indexOf("const PLATES = {"), pe) + "\n" + cut("/* STATS-BEGIN", "/* STATS-END */") + "\n" + cut("/* NB-BEGIN", "/* NB-END */")
  + "\nreturn { parseSpk, electrodeStats, nbWell, nbIsinValues, nbIsinHist, nbValley };")();
const b0 = fs.readFileSync(spk), F = lib.parseSpk(b0.buffer.slice(b0.byteOffset, b0.byteOffset + b0.byteLength));
const dur = F.header.Duration || F.times[F.n - 1], a = t0s ? +t0s : 0, T = (t1s ? +t1s : dur) - a;
// page defaults (ms and % as typed in the panel)
const prm = { minEl: 3, isinN: 10, isinThr: 25 / 1000, isinMinSpikes: 10, burstIsi: 100 / 1000, burstMinSpikes: 5, onsetGap: 100 / 1000, groupPct: 20, recBin: 50 / 1000, recPct: 25 };
const per = new Map();
for (let i = 0; i < F.n; i++) { const c = F.ch[i]; if (c < 0) continue; const t = F.times[i]; if (t < a || t > a + T) continue; if (!per.has(c)) per.set(c, []); per.get(c).push(t); }
const wells = new Map();
for (const c of [...per.keys()].sort((p, q) => p - q)) {
  const t = Float64Array.from(per.get(c)).sort();
  if (!lib.electrodeStats(t, a, T, { activePerMin: 5, fanoBin: 1, shortIsi: 0.0015, minIsi: 20 }).active) continue;
  const w = F.channels[c].well; if (!wells.has(w)) wells.set(w, []); wells.get(w).push({ c, t });
}
const res = {};
for (const [w, trains] of wells) {
  const r = lib.nbWell(trains, a, T, prm), ih = lib.nbIsinHist(lib.nbIsinValues(r.s, prm.isinN), 10, 1e-4, 100);
  res[w] = { active: r.active, spikes: r.spikes, bursts: r.bursts, summary: r.summary, agreement: r.agreement, valley: lib.nbValley(ih.h, 1e-4, 10), hist: Array.from(ih.h) };
}
fs.writeFileSync(out, JSON.stringify(res, (k, v) => (typeof v === "number" && !Number.isFinite(v)) ? null : v));
console.log(`page network bursts for ${Object.keys(res).length} wells with active electrodes -> ${out}`);
