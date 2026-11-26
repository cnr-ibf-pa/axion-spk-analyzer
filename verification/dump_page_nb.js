// Detect network bursts with the page's own code (js/spk_reader.js + js/analysis.js) in Node,
// for compare_network_bursts.py.   node dump_page_nb.js .. recording.spk out.json [t0 t1]   (.. = the app folder)
const fs = require("fs");
const [appArg, spk, out, t0s, t1s] = process.argv.slice(2);
// The page's analysis code: js/spk_reader.js + js/analysis.js of the app folder (or of the folder of mea_pipeline.html).
const path = require("path");
const app = /\.html?$/i.test(appArg) ? path.dirname(appArg) : appArg;
const src = ["spk_reader.js", "analysis.js"].map(f => fs.readFileSync(path.join(app, "js", f), "utf8")).join("\n");
const lib = new Function(src + "\nreturn { parseSpk, electrodeStats, nbWell, nbIsinValues, nbIsinHist, nbValley };")();
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
