// Compute the firing statistics with the page's own code (js/spk_reader.js + js/analysis.js) in Node,
// for compare_firing_stats.py.   node dump_page_stats.js .. recording.spk out.json [t0 t1]   (.. = the app folder)
const fs = require("fs");
const [appArg, spk, out, t0s, t1s] = process.argv.slice(2);
// The page's analysis code: js/spk_reader.js + js/analysis.js of the app folder (or of the folder of mea_pipeline.html).
const path = require("path");
const app = /\.html?$/i.test(appArg) ? path.dirname(appArg) : appArg;
const src = ["spk_reader.js", "analysis.js"].map(f => fs.readFileSync(path.join(app, "js", f), "utf8")).join("\n");
const lib = new Function(src + "\nreturn { parseSpk, electrodeStats, autocorrelogram, isiHistogram };")();
const b0 = fs.readFileSync(spk), F = lib.parseSpk(b0.buffer.slice(b0.byteOffset, b0.byteOffset + b0.byteLength));
const dur = F.header.Duration || F.times[F.n - 1], a = t0s ? +t0s : 0, T = (t1s ? +t1s : dur) - a;
const prm = { activePerMin: 5, fanoBin: 1, shortIsi: 0.0015, minIsi: 20 };
const per = new Map();
for (let i = 0; i < F.n; i++) { const c = F.ch[i]; if (c < 0) continue; const t = F.times[i]; if (t < a || t > a + T) continue; if (!per.has(c)) per.set(c, []); per.get(c).push(t); }
const res = {};
for (const [c, t] of per) { const s = lib.electrodeStats(Float64Array.from(t), a, T, prm); const acg = lib.autocorrelogram(Float64Array.from(t), 0.1, 0.001); s.acgSum = acg.c.reduce((x, y) => x + y, 0); res[F.channels[c].label] = s; }
fs.writeFileSync(out, JSON.stringify(res, (k, v) => (typeof v === "number" && !Number.isFinite(v)) ? null : v));
console.log(`page statistics for ${Object.keys(res).length} electrodes with spikes -> ${out}`);
