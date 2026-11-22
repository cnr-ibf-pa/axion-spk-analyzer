// Compute the firing statistics with the page's own code (parser + STATS block of mea_pipeline.html) in Node,
// for compare_firing_stats.py.   node dump_page_stats.js ../mea_pipeline.html recording.spk out.json [t0 t1]
const fs = require("fs");
const [html, spk, out, t0s, t1s] = process.argv.slice(2);
const page = fs.readFileSync(html, "utf8");
const pa = page.indexOf("const PLATES = {"), pb = page.indexOf("\nfunction parseSpk("), pe = page.indexOf("\n}\n", pb) + 3;
const sa = page.indexOf("/* STATS-BEGIN"), se = page.indexOf("/* STATS-END */");
if (pa < 0 || pb < 0 || sa < 0 || se < 0) throw new Error("parser or STATS block not found");
const lib = new Function(page.slice(pa, pe) + "\n" + page.slice(sa, se) + "\nreturn { parseSpk, electrodeStats, autocorrelogram, isiHistogram };")();
const b0 = fs.readFileSync(spk), F = lib.parseSpk(b0.buffer.slice(b0.byteOffset, b0.byteOffset + b0.byteLength));
const dur = F.header.Duration || F.times[F.n - 1], a = t0s ? +t0s : 0, T = (t1s ? +t1s : dur) - a;
const prm = { activePerMin: 5, fanoBin: 1, shortIsi: 0.0015, minIsi: 20 };
const per = new Map();
for (let i = 0; i < F.n; i++) { const c = F.ch[i]; if (c < 0) continue; const t = F.times[i]; if (t < a || t > a + T) continue; if (!per.has(c)) per.set(c, []); per.get(c).push(t); }
const res = {};
for (const [c, t] of per) { const s = lib.electrodeStats(Float64Array.from(t), a, T, prm); const acg = lib.autocorrelogram(Float64Array.from(t), 0.1, 0.001); s.acgSum = acg.c.reduce((x, y) => x + y, 0); res[F.channels[c].label] = s; }
fs.writeFileSync(out, JSON.stringify(res, (k, v) => (typeof v === "number" && !Number.isFinite(v)) ? null : v));
console.log(`page statistics for ${Object.keys(res).length} electrodes with spikes -> ${out}`);
