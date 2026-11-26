// Run the web app's own .spk parser (js/spk_reader.js) in Node and dump every spike field,
// for compare_browser_parser.py. The first argument is the app folder (or a file in it).
//   node dump_browser_parser.js .. recording.spk /tmp/out
const fs = require("fs"), path = require("path");
const [appArg, inp, out] = process.argv.slice(2);
const app = fs.statSync(appArg).isDirectory() ? appArg : path.dirname(appArg);
const parseSpk = new Function(fs.readFileSync(path.join(app, "js", "spk_reader.js"), "utf8") + "\nreturn parseSpk;")();
const buf0 = fs.readFileSync(inp), buf = buf0.buffer.slice(buf0.byteOffset, buf0.byteOffset + buf0.byteLength);
const F = parseSpk(buf), ns = F.header.NumSamplesPerBlock;
const wf = new Int16Array(F.n * ns);
for (let i = 0; i < F.n; i++) { const o = F.waveOffset(i); for (let k = 0; k < ns; k++) wf[i * ns + k] = F.dv.getInt16(o + 2 * k, true); }
const w = (ext, arr) => fs.writeFileSync(out + ext, Buffer.from(arr.buffer, arr.byteOffset, arr.byteLength));
w(".times", Float64Array.from(F.times)); w(".sd", Float64Array.from(F.sd)); w(".thr", Float64Array.from(F.thr));
w(".trig", Int32Array.from(F.trig)); w(".wf", wf);
fs.writeFileSync(out + ".lab", Array.from(F.ch, c => c >= 0 ? F.channels[c].label : "?").join("\n"));
fs.writeFileSync(out + ".hdr.json", JSON.stringify({ n: F.n, ns, fs: F.fs, scale: F.header.VoltageScale, duration: F.header.Duration,
  plateType: F.plateType, nch: F.channels.length, sortedOnLoad: F.warnings.some(x => /sorted on load/.test(x)) }));
console.log(`browser parser: ${F.n} spikes`);
