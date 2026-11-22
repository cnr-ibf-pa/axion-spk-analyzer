// Run the browser page's own .spk parser (parseSpk, extracted verbatim from the HTML) in Node
// and dump every spike field, for compare_browser_parser.py.
//   node dump_browser_parser.js ../spk_viewer.html recording.spk /tmp/out
const fs = require("fs");
const [html, inp, out] = process.argv.slice(2);
const page = fs.readFileSync(html, "utf8");
const a = page.indexOf("const PLATES = {"), b = page.indexOf("\nfunction parseSpk(");
const end = page.indexOf("\n}\n", b) + 3;                       // end of parseSpk
if (a < 0 || b < 0) throw new Error("parseSpk not found in " + html);
const parseSpk = new Function(page.slice(a, end) + "\nreturn parseSpk;")();
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
