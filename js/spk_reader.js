"use strict";
/* =====================================================================
 * Axion .spk reader, verified against Axion's Data Export Tool and the
 * Python reader (verification/README.md). Originally the reader of the former
 * spk_viewer.html; the detector noise SD and threshold multiplier are kept in
 * double precision here.
 * ===================================================================== */
const PLATES = {
  0x0400000:"LinearSingleWell64",0x0400001:"P200D30S",0x0800000:"MEACreatorKit",0x1800000:"LinearSixWell",
  0x1800001:"TwentyFourWell",0x1800002:"TwentyFourWellCircuit",0x1800003:"SixWell",0x1800004:"TwentyFourWellLumos",
  0x1800005:"TwentyFourWellOptiClear",0x3000000:"LinearTwelveWell",0x3000001:"TwelveWell",0x3000002:"FortyEightWell",
  0x3000003:"NinetySixWell",0x3000004:"FortyEightWellTransparent",0x3000006:"FortyEightWellLumos",
  0x3000007:"FortyEightWellEStimPlus",0x3000008:"NinetySixWellCircuit",0x3000009:"NinetySixWellTransparent",
  0x300000A:"NinetySixWellLumos",0x300000B:"FortyEightWellCircuit",0x300000C:"FortyEightWellAccuSpot",
  0x300000D:"FortyEightWellOptiClear",0x300000E:"NinetySixWellOptiClear",0x3000015:"FortyEightWellOrganoid",
  0x3000016:"TwelveWellTransparent"};
const DATA_TYPES = ["Raw_v1","Spike_v1","NamedContinuousData"];
const CRC_TABLE = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
function crc32(u8, a, b) { let c = 0xFFFFFFFF; for (let i = a; i < b; i++) c = CRC_TABLE[(c ^ u8[i]) & 0xFF] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }
const utf8 = new TextDecoder("utf-8");

function parseSpk(buf) {
  const dv = new DataView(buf), u8 = new Uint8Array(buf), warnings = [];
  let p = 0;
  const rd = {
    u8: () => dv.getUint8(p++),
    u16: () => { const v = dv.getUint16(p, true); p += 2; return v; },
    u32: () => { const v = dv.getUint32(p, true); p += 4; return v; },
    i32: () => { const v = dv.getInt32(p, true); p += 4; return v; },
    i64: () => { const v = dv.getUint32(p, true) + dv.getInt32(p + 4, true) * 4294967296; p += 8; return v; },
    f64: () => { const v = dv.getFloat64(p, true); p += 8; return v; },
    str: n => { const s = utf8.decode(u8.subarray(p, p + n)); p += n; return s; },
    lstr: () => rd.str(rd.i32()),
    date: () => { const a = []; for (let i = 0; i < 7; i++) a.push(rd.u16()); return a[0] ? new Date(a[0], a[1]-1, a[2], a[3], a[4], a[5], a[6]) : null; },
  };
  const records = n => { const out = []; for (let i = 0; i < n; i++) {
      const lo = dv.getUint32(p, true), hiw = dv.getUint32(p + 4, true); p += 8;
      const type = hiw >>> 24, hi = hiw & 0xFFFFFF;
      out.push({ type: type <= 7 ? type : 0xFF, len: (hi === 0xFFFFFF && lo === 0xFFFFFFFF) ? Infinity : hi * 4294967296 + lo });
    } return out; };

  if (rd.str(8) !== "AxionBio") throw new Error("Not an Axion file: missing 'AxionBio' magic word.");
  rd.u16(); const vmaj = rd.u16(), vmin = rd.u16(); p += 8;
  if (rd.u32() !== 600) throw new Error("Incorrect legacy notes length field.");
  if (vmaj !== 1) throw new Error(`Unsupported header version ${vmaj}.${vmin} (only 1.x is implemented).`);
  const entriesStart = rd.i64();
  let recs = records(123);
  if (dv.getUint32(1018, true) !== crc32(u8, 0, 1018)) throw new Error("File header checksum is incorrect.");

  const F = { version: `${vmaj}.${vmin}`, notes: [], tags: [], channels: [], channelIds: [], header: null, warnings };
  p = entriesStart;
  let done = false;
  while (!done) {
    for (const r of recs) {
      const start = p;
      if (r.type === 0) { done = true; break; }
      if (r.type === 1) {            // NotesArray
        for (let i = 0; i < Math.floor(r.len / 618); i++) {
          const b = p, cs = (o, n) => utf8.decode(u8.subarray(b + o, b + o + n)).split("\0")[0].replace(/\r/g, "").trim();
          F.notes.push({ Investigator: cs(0, 50), RecordingName: cs(50, 50), Description: cs(100, 500), Revision: dv.getUint32(b + 600, true) });
          p += 618;
        }
      } else if (r.type === 2) {     // ChannelArray
        F.plateType = rd.u32(); const n = rd.u32();
        for (let i = 0; i < n; i++) {
          const c = { wc: rd.u8(), wr: rd.u8(), ec: rd.u8(), er: rd.u8(), achk: rd.u8(), idx: rd.u8() }; rd.u16();
          c.well = String.fromCharCode(64 + c.wr) + c.wc; c.label = `${c.well}_${c.ec}${c.er}`; c.i = i;
          F.channels.push(c);
        }
      } else if (r.type === 7) {     // CombinedBlockVectorHeader
        const h = {}; const p0 = p;
        h.VersionMajor = rd.u16(); h.VersionMinor = rd.u16(); h.DataType = rd.u16(); h.SampleType = rd.u16();
        h.SamplingFrequency = rd.f64(); h.VoltageScale = rd.f64();
        h.NumChannelsPerBlock = rd.u32(); h.NumDataSetsPerBlock = rd.u32(); h.NumSamplesPerBlock = rd.u32();
        h.VectorHeaderSize = rd.u32(); h.BlockHeaderSize = rd.u32();
        h.BlockVectorStartTime = rd.date(); h.ExperimentStartTime = rd.date(); h.AddedDate = rd.date(); h.ModifiedDate = rd.date();
        h.Duration = (h.VersionMajor > 1 || h.VersionMinor >= 1) ? rd.f64() : null;
        h.Name = rd.lstr(); h.Description = rd.lstr();
        h.DataRegionStart = rd.i64(); h.DataRegionLength = rd.i64();
        if (rd.u32() !== crc32(u8, p0, p - 4)) throw new Error("Block-vector header checksum is incorrect.");
        if (h.DataType === 1) {
          const p1 = p, n = rd.u32();
          for (let i = 0; i < n; i++) { const id = rd.u16(); F.channelIds.push([(id >> 8) & 0xFF, id & 0xFF]); }
          h.DataSetName = rd.lstr();
          if (rd.u32() !== crc32(u8, p1, p - 4)) throw new Error("Spike header checksum is incorrect.");
        }
        if (!F.header) F.header = h; else warnings.push("More than one block-vector set; only the first is shown.");
      } else if (r.type === 6) {     // Tag
        const t = { type: rd.u16() }; rd.date(); p += 16; t.rev = rd.u32();
        try {
          if (t.type === 11) { t.key = rd.lstr(); t.value = rd.lstr(); }
          else if (t.type === 1) {
            t.wc = rd.u8(); t.wr = rd.u8(); p += 2; const wt = rd.u8(); t.on = !!(wt & 1); t.control = !!(wt & 2);
            if (start + r.len - p >= 12) { t.color = "#" + [rd.u8(), rd.u8(), rd.u8()].map(x => x.toString(16).padStart(2, "0")).join(""); t.what = rd.lstr(); t.info = rd.lstr(); }
          }
        } catch (e) { /* malformed tag payload: ignore */ }
        F.tags.push(t);
        p = start + r.len;
      } else {                       // BlockVectorData and anything else: skip
        p = r.len === Infinity ? buf.byteLength : start + r.len;
      }
    }
    if (!done) {
      if (p >= buf.byteLength) break;
      const s = p;
      if (rd.str(8) !== "AxionBio") throw new Error("Bad sub-header magic number at byte " + s);
      recs = records(126);
      if (dv.getUint32(s + 1016, true) !== crc32(u8, s, s + 1016)) throw new Error("Bad sub-header checksum.");
      p = s + 1024;
    }
  }
  const h = F.header;
  if (!h) throw new Error("No block-vector header found. Is this a spike (.spk) file?");
  if (h.DataType !== 1) throw new Error(`This file holds ${DATA_TYPES[h.DataType] || h.DataType} data, not spikes (Spike_v1).`);
  if (h.NumChannelsPerBlock !== 1 || h.NumDataSetsPerBlock !== 1 || h.BlockHeaderSize < 30) throw new Error("Invalid Spike_v1 block layout.");

  // metadata: key/value tags (latest revision) else notes
  F.meta = {};
  F.tags.filter(t => t.key !== undefined).sort((a, b) => a.rev - b.rev).forEach(t => F.meta[t.key] = t.value);
  if (!Object.keys(F.meta).length && F.notes.length) { const n = F.notes[F.notes.length - 1]; F.meta = { RecordingName: n.RecordingName, Investigator: n.Investigator, Description: n.Description }; }
  F.wells = {}; F.tags.filter(t => t.type === 1).forEach(t => F.wells[String.fromCharCode(64 + t.wr) + t.wc] = t);

  // spikes
  const lut = new Int32Array(65536).fill(-1);
  F.channels.forEach(c => lut[c.achk * 256 + c.idx] = c.i);
  const bs = h.BlockHeaderSize + 2 * h.NumSamplesPerBlock, n = Math.floor(h.DataRegionLength / bs), fs = h.SamplingFrequency;
  const times = new Float64Array(n), ch = new Int32Array(n), trig = new Int32Array(n), sd = new Float64Array(n), thr = new Float64Array(n);
  let o = h.DataRegionStart;
  for (let i = 0; i < n; i++, o += bs) {
    const start = dv.getUint32(o, true) + dv.getInt32(o + 4, true) * 4294967296;
    trig[i] = dv.getInt32(o + 10, true);
    times[i] = (start + trig[i]) / fs;
    ch[i] = lut[dv.getUint8(o + 9) * 256 + dv.getUint8(o + 8)];
    sd[i] = dv.getFloat64(o + 14, true); thr[i] = dv.getFloat64(o + 22, true);
  }
  // the page (raster, replay) relies on time-ordered spikes; sort if the file is not
  let order = null;
  for (let i = 1; i < n; i++) if (times[i] < times[i - 1]) { order = Int32Array.from({ length: n }, (_, k) => k).sort((a, b) => times[a] - times[b]); break; }
  const perm = a => { if (!order) return a; const b = new a.constructor(n); for (let k = 0; k < n; k++) b[k] = a[order[k]]; return b; };
  if (order) warnings.push("Spikes were not stored in time order; they were sorted on load.");
  Object.assign(F, { dv, n, times: perm(times), ch: perm(ch), trig: perm(trig), sd: perm(sd), thr: perm(thr), blockSize: bs, fs });
  F.waveOffset = i => h.DataRegionStart + (order ? order[i] : i) * bs + h.BlockHeaderSize;
  return F;
}
