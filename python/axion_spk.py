"""
axion_spk.py - pure-Python (numpy) reader for Axion BioSystems AxIS / AxIS Navigator
spike files (.spk).

The binary layout follows Axion's official MATLAB "AxionFileLoader"
(AxisFile.m, EntryRecord.m, CombinedBlockVectorHeaderEntry.m,
DiscontinuousBlockVectorHeaderEntry.m, SpikeDataSet.m, ChannelArray.m, ...),
re-implemented here so that .spk files can be opened without MATLAB.

Usage
-----
    from axion_spk import read_spk
    spk = read_spk("recording.spk")
    spk.summary()
    spk.spikes            # numpy structured array (one row per spike)
    spk.waveforms_uv      # (n_spikes, n_samples) float32, microvolts
    spk.spike_times_s     # trigger time of each spike, seconds

CLI
---
    python axion_spk.py recording.spk                  # print a summary
    python axion_spk.py recording.spk --csv out.csv    # export spike table
    python axion_spk.py recording.spk --npz out.npz    # export everything
"""
from __future__ import annotations

import argparse
import binascii
import datetime as _dt
import struct
import sys
from dataclasses import dataclass, field

import numpy as np

MAGIC = b"AxionBio"
PRIMARY_HEADER_CRCSIZE = 1018
SUBHEADER_CRCSIZE = 1016
PRIMARY_HEADER_MAXENTRIES = 123
SUBHEADER_MAXENTRIES = 126
EXPECTED_NOTES_LENGTH_FIELD = 600

# EntryRecordID
E_TERMINATE, E_NOTES, E_CHANNELS, E_BV_HEADER, E_BV_DATA, E_BV_HEADER_EXT, E_TAG, E_COMBINED = range(8)
E_SKIP = 0xFF
ENTRY_NAMES = {0: "Terminate", 1: "NotesArray", 2: "ChannelArray", 3: "BlockVectorHeader",
               4: "BlockVectorData", 5: "BlockVectorHeaderExtension", 6: "Tag",
               7: "CombinedBlockVectorHeader", 0xFF: "Skip"}

DATA_TYPES = {0: "Raw_v1", 1: "Spike_v1", 2: "NamedContinuousData"}
SAMPLE_TYPES = {0: ("Short", "<i2"), 1: ("Int", "<i4"), 2: ("Float", "<f4"), 3: ("Double", "<f8")}
TAG_TYPES = {0: "Deleted", 1: "WellTreatment", 2: "UserAnnotation", 3: "SystemAnnotation",
             4: "DataLossEvent", 5: "StimulationEvent", 6: "StimulationChannelGroup",
             7: "StimulationWaveform", 8: "CalibrationTag", 9: "StimulationLedGroup",
             10: "DoseEvent", 11: "StringDictionaryKeyPair", 12: "LeapInductionEvent",
             13: "ViabilityImpedanceEvent"}

# PlateTypes.m  ->  (name, well rows, well cols, electrode cols, electrode rows)
PLATE_TYPES = {
    0x0000000: ("Empty", 0, 0, 0, 0),
    0x0400000: ("LinearSingleWell64", 1, 1, 8, 8),
    0x0400001: ("P200D30S", 1, 1, 8, 8),
    0x0800000: ("MEACreatorKit", 1, 2, 8, 8),
    0x1800000: ("LinearSixWell", 2, 3, 8, 8),
    0x1800001: ("TwentyFourWell", 4, 6, 4, 4),
    0x1800002: ("TwentyFourWellCircuit", 4, 6, 4, 4),
    0x1800003: ("SixWell", 2, 3, 8, 8),
    0x1800004: ("TwentyFourWellLumos", 4, 6, 4, 4),
    0x3000000: ("LinearTwelveWell", 3, 4, 8, 8),
    0x3000001: ("TwelveWell", 3, 4, 8, 8),
    0x3000002: ("FortyEightWell", 6, 8, 4, 4),
    0x3000003: ("NinetySixWell", 8, 12, 3, 3),
    0x3000004: ("FortyEightWellTransparent", 6, 8, 4, 4),
    0x3000006: ("FortyEightWellLumos", 6, 8, 4, 4),
    0x3000007: ("FortyEightWellEStimPlus", 6, 8, 4, 4),
    0x3000008: ("NinetySixWellCircuit", 8, 12, 3, 3),
    0x3000009: ("NinetySixWellTransparent", 8, 12, 3, 3),
    0x300000A: ("NinetySixWellLumos", 8, 12, 3, 3),
    0x300000B: ("FortyEightWellCircuit", 6, 8, 4, 4),
    0x300000C: ("FortyEightWellAccuSpot", 6, 8, 4, 4),
    0x3000016: ("TwelveWellTransparent", 3, 4, 8, 8),
}

SPIKE_HEADER_SIZE = 30  # Spike_v1.LOADED_HEADER_SIZE


class AxisFormatError(Exception):
    pass


def crc32_axis(data: bytes) -> int:
    """Axion uses the standard reflected CRC-32 (poly 0xEDB88320, seed 0xFFFFFFFF)."""
    return binascii.crc32(data) & 0xFFFFFFFF


@dataclass
class EntryRecord:
    type: int
    length: float  # may be inf

    @property
    def name(self):
        return ENTRY_NAMES.get(self.type, f"Unknown({self.type})")

    @staticmethod
    def from_u64(values):
        out = []
        for v in values:
            v = int(v)
            t = v >> 56
            hi = (v >> 32) & 0xFFFFFF
            lo = v & 0xFFFFFFFF
            length = float("inf") if (hi == 0xFFFFFF and lo == 0xFFFFFFFF) else (hi << 32) | lo
            if t not in ENTRY_NAMES:
                t = E_SKIP
            out.append(EntryRecord(t, length))
        return out


class _Reader:
    def __init__(self, buf: bytes):
        self.buf = buf
        self.pos = 0

    def seek(self, p):
        self.pos = int(p)

    def read(self, n):
        b = self.buf[self.pos:self.pos + n]
        if len(b) != n:
            raise AxisFormatError(f"unexpected EOF at {self.pos}")
        self.pos += n
        return b

    def unpack(self, fmt):
        s = struct.calcsize(fmt)
        return struct.unpack(fmt, self.read(s))

    def u8(self): return self.unpack("<B")[0]
    def u16(self): return self.unpack("<H")[0]
    def u32(self): return self.unpack("<I")[0]
    def i32(self): return self.unpack("<i")[0]
    def i64(self): return self.unpack("<q")[0]
    def u64(self): return self.unpack("<Q")[0]
    def f64(self): return self.unpack("<d")[0]

    def datetime(self):
        y, mo, d, h, mi, s, ms = self.unpack("<7H")
        try:
            return _dt.datetime(y, mo, d, h, mi, s, ms * 1000)
        except ValueError:
            return None

    def lstring(self):
        n = self.i32()
        return self.read(n).decode("utf-8", errors="replace")


def _cstr(b: bytes) -> str:
    return b.split(b"\x00", 1)[0].decode("latin-1").replace("\r", "").strip()


@dataclass
class SpikeFile:
    path: str
    primary_data_type: int = 0
    version: tuple = (0, 0)
    notes: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    plate_type: int = 0
    channels: np.ndarray | None = None      # structured array, ChannelArray
    header: dict = field(default_factory=dict)
    channel_ids: list = field(default_factory=list)
    dataset_name: str = ""
    spikes: np.ndarray | None = None        # structured array, one row per spike
    tags: list = field(default_factory=list)
    entries: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    # ---- derived ---------------------------------------------------------
    @property
    def sampling_frequency(self):
        return self.header.get("SamplingFrequency")

    @property
    def voltage_scale(self):
        return self.header.get("VoltageScale")

    @property
    def plate(self):
        return PLATE_TYPES.get(self.plate_type, (f"Unknown(0x{self.plate_type:07X})", 0, 0, 0, 0))

    @property
    def plate_name(self):
        """Plate as AxIS names it (the 'Plate Type' recording setting, e.g. 'CytoView MEA 6'),
        else the plate-type code name from PlateTypes.m."""
        for line in (self.header.get("Description") or "").splitlines():
            k, _, v = line.partition(",")
            if k.strip() == "Plate Type" and v.strip():
                return v.strip()
        return self.plate[0]

    @property
    def spike_times_s(self):
        s = self.spikes
        return (s["startingSample"].astype(np.float64) + s["triggerSample"]) / self.sampling_frequency

    @property
    def waveforms_uv(self):
        return self.spikes["data"].astype(np.float32) * np.float32(self.voltage_scale * 1e6)

    def electrode_labels(self):
        """Axion-style label 'A1_12' (well row letter + well col _ electrode col + row) per spike."""
        lut = self._channel_lut()
        idx = lut[self.spikes["chip"].astype(int) * 256 + self.spikes["channel"].astype(int)]
        ch = self.channels
        return np.array([self._label(ch[i]) if i >= 0 else "?" for i in idx])

    def sort(self, electrodes=None, **params):
        """SAMS-style automatic spike sorting (see axion_sort.py); same results as the viewer."""
        from axion_sort import sort_spikes
        return sort_spikes(self, electrodes=electrodes, **params)

    def _channel_lut(self):
        lut = np.full(256 * 256, -1, dtype=np.int64)
        for i, c in enumerate(self.channels):
            lut[int(c["ChannelAchk"]) * 256 + int(c["ChannelIndex"])] = i
        return lut

    @staticmethod
    def _label(c):
        return f"{chr(64 + int(c['WellRow']))}{int(c['WellColumn'])}_{int(c['ElectrodeColumn'])}{int(c['ElectrodeRow'])}"

    def summary(self, out=sys.stdout):
        h = self.header
        p = self.plate
        w = out.write
        w(f"File            : {self.path}\n")
        w(f"Header version  : {self.version[0]}.{self.version[1]}\n")
        for k, v in self.metadata.items():
            w(f"{k:<16}: {v}\n")
        w(f"Plate           : {self.plate_name} ({p[1]}x{p[2]} wells, {p[3]}x{p[4]} electrodes)\n")
        w(f"Channels        : {len(self.channels) if self.channels is not None else 0}\n")
        w(f"Data set        : {h.get('Name')} / {self.dataset_name} ({DATA_TYPES.get(h.get('DataType'))})\n")
        w(f"Sampling freq   : {h.get('SamplingFrequency')} Hz\n")
        w(f"Voltage scale   : {h.get('VoltageScale')} V/LSB\n")
        w(f"Samples/spike   : {h.get('NumSamplesPerBlock')}\n")
        w(f"Experiment start: {h.get('ExperimentStartTime')}\n")
        w(f"Duration        : {h.get('Duration')} s\n")
        w(f"Spikes          : {len(self.spikes) if self.spikes is not None else 0}\n")
        for wmsg in self.warnings:
            w(f"warning: {wmsg}\n")


def read_spk(path: str, check_crc: bool = True) -> SpikeFile:
    with open(path, "rb") as fh:
        buf = fh.read()
    r = _Reader(buf)
    sf = SpikeFile(path=path)

    if r.read(8) != MAGIC:
        raise AxisFormatError("Not a modern Axis file (missing 'AxionBio' magic word)")
    sf.primary_data_type = r.u16()
    vmaj, vmin = r.u16(), r.u16()
    sf.version = (vmaj, vmin)
    _notes_start = r.u64()
    notes_len = r.u32()
    if notes_len != EXPECTED_NOTES_LENGTH_FIELD:
        raise AxisFormatError("Incorrect legacy notes length field")
    if vmaj != 1:
        raise AxisFormatError(f"Unsupported header version {vmaj}.{vmin} (only 1.x implemented)")
    entries_start = r.i64()
    records = EntryRecord.from_u64(r.unpack(f"<{PRIMARY_HEADER_MAXENTRIES}Q"))
    if check_crc:
        read_crc = struct.unpack_from("<I", buf, PRIMARY_HEADER_CRCSIZE)[0]
        if read_crc != crc32_axis(buf[:PRIMARY_HEADER_CRCSIZE]):
            raise AxisFormatError("File header checksum was incorrect")

    r.seek(entries_start)
    combined = None
    data_entry = None
    terminated = False
    while not terminated:
        for rec in records:
            start = r.pos
            sf.entries.append({"type": rec.name, "start": start, "length": rec.length})
            if rec.type == E_TERMINATE:
                terminated = True
                break
            elif rec.type == E_NOTES:
                for _ in range(int(rec.length) // 618):
                    nb = r.read(618)
                    sf.notes.append({
                        "Investigator": _cstr(nb[0:50]),
                        "RecordingName": _cstr(nb[50:100]),
                        "Description": _cstr(nb[100:600]),
                        "Revision": struct.unpack_from("<I", nb, 600)[0],
                    })
            elif rec.type == E_CHANNELS:
                sf.plate_type = r.u32()
                n = r.u32()
                dt = np.dtype([("WellColumn", "u1"), ("WellRow", "u1"), ("ElectrodeColumn", "u1"),
                               ("ElectrodeRow", "u1"), ("ChannelAchk", "u1"), ("ChannelIndex", "u1"),
                               ("AuxData", "<u2")])
                sf.channels = np.frombuffer(r.read(8 * n), dtype=dt).copy()
            elif rec.type == E_COMBINED:
                combined = _read_combined_header(r, sf, check_crc)
            elif rec.type == E_BV_DATA:
                data_entry = (start, rec.length)
                if rec.length == float("inf"):
                    r.seek(len(buf))
                else:
                    r.seek(start + rec.length)
            elif rec.type == E_TAG:
                sf.tags.append(_read_tag(r, start, int(rec.length)))
                r.seek(start + rec.length)
            elif rec.type in (E_BV_HEADER, E_BV_HEADER_EXT):
                sf.warnings.append(f"Legacy {rec.name} entry skipped (not used by modern .spk files)")
                r.seek(start + rec.length)
            else:
                r.seek(start + rec.length)
        if not terminated:
            if r.pos >= len(buf):
                break
            sub_start = r.pos
            if r.read(8) != MAGIC:
                raise AxisFormatError(f"Bad sub header magic number at {sub_start}")
            records = EntryRecord.from_u64(r.unpack(f"<{SUBHEADER_MAXENTRIES}Q"))
            read_crc = struct.unpack_from("<I", buf, sub_start + SUBHEADER_CRCSIZE)[0]
            if check_crc and read_crc != crc32_axis(buf[sub_start:sub_start + SUBHEADER_CRCSIZE]):
                raise AxisFormatError("Bad sub header checksum")
            r.seek(sub_start + SUBHEADER_CRCSIZE + 4 + 4)

    # metadata: key/value tags (newest revision wins) else legacy notes
    kv = {}
    for t in sorted((t for t in sf.tags if t.get("Key") is not None), key=lambda t: t["Revision"]):
        kv[t["Key"]] = t["Value"]
    if kv:
        sf.metadata = kv
    elif sf.notes:
        n = sorted(sf.notes, key=lambda x: x["Revision"])[-1]
        sf.metadata = {"RecordingName": n["RecordingName"], "Investigator": n["Investigator"],
                       "Description": n["Description"]}

    if combined is None:
        raise AxisFormatError("No CombinedBlockVectorHeader found - is this a spike (.spk) file?")
    if combined["DataType"] != 1:
        raise AxisFormatError(f"Data set is {DATA_TYPES.get(combined['DataType'])}, not Spike_v1")
    sf.spikes = _map_spikes(buf, combined)
    return sf


def _read_combined_header(r: _Reader, sf: SpikeFile, check_crc: bool) -> dict:
    p0 = r.pos
    h = {}
    h["VersionMajor"], h["VersionMinor"], h["DataType"], h["SampleType"] = r.unpack("<4H")
    h["SamplingFrequency"], h["VoltageScale"] = r.unpack("<2d")
    (h["NumChannelsPerBlock"], h["NumDataSetsPerBlock"], h["NumSamplesPerBlock"],
     h["VectorHeaderSize"], h["BlockHeaderSize"]) = r.unpack("<5I")
    h["BlockVectorStartTime"] = r.datetime()
    h["ExperimentStartTime"] = r.datetime()
    h["AddedDate"] = r.datetime()
    h["ModifiedDate"] = r.datetime()
    h["Duration"] = r.f64() if (h["VersionMajor"] > 1 or h["VersionMinor"] >= 1) else None
    h["Name"] = r.read(r.i32()).decode("utf-8", errors="replace")
    h["Description"] = r.read(r.i32()).decode("utf-8", errors="replace")
    h["DataRegionStart"] = r.i64()
    h["DataRegionLength"] = r.i64()
    crc_calc = crc32_axis(r.buf[p0:r.pos])
    if check_crc and r.u32() != crc_calc:
        raise AxisFormatError("BlockVectorMetaData checksum was incorrect")
    elif not check_crc:
        r.u32()
    if h["DataType"] == 1:  # Spike_v1 -> DiscontinuousBlockVectorHeaderEntry
        p1 = r.pos
        n = r.u32()
        ids = r.unpack(f"<{n}H")
        sf.channel_ids = [(i >> 8 & 0xFF, i & 0xFF) for i in ids]
        sf.dataset_name = r.read(r.i32()).decode("latin-1")
        crc_calc = crc32_axis(r.buf[p1:r.pos])
        if r.u32() != crc_calc and check_crc:
            raise AxisFormatError("Discontinuous header checksum was incorrect")
    h["SampleTypeName"] = SAMPLE_TYPES.get(h["SampleType"], ("?",))[0]
    sf.header = h
    return h


def _map_spikes(buf: bytes, h: dict) -> np.ndarray:
    if h["NumChannelsPerBlock"] != 1 or h["NumDataSetsPerBlock"] != 1:
        raise AxisFormatError("Invalid header for SPIKE file (channels/datasets per block != 1)")
    reserved = h["BlockHeaderSize"] - SPIKE_HEADER_SIZE
    if reserved < 0:
        raise AxisFormatError("Block header size too small")
    ns = h["NumSamplesPerBlock"]
    fields = [("startingSample", "<i8"), ("channel", "u1"), ("chip", "u1"), ("triggerSample", "<i4"),
              ("standardDeviation", "<f8"), ("thresholdMultiplier", "<f8")]
    if reserved:
        fields.append(("reserved", "u1", (reserved,)))
    fields.append(("data", "<i2", (ns,)))
    dt = np.dtype(fields)
    n = h["DataRegionLength"] // dt.itemsize
    return np.frombuffer(buf, dtype=dt, count=n, offset=h["DataRegionStart"])


def _read_tag(r: _Reader, start: int, length: int) -> dict:
    t = {"Start": start}
    ttype = r.u16()
    t["Type"] = TAG_TYPES.get(ttype, f"Unknown({ttype})")
    t["CreationDate"] = r.datetime()
    g = r.read(16)
    t["Guid"] = g.hex()
    t["Revision"] = r.u32()
    try:
        if ttype == 11:
            t["Key"] = r.lstring()
            t["Value"] = r.lstring()
        elif ttype == 1:
            t["WellColumn"], t["WellRow"] = r.u8(), r.u8()
            r.u8(); r.u8()
            wt = r.u8()
            t["IsOn"], t["IsControl"] = bool(wt & 1), bool(wt & 2)
            if start + length - r.pos >= 12:
                t["Color"] = "#%02x%02x%02x" % tuple(r.read(3))
                t["TreatmentWhat"] = r.lstring()
                t["AdditionalInformation"] = r.lstring()
    except AxisFormatError:
        pass
    return t


def _main(argv=None):
    ap = argparse.ArgumentParser(description="Read an Axion BioSystems .spk spike file")
    ap.add_argument("spk")
    ap.add_argument("--csv", help="write spike table (time, electrode, amplitude...) to CSV")
    ap.add_argument("--npz", help="write spikes, waveforms (uV) and channel map to .npz")
    ap.add_argument("--no-crc", action="store_true", help="skip checksum verification")
    a = ap.parse_args(argv)
    sf = read_spk(a.spk, check_crc=not a.no_crc)
    sf.summary()
    if a.csv:
        labels = sf.electrode_labels()
        wf = sf.waveforms_uv
        trig = sf.spikes["triggerSample"]
        amp = wf[np.arange(len(wf)), np.clip(trig, 0, wf.shape[1] - 1)]
        with open(a.csv, "w") as fh:
            fh.write("time_s,electrode,trigger_amplitude_uV,min_uV,max_uV,std_uV,threshold_multiplier\n")
            for t, lab, am, mn, mx, sd, th in zip(sf.spike_times_s, labels, amp, wf.min(1), wf.max(1),
                                                  sf.spikes["standardDeviation"] * 1e6,
                                                  sf.spikes["thresholdMultiplier"]):
                fh.write(f"{t:.6f},{lab},{am:.3f},{mn:.3f},{mx:.3f},{sd:.3f},{th:.3f}\n")
        print(f"wrote {a.csv}")
    if a.npz:
        np.savez_compressed(a.npz, times_s=sf.spike_times_s, electrode=sf.electrode_labels(),
                            waveforms_uV=sf.waveforms_uv, channels=sf.channels,
                            sampling_frequency=sf.sampling_frequency)
        print(f"wrote {a.npz}")


if __name__ == "__main__":
    _main()
