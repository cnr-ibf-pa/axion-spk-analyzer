"""Minimal NeuroExplorer .nex (v1) reader, per the published NexFileHeader / NexVarHeader layout."""
import struct, numpy as np
TYPES = {0: "neuron", 1: "event", 2: "interval", 3: "waveform", 4: "popvector", 5: "continuous", 6: "marker"}
def read_nex(path):
    b = open(path, "rb").read()
    magic, ver = struct.unpack_from("<4si", b, 0)
    assert magic == b"NEX1", magic
    comment = b[8:264].split(b"\0")[0].decode("latin-1")
    freq, tbeg, tend, nvar = struct.unpack_from("<diii", b, 264)
    hdr = dict(version=ver, comment=comment, freq=freq, tbeg=tbeg, tend=tend, nvar=nvar)
    vars_ = []
    for k in range(nvar):
        o = 544 + 208 * k
        typ, vver = struct.unpack_from("<ii", b, o)
        name = b[o + 8:o + 72].split(b"\0")[0].decode("latin-1")
        off, cnt, wire, unit, gain, filt = struct.unpack_from("<6i", b, o + 72)
        xpos, ypos, wfreq, adtomv = struct.unpack_from("<4d", b, o + 96)
        npw, nmark, mlen = struct.unpack_from("<3i", b, o + 128)
        mvoff = struct.unpack_from("<d", b, o + 140)[0]
        pre = struct.unpack_from("<d", b, o + 148)[0] if vver >= 104 else None
        v = dict(type=TYPES.get(typ, typ), version=vver, name=name, offset=off, count=cnt, wire=wire, unit=unit,
                 wfreq=wfreq, adtomv=adtomv, npw=npw, mvoffset=mvoff, prethr=pre)
        if typ in (0, 1, 3):
            v["ts"] = np.frombuffer(b, "<i4", cnt, off).copy()
        if typ == 3:
            v["wf_raw"] = np.frombuffer(b, "<i2", cnt * npw, off + 4 * cnt).reshape(cnt, npw).copy()
        vars_.append(v)
    return hdr, vars_
