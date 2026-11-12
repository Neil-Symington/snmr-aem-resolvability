#!/usr/bin/env python3
"""aseg_to_geoparquet.py -- an ASEG-GDF2 AEM delivery to GeoParquet.

GA ships AEM as three files sharing one stem: a `.dfn` declaring the record
layout, a fixed-width `.dat` holding the records, and a `.prj` or `.des` with
the projection and survey notes. This reads the `.dfn`, streams the `.dat`, and
writes one GeoParquet file with a WKB point geometry column.

    python3 util/aseg_to_geoparquet.py <stem-or-.dat> -o out.parquet
    python3 util/aseg_to_geoparquet.py <stem> --limit 1000 --schema
    python3 util/aseg_to_geoparquet.py --selftest

Why this exists: the deliveries are 2.8 GB and 4.9 GB, so every question about
them ("which line passes near this site", "what is the gate structure") turned
into a bespoke seek-reading script, and the nav GeoParquet the field workflow
once consumed had been produced by hand outside the repository. This is the step
that was missing. It is survey-data preparation, not part of any inversion
stage, which is why it sits in util/ rather than tools/ and why it outlived the
field workflow it was first written beside.

WHAT IT GETS RIGHT, because each was checked against both real deliveries:

  * Array fields. `LM_X:18F15.5` is EIGHTEEN columns, not one. The two
    deliveries carry 18/23 and 26/38 gate counts, and most of the record is
    these: 195 of the UDF file's columns come from 39 DEFN lines.
  * The record-type prefix. The 2017 Ord-Keep file declares `RT:A4` and every
    data line begins with "DATA"; the 2022 UDF file declares no RT field and
    begins with the first data value. Widths then sum to exactly the line
    length in both: 4184 and 2788. A tool that assumes either shape silently
    mis-slices the other, and every column after the offset is wrong while
    still parsing as a number -- which is why this asserts the sum.
  * Per-field NULL. The older file declares a different NULL for nearly every
    field (-99999.99, -9999999.99999, -999999.999999). They become nulls.
  * Unkeyed descriptions. The 2022 file writes `NAME=Local Date`; the 2017 file
    writes the description with no key at all, and one of them contains a comma
    ("TMI, corrected by IGRF"), so attributes cannot be parsed by splitting on
    commas alone.
  * The CRS is taken from the file, never guessed. The 2022 DFN carries
    `EPSGCODE=7844` on Longitude and `EPSGCODE=7855` on Easting. The 2017 DFN
    carries no EPSG code at all, only "(GDA94)" in free text, so this REFUSES
    to write and prints what it found. Writing a GeoParquet with the wrong CRS
    is worse than not writing one, because it is silent.

Requires numpy and pyarrow. pyproj is optional and only improves the CRS
metadata; without it the geometry column carries an EPSG identifier rather than
full PROJJSON. Nothing on the project's conda env has pyarrow as of Nov 2026:

    python3 -m pip install pyarrow        # numpy is already there
"""
from __future__ import annotations
import argparse, json, os, re, sys
from dataclasses import dataclass

import numpy as np

# ASEG-GDF format token: an optional repeat count, a type letter, a width, and
# for real types an optional decimal count. `18F15.5` -> 18 columns of width 15.
FMT = re.compile(r"^(\d*)([IFEDAifeda])(\d+)(?:\.(\d+))?$")

# Attribute keys ASEG uses. Anything else in an attribute section is treated as
# the free-text description, because the 2017 deliveries write it unkeyed.
KEYS = ("NAME", "UNIT", "UNITS", "NULL", "EPSGCODE", "DATUM", "PROJECTION")

# Record types that are not data. COMM holds the header comment block.
NON_DATA_RT = ("COMM",)


@dataclass
class Field:
    """One DFN field, before array expansion."""
    name: str
    kind: str          # I, F, E, D or A
    width: int
    count: int = 1
    decimals: int | None = None
    unit: str = ""
    null: float | None = None
    epsg: str = ""
    desc: str = ""
    rt: str = ""

    @property
    def span(self) -> int:
        return self.width * self.count

    def columns(self) -> list[str]:
        """Expanded column names. `LM_X` with count 18 -> LM_X_01..LM_X_18."""
        if self.count == 1:
            return [self.name]
        pad = len(str(self.count))
        return [f"{self.name}_{i + 1:0{pad}d}" for i in range(self.count)]


def parse_attrs(sections: list[str]) -> dict:
    """Pull KEY=VALUE pairs out of the colon-separated attribute sections.

    `NULL=-99999.99` and `UNIT=nT,NAME=Local Date` both appear, and so does a
    bare description. NAME and the description take the remainder of their
    section, so a comma inside them survives -- the 2017 file has
    `UNIT=nT,TMI, corrected by IGRF + 48425 nT`, where splitting on commas and
    keeping the last piece would report the description as " corrected by IGRF
    + 48425 nT".
    """
    out: dict[str, str] = {}
    desc: list[str] = []
    for sec in sections:
        rest = sec
        while rest:
            m = re.match(r"\s*([A-Za-z_]+)\s*=\s*", rest)
            key = m.group(1).upper() if m else None
            if key in KEYS:
                rest = rest[m.end():]
                if key in ("NAME",):
                    # takes the remainder: a NAME may contain commas
                    out[key] = rest.strip()
                    rest = ""
                else:
                    val, sep, rest = rest.partition(",")
                    out[key] = val.strip()
            else:
                # unkeyed: the whole remaining section is the description
                desc.append(rest.strip())
                rest = ""
    if desc:
        out.setdefault("DESC", ", ".join(d for d in desc if d))
    return out


def parse_dfn(path: str) -> tuple[list[Field], set[str], bool]:
    """Fields of the data record, the data record-type tokens, and whether the
    record carries a leading RT marker field.

    DEFN lines look like

        DEFN 11 ST=RECD,RT=;Longitude:F12.6:UNIT=degrees,EPSGCODE=7844,NAME=...

    Several DEFNs accumulate into one record. `RT=` (blank) and `RT=DATA` are
    both the data record; `RT=COMM` is the comment block and is skipped.
    `;END DEFN` marks the end, except the 2017 file writes it on its
    second-to-last line as well, so it is stripped rather than trusted.
    """
    fields: list[Field] = []
    data_rts: set[str] = set()
    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line.upper().startswith("DEFN"):
                continue
            head, _, rest = line.partition(";")
            m = re.search(r"RT\s*=\s*([A-Za-z0-9_]*)", head)
            rt = (m.group(1) if m else "").upper()
            if rt in NON_DATA_RT:
                continue
            data_rts.add(rt)
            rest = re.sub(r";?\s*END\s+DEFN\s*$", "", rest, flags=re.I)
            for fd in rest.split(";"):
                if not fd.strip():
                    continue
                parts = fd.split(":")
                if len(parts) < 2:
                    continue
                name = parts[0].strip()
                fm = FMT.match(parts[1].strip())
                if not fm:
                    sys.exit(f"{path}: cannot read the format {parts[1]!r} "
                             f"of field {name!r}")
                a = parse_attrs(parts[2:])
                nul = a.get("NULL")
                fields.append(Field(
                    name=name, kind=fm.group(2).upper(), width=int(fm.group(3)),
                    count=int(fm.group(1) or 1),
                    decimals=int(fm.group(4)) if fm.group(4) else None,
                    unit=a.get("UNIT", a.get("UNITS", "")),
                    null=float(nul) if nul not in (None, "") else None,
                    epsg=a.get("EPSGCODE", ""),
                    desc=a.get("NAME", a.get("DESC", "")), rt=rt))
    if not fields:
        sys.exit(f"{path}: no data-record DEFN lines found")
    # An A-type field literally named RT is the record-type marker, not data.
    has_rt = fields[0].name.upper() == "RT" and fields[0].kind == "A"
    return fields, {r for r in data_rts if r}, has_rt


def record_length(fields: list[Field]) -> int:
    return sum(f.span for f in fields)


def is_data_line(ln: str, has_rt: bool, data_rts: set[str]) -> bool:
    """Whether a raw line is a data record.

    ONE predicate, used by both the width check and the reader. They disagreed
    in the first version: the check read line 1 of the file, which in a delivery
    that carries a COMM header block is a 62-character comment, and compared it
    against the 74-character data record. The result was a refusal to convert a
    perfectly good file. Neither real delivery has COMM lines in its .dat, so
    only a synthetic fixture found it.
    """
    if not ln.strip():
        return False
    head = ln[:4].strip().upper()
    if has_rt:
        return (not data_rts) or head in data_rts
    return head not in NON_DATA_RT


def first_data_line(dat: str, has_rt: bool, data_rts: set[str]) -> str:
    with open(dat, "r", errors="replace") as fh:
        for ln in fh:
            if is_data_line(ln, has_rt, data_rts):
                return ln.rstrip("\r\n")
    sys.exit(f"{dat}: no data records found")


def check_width(dat: str, fields: list[Field], has_rt: bool,
                data_rts: set[str]) -> int:
    """Assert the DFN layout matches the file, and return the record length.

    This is the check worth having. A one-character disagreement means every
    field after it is sliced off by one, and fixed-width numeric fields mostly
    still parse when shifted, so the output looks fine and is wrong.
    """
    want = record_length(fields)
    first = first_data_line(dat, has_rt, data_rts)
    if len(first) != want:
        sys.exit(
            f"{dat}: the DFN describes a {want}-character record but the first "
            f"line is {len(first)} characters.\n"
            f"  Fields: {sum(f.count for f in fields)} columns from "
            f"{len(fields)} DEFN entries.\n"
            f"  This is the mis-slicing trap: shifted fixed-width numbers still "
            f"parse, so a converter that carries on produces a plausible, wrong "
            f"file. Check whether the record carries a leading RT marker the DFN "
            f"does not declare.")
    return want


def decode_numeric(raw: np.ndarray, width: int, name: str,
                   bad: dict) -> np.ndarray:
    """A column of fixed-width bytes to float64, blanks and junk as NaN.

    ASEG writes an overflowed value as a field of asterisks, so this cannot
    assume every token parses; it counts what it could not read per column
    instead of failing the run or hiding it.
    """
    txt = np.char.strip(raw.astype("U"))
    out = np.full(txt.shape[0], np.nan, dtype=np.float64)
    good = txt != ""
    if not good.any():
        return out
    try:
        out[good] = txt[good].astype(np.float64)
    except ValueError:
        # rare: fall back to per-element so one bad token does not lose the
        # column, and report how many.
        n = 0
        idx = np.flatnonzero(good)
        for i in idx:
            try:
                out[i] = float(txt[i])
            except ValueError:
                n += 1
        if n:
            bad[name] = bad.get(name, 0) + n
    return out


def parse_block(lines: list[str], fields: list[Field], reclen: int,
                has_rt: bool, bad: dict) -> dict[str, np.ndarray]:
    """One block of records to column arrays, by slicing a 2-D byte view."""
    buf = bytearray()
    for ln in lines:
        s = ln.rstrip("\r\n").encode("ascii", "replace")
        buf += s.ljust(reclen)[:reclen]
    arr = np.frombuffer(bytes(buf), dtype="S1").reshape(len(lines), reclen)
    cols: dict[str, np.ndarray] = {}
    pos = 0
    for f in fields:
        for name in f.columns():
            sl = np.ascontiguousarray(arr[:, pos:pos + f.width])
            raw = sl.view(f"S{f.width}").ravel()
            pos += f.width
            if has_rt and f is fields[0]:
                continue                      # the RT marker is not data
            if f.kind == "A":
                cols[name] = np.char.strip(raw.astype("U"))
            else:
                v = decode_numeric(raw, f.width, name, bad)
                if f.null is not None:
                    v[np.isclose(v, f.null, rtol=1e-9, atol=0.0)] = np.nan
                cols[name] = v
    return cols


def to_arrow(v: np.ndarray, kind: str):
    """One column to Arrow, with missing values as NULLS rather than NaN.

    A NaN in a float column is not a null: `null_count` reads zero and a
    `WHERE x IS NULL` filter finds nothing, so the sentinel values this tool
    just finished converting would be invisible to the consumer. Everything
    missing -- a blank field, a declared NULL sentinel, an unreadable token --
    becomes a null here, and the DFN's sentinel is recorded in the column
    metadata as `aseg_null` so the original value is not lost.
    """
    import pyarrow as pa
    if kind == "A":
        return pa.array(v, type=pa.string())
    mask = ~np.isfinite(v)
    if kind == "I":
        iv = np.zeros(v.shape[0], dtype=np.int64)
        np.rint(v, out=v, where=~mask)
        iv[~mask] = v[~mask].astype(np.int64)
        return pa.array(iv, type=pa.int64(), mask=mask)
    return pa.array(v, type=pa.float64(), mask=mask)


def wkb_points(x: np.ndarray, y: np.ndarray):
    """A WKB point column, built as one buffer rather than per-row objects.

    Little-endian byte order flag, geometry type 1, then two float64s: 21 bytes
    a point. Rows with a missing coordinate become nulls, not points at NaN.
    """
    import pyarrow as pa
    n = x.shape[0]
    b = np.zeros((n, 21), dtype=np.uint8)
    b[:, 0] = 1                                        # little endian
    b[:, 1] = 1                                        # wkbPoint, uint32 LE
    b[:, 5:13] = np.ascontiguousarray(x, np.float64).view(np.uint8).reshape(n, 8)
    b[:, 13:21] = np.ascontiguousarray(y, np.float64).view(np.uint8).reshape(n, 8)
    valid = np.isfinite(x) & np.isfinite(y)
    offsets = np.arange(0, 21 * (n + 1), 21, dtype=np.int32)
    bufs = [pa.py_buffer(np.packbits(valid, bitorder="little")),
            pa.py_buffer(offsets), pa.py_buffer(b.tobytes())]
    nulls = int((~valid).sum())
    return pa.Array.from_buffers(pa.binary(), n, bufs, null_count=nulls), nulls


def crs_metadata(epsg: str, crs_arg: str | None, xname: str, yname: str,
                 fields: list[Field]):
    """PROJJSON for the geometry column, from the DFN or the command line."""
    code = crs_arg or (f"EPSG:{epsg}" if epsg else "")
    if not code:
        hints = [f"    {f.name}:{' ' * max(0, 12 - len(f.name))}{f.desc}"
                 for f in fields if f.name in (xname, yname)]
        sys.exit(
            "this DFN declares no EPSGCODE on the coordinate fields, so the CRS "
            "is not known.\n" + "\n".join(hints) + "\n"
            "  Pass it explicitly, e.g. --crs EPSG:4283 for GDA94 geographic, "
            "EPSG:7844 for GDA2020,\n"
            "  EPSG:28352 for GDA94 MGA zone 52, EPSG:7855 for GDA2020 MGA "
            "zone 55.\n"
            "  Guessing it from the description text is exactly the kind of "
            "silent error this refuses to make.")
    try:
        import pyproj
        return json.loads(pyproj.CRS.from_user_input(code).to_json()), code, True
    except ImportError:
        auth, _, num = code.partition(":")
        return {"id": {"authority": auth, "code": int(num)}}, code, False


# A fixture small enough to read, carrying every case that is easy to get
# wrong. The COMM record in it found a real defect: the width check used to
# read line 1 of the file and compare a 62-character comment against the
# 74-character data record, refusing to convert a good file. Neither real
# delivery has COMM lines in its .dat, so nothing else would have caught it.
SELFTEST_DFN = """DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76
DEFN 1 ST=RECD,RT=;Line:I10:NAME=Line Number
DEFN 2 ST=RECD,RT=;Longitude:F12.6:UNIT=degrees,EPSGCODE=4283,NAME=Longitude GDA94
DEFN 3 ST=RECD,RT=;Latitude:F12.6:UNIT=degrees,EPSGCODE=4283,NAME=Latitude GDA94
DEFN 4 ST=RECD,RT=;Depth:F10.2:NULL=-99999.99:UNIT=m,A depth, with a comma in its description
DEFN 5 ST=RECD,RT=;Gates:3F10.3:UNIT=pV,NAME=three gates;END DEFN
"""
SELFTEST_DAT = (
    "COMM a comment record, shorter than a data record, and skipped\n"
    "    100001  143.090911  -31.803289     12.34     1.100     2.200     3.300\n"
    "    100002  143.100000  -31.810000 -99999.99     1.000     2.000     3.000\n"
    "    100003                          5.00     0.500     0.600     0.700\n"
    "    100004  143.120000  -31.820000      7.00**********     2.500     2.600\n"
)


def selftest() -> int:
    """Convert the fixture and check every behaviour it exercises."""
    import tempfile
    import pyarrow.parquet as pq
    fails = []

    def ck(label, cond):
        print(f"  {'PASS' if cond else 'FAIL'}  {label}")
        if not cond:
            fails.append(label)

    with tempfile.TemporaryDirectory() as d:
        stem = os.path.join(d, "t")
        open(stem + ".dfn", "w").write(SELFTEST_DFN)
        open(stem + ".dat", "w").write(SELFTEST_DAT)
        fields, rts, has_rt = parse_dfn(stem + ".dfn")
        ck("array field 3F10.3 expands to three columns",
           [c for f in fields for c in f.columns()].count("Gates_1") == 1
           and sum(f.count for f in fields) == 7)
        ck("no RT marker declared on the data record", has_rt is False)
        ck("record length is 74 from the DFN", record_length(fields) == 74)
        ck("the comment record is not mistaken for the data record",
           len(first_data_line(stem + ".dat", has_rt, rts)) == 74)
        d4 = [f for f in fields if f.name == "Depth"][0]
        ck("NULL sentinel read from its own colon section", d4.null == -99999.99)
        ck("a description containing a comma survives intact",
           d4.desc == "A depth, with a comma in its description")
        ck("EPSGCODE read from the coordinate field",
           [f for f in fields if f.name == "Longitude"][0].epsg == "4283")

        out = os.path.join(d, "t.parquet")
        sys.argv = [sys.argv[0], stem, "-o", out]
        main()
        t = pq.read_table(out)
        ck("the comment record is skipped, leaving four rows", t.num_rows == 4)
        ck("I-format stays an integer",
           t.schema.field("Line").type == "int64"
           and t.column("Line")[0].as_py() == 100001)
        ck("the NULL sentinel becomes a null",
           t.column("Depth").to_pylist() == [12.34, None, 5.0, 7.0])
        ck("a blank field becomes a null",
           t.column("Longitude")[2].as_py() is None)
        ck("an unreadable token becomes a null",
           t.column("Gates_1")[3].as_py() is None)
        ck("a row with no coordinate gets a null geometry",
           t.column("geometry").null_count == 1)
        g = t.column("geometry")[0].as_py()
        import struct
        ck("the WKB point is 21 bytes, little-endian, type 1",
           len(g) == 21 and g[0] == 1 and struct.unpack("<I", g[1:5])[0] == 1)
        ck("the WKB coordinates are the file's coordinates",
           struct.unpack("<dd", g[5:21]) == (143.090911, -31.803289))
        geo = json.loads(t.schema.metadata[b"geo"])
        ck("GeoParquet metadata names the WKB point column",
           geo["primary_column"] == "geometry"
           and geo["columns"]["geometry"]["encoding"] == "WKB")
        ck("the unit and description reach the column metadata",
           (t.schema.field("Depth").metadata or {}).get(b"unit") == b"m")

    print(f"\n{'all checks passed' if not fails else str(len(fails)) + ' FAILED'}")
    return 1 if fails else 0


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert an ASEG-GDF2 AEM delivery to GeoParquet.")
    ap.add_argument("stem", help="path stem, or the .dat/.dfn itself")
    ap.add_argument("-o", "--out", help="output .parquet (default: <stem>.parquet)")
    ap.add_argument("--xy", nargs=2, metavar=("X", "Y"),
                    help="coordinate fields (default: Longitude Latitude, "
                         "else Easting Northing)")
    ap.add_argument("--crs", help="override the CRS, e.g. EPSG:4283")
    ap.add_argument("--limit", type=int, help="stop after N records")
    ap.add_argument("--chunk", type=int, default=100_000,
                    help="records per row group (default 100000). This sets "
                         "peak memory: 100000 records of the 196-column UDF "
                         "delivery peaks near 1.5 GB, and the cost scales with "
                         "it, so lower it on a small machine")
    ap.add_argument("--compression", default="zstd")
    ap.add_argument("--schema", action="store_true",
                    help="print the parsed layout and exit, writing nothing")
    ap.add_argument("--selftest", action="store_true",
                    help="convert a built-in fixture and check the result; "
                         "needs no delivery files")
    ap.add_argument("--bbox", action="store_true",
                    help="record the geometry bbox, which needs a first pass "
                         "over the coordinates")
    a = ap.parse_args() if "--selftest" not in sys.argv else None
    if a is None:
        sys.exit(selftest())

    stem = re.sub(r"\.(dat|dfn|des|prj)$", "", a.stem, flags=re.I)
    dfn, dat = stem + ".dfn", stem + ".dat"
    for p in (dfn, dat):
        if not os.path.isfile(p):
            sys.exit(f"not found: {p}")

    fields, data_rts, has_rt = parse_dfn(dfn)
    ncol = sum(f.count for f in fields) - (1 if has_rt else 0)

    if a.schema:
        print(f"{dfn}\n  {len(fields)} DEFN entries -> {ncol} columns"
              f"{'  (leading RT marker, not a column)' if has_rt else ''}")
        print(f"  record types: {sorted(data_rts) or ['(blank)']}")
        print(f"  record length: {record_length(fields)} characters")
        print(f"\n  {'field':22} {'fmt':10} {'cols':>5} {'unit':14} null")
        for f in fields:
            fmt = (f"{f.count if f.count > 1 else ''}{f.kind}{f.width}"
                   + (f".{f.decimals}" if f.decimals is not None else ""))
            print(f"  {f.name:22} {fmt:10} {f.count:5d} {f.unit:14} "
                  f"{'' if f.null is None else f.null}"
                  + (f"   [EPSG:{f.epsg}]" if f.epsg else ""))
        return

    # Everything above here runs on numpy alone, which is why --schema does.
    # Only the write needs pyarrow, and a bare ModuleNotFoundError traceback is
    # not an answer to "why did this not work".
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        sys.exit(
            "writing Parquet needs pyarrow, which this python does not have:\n"
            f"    {sys.executable}\n"
            "  Install it:  python3 -m pip install pyarrow\n"
            "  The project's conda env has no pyarrow either, as of Nov 2026.\n"
            "  --schema works without it, if the layout is all you need.")

    names = [c for f in fields for c in f.columns()]
    if has_rt:
        names = names[1:]
    if a.xy:
        xname, yname = a.xy
    elif "Longitude" in names and "Latitude" in names:
        xname, yname = "Longitude", "Latitude"
    elif "Easting" in names and "Northing" in names:
        xname, yname = "Easting", "Northing"
    else:
        sys.exit("no Longitude/Latitude or Easting/Northing field; pass --xy")
    for n in (xname, yname):
        if n not in names:
            sys.exit(f"--xy field {n!r} is not in this file. Columns: "
                     + ", ".join(names[:12]) + " ...")

    byname = {f.name: f for f in fields}
    epsg = byname.get(xname, Field("", "F", 0)).epsg or \
        byname.get(yname, Field("", "F", 0)).epsg
    crs_json, crs_code, exact = crs_metadata(epsg, a.crs, xname, yname, fields)

    out = a.out or (stem + ".parquet")
    reclen = check_width(dat, fields, has_rt, data_rts)

    bbox = None
    if a.bbox:
        bbox = bbox_pass(dat, fields, reclen, has_rt, xname, yname,
                         data_rts, a.limit)

    # Column metadata from the DFN, so the parquet is self-describing.
    meta_for = {}
    for f in fields:
        for c in f.columns():
            m = {}
            if f.unit:
                m[b"unit"] = f.unit.encode()
            if f.desc:
                m[b"description"] = f.desc.encode()
            if f.null is not None:
                m[b"aseg_null"] = repr(f.null).encode()
            meta_for[c] = m or None

    # I-format fields stay integers. Emitting them as float64 made
    # `Date_Local` read 20220329.0 and `Line` 100001.0, which is exact but
    # wrong to look at and wrong to join on.
    kind_of: dict[str, str] = {}
    arrow_fields = []
    for f in fields:
        for c in f.columns():
            if has_rt and c == fields[0].name:
                continue
            kind_of[c] = f.kind
            t = (pa.string() if f.kind == "A"
                 else pa.int64() if f.kind == "I" else pa.float64())
            arrow_fields.append(pa.field(c, t, metadata=meta_for.get(c)))
    geo = {"version": "1.1.0", "primary_column": "geometry",
           "columns": {"geometry": {"encoding": "WKB",
                                    "geometry_types": ["Point"],
                                    "crs": crs_json}}}
    if bbox:
        geo["columns"]["geometry"]["bbox"] = list(bbox)
    arrow_fields.append(pa.field("geometry", pa.binary()))
    schema = pa.schema(arrow_fields, metadata={
        b"geo": json.dumps(geo).encode(),
        b"aseg_dfn": os.path.basename(dfn).encode(),
        b"aseg_source": os.path.basename(dat).encode()})

    bad: dict[str, int] = {}
    nrows = geom_nulls = 0
    writer = pq.ParquetWriter(out, schema, compression=a.compression)
    try:
        for block in blocks(dat, reclen, has_rt, data_rts, a.chunk, a.limit):
            cols = parse_block(block, fields, reclen, has_rt, bad)
            g, nn = wkb_points(cols[xname], cols[yname])
            geom_nulls += nn
            arrays = []
            for fl in schema:
                if fl.name == "geometry":
                    arrays.append(g)
                else:
                    arrays.append(to_arrow(cols[fl.name], kind_of[fl.name]))
            writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
            nrows += len(block)
            print(f"\r  {nrows:,} records", end="", file=sys.stderr, flush=True)
    finally:
        writer.close()
    print(file=sys.stderr)

    size = os.path.getsize(out)
    print(f"wrote {out}")
    print(f"  {nrows:,} records x {ncol + 1} columns "
          f"({size / 1e6:.1f} MB, {a.compression})")
    print(f"  geometry: WKB point from {xname}/{yname}, {crs_code}"
          f"{'' if exact else '  (EPSG id only: pyproj absent)'}")
    if geom_nulls:
        print(f"  {geom_nulls:,} records had no coordinate and carry a null "
              f"geometry")
    if bbox:
        print(f"  bbox: {bbox[0]:.6f} {bbox[1]:.6f} {bbox[2]:.6f} {bbox[3]:.6f}")
    for k, v in bad.items():
        print(f"  WARNING {v:,} unreadable values in {k}, written as null")


def blocks(dat: str, reclen: int, has_rt: bool, data_rts: set[str],
           chunk: int, limit: int | None):
    """Yield lists of data lines, skipping comment records."""
    n = 0
    buf: list[str] = []
    with open(dat, "r", errors="replace") as fh:
        for ln in fh:
            if not is_data_line(ln, has_rt, data_rts):
                continue
            buf.append(ln)
            n += 1
            if len(buf) >= chunk:
                yield buf
                buf = []
            if limit and n >= limit:
                break
    if buf:
        yield buf


def bbox_pass(dat, fields, reclen, has_rt, xname, yname, data_rts, limit):
    """Minimum and maximum of the two coordinate columns, read on their own."""
    off = {}
    pos = 0
    for f in fields:
        for c in f.columns():
            off[c] = (pos, f.width)
            pos += f.width
    (xs, xw), (ys, yw) = off[xname], off[yname]
    lo = np.array([np.inf, np.inf])
    hi = np.array([-np.inf, -np.inf])
    bad: dict = {}
    for block in blocks(dat, reclen, has_rt, data_rts, 200_000, limit):
        for (s, w), k in (((xs, xw), 0), ((ys, yw), 1)):
            raw = np.array([ln[s:s + w] for ln in block], dtype=f"S{w}")
            v = decode_numeric(raw, w, "bbox", bad)
            v = v[np.isfinite(v)]
            if v.size:
                lo[k] = min(lo[k], v.min())
                hi[k] = max(hi[k], v.max())
    if not np.isfinite(lo).all():
        return None
    return (lo[0], lo[1], hi[0], hi[1])


if __name__ == "__main__":
    main()
