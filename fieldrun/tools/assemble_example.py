#!/usr/bin/env python3
"""assemble_example.py -- build fieldrun/example/ from the GA source archives.

This runs on a machine that holds the Upper Darling Floodplain deliveries. It
is here for provenance, not for a reader: the sources it reads are not
published, and the two-site extract it writes is what ships instead.

What it does NOT do: subset the 4.9 GB ASEG-GDF EM delivery. The per-site
sounding CSVs were produced by udf_realdata/export_udf_soundings.py, which
seek-reads chosen records by number; this script copies its output. The nav
GeoParquet IS clipped here, from 1,764,231 rows to the points within RADIUS_M
of the two sites, so a reader can see how a site was matched to a sounding.

    python3 fieldrun/tools/assemble_example.py [--sites ID ID] [--radius 1000]
"""
from __future__ import annotations
import argparse, csv, hashlib, os, shutil, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "example"))
UDF = "/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/udf_realdata"
SNMR_ARCHIVE = "/home/symo/Documents/GA/data/SNMR_dat/Field_Inversions"
AEM_SYSTEM = "/home/symo/Documents/GA/data/AEM/UDF_AEM/01_EM/UDF_AEM_Data"

# Vista Clara writes several products per site; the operator's preferred one is
# the most processed. Same order as udf_realdata/udf_smr.jl udf_matfile().
MAT_PREFERENCE = ("harmonic_removed_noemi", "noemi", "harmonic_removed", "")

# These carry absolute paths into the author's archive and are dropped.
META_DROP = ("pathname", "gmr_inv_path")


def pick_mat(field_id):
    d = os.path.join(SNMR_ARCHIVE, field_id, "Processed")
    if not os.path.isdir(d):
        sys.exit(f"no Processed directory for {field_id}: {d}")
    files = sorted(f for f in os.listdir(d) if f.endswith(".mat"))
    for pat in MAT_PREFERENCE:
        for f in files:
            if pat in f:
                return os.path.join(d, f), f
    sys.exit(f"no .mat for {field_id}")


def subset_csv(src, dst, id_col, ids, drop=()):
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh))
    keep = [r for r in rows if r[id_col] in ids]
    missing = set(ids) - {r[id_col] for r in keep}
    if missing:
        sys.exit(f"{src}: no row for {sorted(missing)}")
    cols = [c for c in rows[0] if c not in drop]
    # lineterminator="\n" is what controls it. The previous attempt set
    # newline="\n" and claimed to have fixed the CRLF endings; a reviewer ran
    # this function and got output BYTE-IDENTICAL to the CRLF file, because
    # csv.DictWriter writes its own lineterminator regardless of newline=.
    # A comment claiming a fix that changed nothing is worse than no comment.
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in sorted(keep, key=lambda r: r[id_col]):
            w.writerow(r)
    return len(keep), cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", nargs="+", default=["SMR_UD21", "SMR_UD24"])
    ap.add_argument("--radius", type=float, default=1000.0)
    # The three source roots were hardcoded to the author's machine with no way
    # to override them, so a reader holding the deliveries had to edit the file.
    ap.add_argument("--udf", default=UDF, help="the udf_realdata working directory")
    ap.add_argument("--snmr-archive", default=SNMR_ARCHIVE)
    ap.add_argument("--aem-system", default=AEM_SYSTEM, help="path stem, no extension")
    a = ap.parse_args()
    ids = set(a.sites)
    global UDF, SNMR_ARCHIVE, AEM_SYSTEM
    UDF, SNMR_ARCHIVE, AEM_SYSTEM = a.udf, a.snmr_archive, a.aem_system

    for d in ("snmr", "aem/soundings", "aem/system"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)

    written = []

    # 1. site table and SNMR metadata, two rows each
    n, _ = subset_csv(os.path.join(UDF, "sites_udf.csv"),
                      os.path.join(OUT, "sites.csv"), "Field_ID", ids)
    written.append(("sites.csv", f"{n} rows of sites_udf.csv"))

    n, cols = subset_csv(os.path.join(UDF, "udf_meta.csv"),
                         os.path.join(OUT, "snmr_metadata.csv"), "Field_ID", ids,
                         drop=META_DROP)
    written.append(("snmr_metadata.csv", f"{n} rows, {len(cols)} cols, dropped {META_DROP}"))

    # 2. the processed SNMR sounding each site is inverted from
    for sid in sorted(ids):
        src, name = pick_mat(sid)
        shutil.copy2(src, os.path.join(OUT, "snmr", name))
        written.append((f"snmr/{name}", f"{os.path.getsize(src)} B, chosen by preference order"))

    # 3. the extracted AEM soundings and their manifest
    # The manifest's `csv` column holds the AUTHOR's source-tree path,
    # aem_soundings/<ID>.csv, while the shipped files are at aem/soundings/.
    # Copying the rows unchanged left a path that resolves to nothing, which a
    # reviewer found. Rewrite it to where the file actually is.
    man_dst = os.path.join(OUT, "aem", "soundings", "manifest.csv")
    n, _ = subset_csv(os.path.join(UDF, "aem_soundings", "manifest.csv"),
                      man_dst, "field_id", ids)
    with open(man_dst, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        if r.get("csv"):
            r["csv"] = "aem/soundings/" + os.path.basename(r["csv"])
    with open(man_dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    written.append(("aem/soundings/manifest.csv",
                    f"{n} rows, csv column repointed at aem/soundings/"))
    for sid in sorted(ids):
        src = os.path.join(UDF, "aem_soundings", f"{sid}.csv")
        dst = os.path.join(OUT, "aem", "soundings", f"{sid}.csv")
        shutil.copy2(src, dst)
        with open(dst) as fh:
            ngate = sum(1 for _ in fh) - 1
        written.append((f"aem/soundings/{sid}.csv", f"{ngate} gates"))

    # 4. the nav GeoParquet, clipped
    sites = pd.read_csv(os.path.join(UDF, "sites_udf.csv")).set_index("Field_ID")
    nav = pd.read_parquet(os.path.join(UDF, "udf_nav.parquet"))
    keep = np.zeros(len(nav), bool)
    for sid in sorted(ids):
        d = np.hypot(nav.Easting - sites.loc[sid, "easting"],
                     nav.Northing - sites.loc[sid, "northing"])
        keep |= d < a.radius
    sub = nav[keep].reset_index(drop=True)
    p = os.path.join(OUT, "aem", "nav_clip.parquet")
    sub.to_parquet(p, index=False)
    # Plain Parquet, not GeoParquet: the file carries pandas metadata and no
    # `geo` key, which a reviewer checked. It is also not byte-reproducible --
    # pyarrow and pandas versions are embedded -- so its checksum is a record of
    # what shipped, not a target to reproduce.
    written.append(("aem/nav_clip.parquet",
                    f"{len(sub)} of {len(nav)} rows, {sub.Line.nunique()} lines, "
                    f"within {a.radius:.0f} m"))

    # 5. the SkyTEM system description -- gates, ramps, geometry, low-pass
    for ext in (".des", ".dfn", ".prj"):
        src = AEM_SYSTEM + ext
        shutil.copy2(src, os.path.join(OUT, "aem", "system", os.path.basename(src)))
        written.append((f"aem/system/{os.path.basename(src)}", f"{os.path.getsize(src)} B"))

    # 6. checksums over everything shipped
    # PROVENANCE.md is written by hand, not by this script, and it lives in the
    # same directory -- so a naive walk hashed it and a re-run produced 12 lines
    # against the 11 shipped. A reviewer tried to reproduce the checksums and
    # could not. Only the files this script writes are hashed.
    SKIP = {"SHA256SUMS", "PROVENANCE.md"}
    lines, total = [], 0
    for root, _, files in os.walk(OUT):
        for f in sorted(files):
            if f in SKIP:
                continue
            fp = os.path.join(root, f)
            rel = os.path.relpath(fp, OUT)
            h = hashlib.sha256(open(fp, "rb").read()).hexdigest()
            lines.append(f"{h}  {rel}")
            total += os.path.getsize(fp)
    with open(os.path.join(OUT, "SHA256SUMS"), "w") as fh:
        fh.write("\n".join(sorted(lines, key=lambda s: s.split("  ", 1)[1])) + "\n")

    for name, note in written:
        print(f"  {name:38s} {note}")
    print(f"\n  {len(lines)} files, {total/1e6:.2f} MB total, SHA256SUMS written")


if __name__ == "__main__":
    main()
