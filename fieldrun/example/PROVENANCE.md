# Provenance of the two-site example

These eleven files, 1.59 MB, are a two-site extract from Geoscience Australia's
Upper Darling Floodplain surveys. They are inputs, not results. Nothing here was
produced by this repository.

**Release status.** Publication of these files is a Geoscience Australia
decision. If you are reading this in a public clone, that decision has been
made; if the directory is absent, it has not, and you should ask the authors.

## The two sites

| | SMR_UD21 | SMR_UD24 |
|---|---|---|
| Area | Louth | Salt interception scheme |
| MGA zone 55 easting, northing | 315059, 6622677 | 371232, 6654668 |
| Nearest AEM sounding | 67.4 m, line 238201, record 771528 | 73.4 m, line 320101, record 1474272 |
| Control bore | GW036942.2.2, 282 m away, 50.25 m deep | GW098204.1.1, 499 m away, 65.69 m deep |
| Water table (bore) | 10.25 m | 8.69 m |
| SNMR loop | 100 m square, FID, 70 ms pulse | 100 m square, FID, 60 ms pulse |
| Inclination, declination | -61.92°, 9.42° | -61.57°, 9.54° |
| GA Record 2024/29 | interbedded silt/clay, low mobile water | salt interception scheme area, Figure 5.8 |

The pair is deliberate. UD21 is a site where the honest answer is likely *not
resolved*; UD24 is the site that anchors the paired-conductance comparison. One
of each verdict is worth more as an example than two passes, and both have bore
control a reader can check against.

## What each file is, and where it came from

| Shipped file | Source | How |
|---|---|---|
| `sites.csv` | `udf_realdata/sites_udf.csv` (48 rows) | two rows, verbatim |
| `snmr_metadata.csv` | `udf_realdata/udf_meta.csv` | two rows; `pathname` and `gmr_inv_path` dropped because they are absolute paths into the author's archive |
| `snmr/SMR_UD21_noemi.mat` | `.../Field_Inversions/SMR_UD21/Processed/` | Vista Clara GMR product, copied byte for byte |
| `snmr/SMR_UD24_harmonic_removed_noemi.mat` | `.../Field_Inversions/SMR_UD24/Processed/` | as above |
| `aem/soundings/*.csv` | `udf_realdata/aem_soundings/` | 41 gates each: LM 9–26, HM 16–38, every gate the delivery carries |
| `aem/soundings/manifest.csv` | same | two rows |
| `aem/nav_clip.parquet` | `udf_realdata/udf_nav.parquet` | clipped to the 1,336 of 1,764,231 points within 1 km of either site — 13 flight lines |
| `aem/system/UDF_AEM_Data.des/.dfn/.prj` | `data/AEM/UDF_AEM/01_EM/` | verbatim; the SkyTEM system description from which gates, ramps, geometry and low-pass are parsed |

`tools/assemble_example.py` writes all of it, and `SHA256SUMS` covers every
file. Re-running it on a machine with the deliveries should reproduce the
checksums exactly.

## What is deliberately not here

**The EM delivery.** `UDF_AEM_Data.dat` is 4.9 GB. The per-site sounding CSVs
were extracted from it by `udf_realdata/export_udf_soundings.py`, which
seek-reads the chosen records by number rather than materialising the file. That
script needs the delivery and is not part of the batch run. The nav clip is here
so a reader can see *how* a site was matched to a sounding — which line, which
record, how far away — without it.

**The SNMR archive.** 727 MB across 30 sites. Two preferred products ship.

**Bore logs.** `udf_realdata/bore/udf_borehole_wc.csv` is control data used to
check results, not an input to the inversion, and is left for the report rather
than the run.

## Noise model

AEM, from AUS_10087_report.pdf p.33: *"a nominal uncertainty of 3% was applied
to each datum, as well as the calculated relative uncertainty as provided with
the EM data"*, so `sigma_i = |d_i| * sqrt(0.03^2 + RUNC_i^2)`. The `runc` column
in each sounding CSV is that per-datum term; `sigma_pV` is the combination, so a
reader can check the arithmetic rather than trust it.

SNMR noise is not shipped as a number. It is estimated per pulse moment from the
FID residuals at run time (`sigma_hat_branch = "fid_residual_per_q"`), and a
hierarchical scale factor `s` is inferred alongside the model. The records are
effectively noiseless as delivered, which is why that scale factor is not
optional.
