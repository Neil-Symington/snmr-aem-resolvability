# The synthetic columns, and which figure needs which

Split out of `README.md` on 12 Nov 2026. The README follows one worked example,
`2P`, because that is all a reader needs to run the thing. Everything here is
per-column detail that only matters once you are rebuilding the paper's figures,
and most of it was measured the hard way.

## The columns

| tag | what it is |
|---|---|
| `3P` | test one, the unconfined Loxton-Parilla sands |
| `2P` | test two, the confined aquifer pair, and the README's worked example |
| `3S`, `2S` | saline controls of `3P` and `2P` |
| `3M`, `2M` | the paper's Stage 1 tags for `3P` and `2P`, see below |

A saline control shares its parent's clay, porosity and water content bit for
bit and differs *only* in pore-water conductivity, so any difference between a
column and its control can be nothing but the water. A control names its parent
in `[data].parent` and the validator refuses one that does not: comparing a
control against the wrong parent has already produced one retracted finding in
this project.

`synthetics/config/paper_3P.toml` is the paper's test-one column and the only one
of the six that ships. With it and `confined_aquifer.toml`, which is `2P`, you
have the two columns `fig04`, `fig05` and `fig06` need.

Note that test one is the UNCONFINED column. Figures and notes written before
13 Oct 2026 say the opposite.

## Which script reads which column

Written down nowhere until a reviewer went looking, and then wrong for three of
its rows in its first version: `fig02` and `figS8` read four columns, not one or
two, and `fig06` reads two. Counts below are from the tags each script actually
references.

| script | columns | also needs |
|---|---|---|
| `fig02` | 3P, 3S, 2P, 2S | `data/truth_murray<tag>.jld2` |
| `fig03` | 3P, 2P (+ 2PJ if present) | `tools/extract_full_draws.jl`, `tools/aem_sensitivity.jl` outputs |
| `fig04`, `fig05` | 3P, 2P | truth files and Stage 3 outputs |
| `fig06` | 3P, 2P (+ 2PJ if present) | Stage 5 outputs |
| `figS8` | 3P, 3S, 2P, 2S | Stage 4 outputs |
| `fig07`, `export_numbers_v2` | 3P, 3S, 2P, 2S | Stage 4 post, and `singlepass_archive/` if present |
| `fig08` | 2P | `tools/section6_bore_value.jl` |

`singlepass_archive/` is read by `fig07` and **written by nothing**. It is an
archive of an earlier grid mode and its absence is handled.

## Building the 3P column, and the two things that will stop you

Both measured by a reviewer on a clean clone.

**1. The smoke gate grid does not bracket 3P.** `smoke_2P.toml`'s
`n_nodes_coarse = 4, n_nodes_fine = 3` fails on a `3P` column after 19 minutes
with *"the coarse pass found 1 reachable node(s) of 4. The fine grid cannot be
bracketed."* `3P`'s reachable range is roughly 0.0037 to 0.0215 S/m, 0.57 dex,
and a 4-node grid 1.33 dex apart lands one node in it at best. Use:

```toml
[gate]
n_nodes_coarse = 8      # 4 does NOT bracket 3P
n_nodes_fine   = 4
```

That brackets, and takes 53 minutes to Stage 4 post on smoke ensembles. The grid
is column-specific, which is why the smoke config is named for `2P`.

**2. Stage 5 on 3P needs `data/aem_draws_3M.jld2`, which nothing builds.**
`configs/synth_3P.toml` points at the `3M` tag for the same reason
`synth_2P.toml` points at `2M`: the paper's Stage 1 for that column ran under the
older tag, and σ_bulk is identical. The chain builds `aem_draws_3P`, so:

```bash
cp data/aem_draws_3P.jld2 data/aem_draws_3M.jld2
```

`run_site.jl` prints `STAGE 1 INPUT MISSING: data/aem_draws_3M.jld2` early, but
it is non-fatal and the chain hint never mentions `3M`, so the run proceeds
through Stages 3 and 4 and then dies in Stage 5 an hour later. A reviewer lost
that hour.

## Why the label check exists

All eight figures carrying column tags call `FigStyle.assert_column`. When that
was first claimed four were unwired, including `figS8`, which was one of three a
reviewer drew from fabricated columns with no message at all.

The check compares the tag against the `column` the Stage 0 artefact recorded,
and the artefact's `truth_file_sha` against the hash of
`data/truth_murray<tag>.jld2`. The second is needed because `column` is copied
from the config: a config declaring `column = "3P"` while pointing at another
column's truth is sealed cleanly by the real Stage 0, and a reviewer did exactly
that and drew two figures from it.
