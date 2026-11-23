# Figure inputs and outputs

`workflow/figures.py` writes numbered PNG, SVG and PDF files in the selected output folder. It also writes SHA-256 digests and keeps renderer intermediates under `runs/`. Each group needs a new output run directory. Renderings are made from numerical inputs. Page composition, typefaces, map baselayers, and manual vector edits can differ from the submitted artwork.

| Figure numbers | Group | Inputs for the redraw |
| --- | --- | --- |
| 1, 5, S1–S5, S11, S14, S16 | `workbooks` | S1, S3 and S4 workbook sheets read by `plots/workbook_panels.py`. Figure 1 redraws station coordinates and study design. Figure 5 refits the line through the ten lake means. |
| 2 | `community` | Seven source tables derived from S1 and S2 by `paper_inputs.py`. Includes the paired model, the 21-library composition display and the occurrence tables. |
| 3 | `carriers` | Five final panel tables from S5. Also regenerable with `carrier_attribution.py` and `figure3_tables.py` using the local source matrices. Panel C includes all 12 displayed module columns. |
| 4, S9 | `traits` | S4, the midpoint tree in S5 and equal-MAG-weight pI grids. Supply archived grids or generate them from per-protein pI with `workflow/protein_density.R`. The submitted Figure 4 includes manual Illustrator work. |
| S6, S7, S10, S12, S13, S15, S19 | `traits` (also included in `paper`) | S4 traits, marker tests, lake/DOM diagnostics and the S3 Figure 3 addendum for ion and carrier tables. |
| S8, S17, S18 | `display` | S5 `inputs/si_display_sources`, including the independently archived S8 density table. |
| S20 | `coverage` | S5 Nonpareil sample estimates, observed points and fitted curves. |
| S21, S22, S23 | `assembly` | S5 catalogue recovery, paired response estimates, dispersion and directional sensitivity results. |
| S24 | `tree` | S5 protein/nucleotide tree comparisons and within-lineage model outputs. |
| S25 | `covariates` | S5 sequencing-coverage covariate diagnostics and omission models. |
| S26 | `matches` | S5 reciprocal-ANI population matches, abundances and concordance summaries. |
| S27, S28 | `fcs_quality`, `fcs_inference` | S6 review-stage quality and fixed-reference inference tables. Both R plotting scripts are in `plots/`. For S28, optional arguments select freshly refitted ecological and phylogenetic models. |

The optional Figure 1a cartography program uses the archived regional-map drawing code and joins the station coordinates from S1. Its local layer directory contains the four compressed SRTM tiles, six regional relief tiles, lake and stream geometries, grassland masks, the provincial boundary layer and the map settings. The thematic locator used the GS(2023)2766 standard map as its boundary reference. The watercourse arrows use the direction of decreasing elevation. The public script reads these layers from disk. The workbook-only Figure 1 route requires no GIS software or map layers.

```bash
python -m pip install 'shapely>=2' pyproj Pillow
python plots/site_map.py --s1 ../books/Supplementary_Data_S1.xlsx --layers "$LOCAL/map_layers" --out ../cartography
```

Map inputs originate from NASA SRTMGL1 and Terrarium relief tiles, ESA WorldCover 2021 class 30, OpenStreetMap lake relations 11763855 and 2309319 and associated waterways, and the DataV provincial geometry. The local input pack retains the layers used for the redraw. [Elevation tiles](https://registry.opendata.aws/terrain-tiles/) · [ESA WorldCover](https://esa-worldcover.org/en) · [OpenStreetMap](https://www.openstreetmap.org/copyright)

With S1–S5 and the S3 addendum, `--group paper` draws 29 figure numbers. Adding S6 enables S27 and S28, bringing the total to 31. To render all 33, supply both grids below in the prepared input directory and use `--group all`. Each is a TSV with columns `group`, `x`, `y`, `lo` and `hi`. The curves weight genomes equally within a lineage group. Replacing them with a density of MAG medians would describe a different distribution.

```text
from-supplements/
  extra/all_MAG_density.tsv
  extra/HQ_density.tsv
```

```bash
Rscript --vanilla workflow/protein_density.R "$W/traits/protein_pI.tsv" \
  ../books/Supplementary_Data_S4.xlsx ../from-supplements/extra
python workflow/figures.py --group all --books ../books --project ../zabuye-replay \
  --local-inputs ../from-supplements --out ../all-figures
```

The R command builds the grids from the `mag_id` and `pI` columns of the primary-traits output. It also accepts the two-column gzipped checkpoint. Use S4 for the archived MAG identities or the corresponding catalogue TSV for fresh identities. Each Zabuye water MAG has its own density on a common 512-point grid with one bandwidth estimated from all included proteins. The interval resamples MAGs 1,999 times within each group. The same bandwidth and the archived pre-R-4.4 FFT grid convention apply to the full, HQ and protein-pooled curves. The third output, `pooled_protein_density.tsv`, can be compared with the S8 source in S5.

For a visual comparison against the submitted artwork, supply the clean, author-held Word files to the separate exporter. It matches each figure caption to its embedded image, checks all expected figure numbers and records image digests. This operation exports existing artwork and performs no analysis.

```bash
python workflow/artwork.py --main "$MAIN_DOCX" --si "$SI_DOCX" --out ../embedded-artwork
```
