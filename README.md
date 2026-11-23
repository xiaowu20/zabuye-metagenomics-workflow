# Zabuye metagenomics workflow

Code for *Multilevel metagenomics distinguishes community turnover, genome-resolved functional redistribution and lineage-associated proteome acidity in a hypersaline lake on the Tibetan Plateau*. The analysis starts with 35 physical libraries from Zabuye Salt Lake and Dangxiong Co. Four grouped assemblies yielded the study catalogue. The primary within-Zabuye comparisons use seven stations sampled on both 2 April and 2 October 2023. Summer and winter basin composites have a descriptive role.

The repository contains study code, software specifications and run instructions under the MIT licence. Study data, databases, plot inputs and figures are distributed separately. Use a sibling directory for inputs, outputs and compute workspaces. [中文使用说明](docs/quickstart-zh.md)

## Start from the paper supplements

For the revision analyses alone, start with `Supplementary_Data_S5.zip`. Add S6 for the contamination-screening sensitivity figures. The main figures also use workbooks S1–S4 and the Figure 3 addendum. Download the authors’ supplements in their original formats. Staging checks the SHA-256 digest of each supplied S5/S6 ZIP and creates a new workspace.

```bash
git clone --branch v1.1.0 https://github.com/xiaowu20/zabuye-metagenomics-workflow.git
cd zabuye-metagenomics-workflow
conda env create -f env/replay.yml
conda activate zabuye-replay
python -m unittest discover -s tests -v

# Set BOOKS to the directory containing the study supplements.
BOOKS=../books
P=../zabuye-replay
python workflow/stage.py --mode replay --project "$P" \
  --s5 "$BOOKS/Supplementary_Data_S5.zip" --s6 "$BOOKS/Supplementary_Data_S6_FCS.zip"
python workflow/paper_inputs.py --books "$BOOKS" --project "$P" --out ../from-supplements
python workflow/figures.py --books "$BOOKS" --local-inputs ../from-supplements \
  --project "$P" --group paper --out ../figures-from-paper
```

With S1–S6 and the S3 addendum this route draws 31 figure numbers and writes a PDF, SVG and PNG for each. Figure 4 and Figure S9 use two additional equal-MAG-weight pI density grids. Supply the author-held grids or regenerate them from protein-level pI and S4 with `workflow/protein_density.R` to draw all 33. The [figure index](docs/figures.md) gives the exact command.

For a quick S5-only example, the command below draws the 14-library Nonpareil diagnostic from its archived curve tables. The group runner accepts paths for the books and plotting inputs even when a selected group reads only S5.

```bash
python workflow/stage.py --s5 "$BOOKS/Supplementary_Data_S5.zip" --project ../s5-only
python workflow/figures.py --group coverage --project ../s5-only --books "$BOOKS" \
  --local-inputs ../from-supplements --out ../coverage-figure
```

The S5 replay refits the seven catalogue specifications, sequencing-coverage sensitivity and phylogenetic models inside the **newly staged workspace**. Figures S22 and S23 read the catalogue-specific fit tables directly. To include the external lake models, supply the S4 workbook directory and an output directory. Each command below defines a complete run. The default prints the command plan. Add `--execute` to perform the calculations.

```bash
# S5 alone
python workflow/run.py replay --project ../s5-only --execute
# S5 and the external S4 fits
python workflow/run.py replay --project "$P" --books "$BOOKS" --out ../model-fits --execute
```

## Start from the sequencing runs

The NCBI RunInfo metadata provide the physical library names and SRA accessions. The manifest builder checks mate filenames and input order against the four MEGAHIT option files in S5, reads the historical input count for each library from the corresponding assembly log, and verifies the paired subset. [NCBI BioProject PRJNA1531265](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1531265)

```bash
curl -fL 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1531265' -o ../RunInfo.csv
python workflow/manifest.py --s5 "$BOOKS/Supplementary_Data_S5.zip" --runs ../RunInfo.csv \
  --sample-column LibraryName --run-column Run --out ../samples.tsv
P=../zabuye-fresh
W=../compute
python workflow/stage.py --mode fresh --project "$P" --conda "$(conda info --base)" \
  --s5 "$BOOKS/Supplementary_Data_S5.zip"
python workflow/run.py reads --samples ../samples.tsv --project "$P" --work "$W" --threads 8
```

The `reads` branch uses fastp 1.3.2 without overlap correction for profiling, binning depth and the new 14-library assemblies. The four historical coassemblies used fastp 0.22.0 with paired-end adapter detection, a Phred threshold of 20 for qualified bases, a minimum read length of 50 bp and overlap correction. The authors’ 9 April log for D1 records the output filenames later read by MEGAHIT. Reprocessing the 35 author-archived raw libraries with these settings reproduced the input read count for every library in the four MEGAHIT logs. Historical FASTQ checksums were unavailable. This comparison used the authors’ raw files and assessed read counts through filtering.

After running `reads`, choose one of the following historical assembly inputs. The default uses the raw pairs in `work/fetch`. To use the authors’ raw pairs instead, add `--historical-raw DIR`. With a newly generated manifest, the workflow checks each reconstructed library against S5 before assembly. These commands print the plan. Add `--execute` in a Linux compute allocation. [Software setup and input filenames](docs/compute.md#executing-the-raw-data-branches)

```bash
python workflow/run.py historical --samples ../samples.tsv --project "$P" --work "$W" \
  --historical-fastp "$FASTP_0220" --assembly-threads 48
# Alternative for users with the archived clean pairs in 00_qc_raw
python workflow/run.py historical --samples ../samples.tsv --project "$P" --work "$W" \
  --assembly-reads "$ASSEMBLY_READS" --assembly-threads 48
```

| Scientific task | Entry point | Main products |
| --- | --- | --- |
| Paired read processing and taxonomic profiling | `run.py reads`, `taxonomy` | fastp reads, two-pass Kraken2/Bracken counts, paired taxonomic models |
| Community gene profiling | `run.py community-ko` | HUMAnN tables, KEGG orthologue abundance and primary-control fits |
| Historical MAG recovery | `run.py historical` | Version-checked read filtering or archived clean pairs, four coassemblies, refined bins, quality scores, dRep representatives and r226 taxonomy |
| Genome-linked traits | `run.py primary-traits`, `reference_steps.py`, `markers.py` | Competitive recruitment, graded module potential, strict carrier abundance, detected-carrier prevalence, pI and marker summaries |
| Identical-input assembly, binning and coverage comparisons | `run.py sensitivity` | 29 assembly tasks, three-binner runs, strategy-specific catalogues, matched groups and Nonpareil diagnostics |
| Organismal phylogeny | `run.py trees`, `run.py replay` | Protein and codon trees, Pagel’s λ and within-lineage PGLS |
| Environmental, carrier and external analyses | `run.py study-models`, `fresh_model_inputs.py`, `external_models.R` | Paired fits, salinity conditioning, additive carrier contrasts and external lake-level models |
| Deposition sensitivity | `deposition_qc.py`, `fcs_matrices.py`, `native/revision/fcs_phylogenetic_sensitivity.R` | FCS-GX reports, fixed-reference ecological models and pruned-tree fits using the reviewed S6 inputs |
| Display and artwork comparison | `figures.py`, `artwork.py` | Numbered redraws, local checksums and optional export of the Word-embedded figures |

The figure runner has a separate source-table entry point. An optional Figure 1a script also draws the two lakes, relief and sampling stations from the archived map layers and S1. Field chemistry, optical measurements, historical reference files and author artwork have their own provenance. The release checks exercised downstream fits and redraws with study inputs. The computationally intensive route from all 35 raw libraries was not rerun end to end for this code release.

For the commands and inputs for each panel, see the [Linux and HPC run guide](docs/compute.md) and [figure index](docs/figures.md). The release page provides the citation for this version. [Release v1.1.0](https://github.com/xiaowu20/zabuye-metagenomics-workflow/releases/tag/v1.1.0) · [Zenodo archive for all versions](https://doi.org/10.5281/zenodo.24455127)
