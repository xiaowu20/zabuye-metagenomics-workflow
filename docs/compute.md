# Computing the study analyses

## Software, references and paths

Use Linux for the sequencing steps. Model replay and table-based figure rendering can also run in a local Python/R installation. The replay environment file specifies dependencies. It is not a solver lock or a container image. Python ≥3.11 and the packages in `env/requirements.txt` are sufficient for the Python-only commands. Run `Rscript --vanilla -e 'sessionInfo()'` and retain software versions, reference releases and the staged checksums with a run.

| Step | Study software or reference |
| --- | --- |
| Retrieval and QC | NCBI SRA Toolkit (`prefetch`, `fasterq-dump`), fastp 0.22.0 for the four historical coassemblies and fastp 1.3.2 for subsequent read processing |
| Assembly and binning | MEGAHIT 1.2.9, metaSPAdes 4.0.0, Bowtie2, SAMtools, MetaBAT2 2.18, MaxBin2 2.2.7, CONCOCT and its auxiliary scripts, DAS Tool 1.1.7 |
| Genome quality and redundancy | CheckM2 1.1.0 with the matching DIAMOND database, dRep 3.6.2, fastANI |
| Read and genome annotation | Kraken2 2.1.3, Bracken 3.0.1, PlusPF `k2_pluspf_20240605`, Prodigal, eggNOG-mapper 2.1.12 with data 5.0.2 and DIAMOND 2.1.24 |
| Community genes and recruitment | HUMAnN 3.9, MetaPhlAn 4.2.4 with `mpa_vJun23_CHOCOPhlAnSGB_202307`, ChocoPhlAn and UniRef90, CoverM 0.7.0 |
| Taxonomy and trees | Study catalogue GTDB-Tk 2.6.1/r226, external classification and fresh codon-marker workflow 2.5.1/r226, older organismal alignment r220, IQ-TREE 3.1.1, MAFFT, PAL2NAL |
| Coverage and protein traits | Nonpareil 3.5.5 with its R package, Biopython 1.87, KOfamScan and HMMER with the KOfam and Pfam-A databases |
| Statistical analysis | R 4.5.1, vegan, permute, nlme, lme4, ape 5.8-1, phangorn 2.12.1 and phytools 2.5-2. The read-based genus analysis also uses MaAsLin2. |

The `primary.py` wrapper resolves programs on `PATH`. Activate the relevant environment for each stage, or expose the required binaries with your site's modules. The revision templates use named environments relative to the **Conda base** provided to `stage.py --conda`. Install tools in the environments below with conda-forge and Bioconda, or point the corresponding paths at equivalent installations. Keep the old and current GTDB reference releases in separate directories.

| Named environment | Binaries required by the revision templates |
| --- | --- |
| `assembly` | MEGAHIT, SPAdes, CoverM, fastANI and Prodigal. For the reference build include Kraken2, Bracken and their database-build dependencies. |
| `binning` | Bowtie2, SAMtools, MetaBAT2 including `jgi_summarize_bam_contig_depths`, MaxBin2, CONCOCT, DAS Tool, dRep and their executable dependencies |
| `checkm2` | CheckM2 and DIAMOND |
| `eggnog` | eggNOG-mapper and DIAMOND |
| `gtdbtk_v25` | GTDB-Tk 2.5.1 and its dependencies |
| `iqtree` | IQ-TREE, MAFFT and PAL2NAL; make the `iqtree` executable available under that name |

On the original compute node CONCOCT had a separate Python package path containing setuptools 80.9.0 for `pkg_resources`. With a compatible environment this package can reside in `binning`. To reproduce the separate path, use the binning interpreter after staging, with `REV="$P/Major_revised_FM/R1_20261110"`.

```bash
"$(conda info --base)/envs/binning/bin/python" -m pip install --no-deps \
  --target "$REV/envs/concoct_compat" setuptools==80.9.0
```

Install databases outside the checkout. Set these variables in the compute allocation and point each one to the actual local installation. Binary databases and third-party reference genomes retain their upstream distribution terms.

| Variable for primary wrapper | Input |
| --- | --- |
| `CHECKM2_DB` | `uniref100.KO.1.dmnd` file |
| `GTDBTK_DATA_PATH` | GTDB r226 database root for the current classification |
| `EGGNOG_DATA_DIR` | eggNOG database directory with the DIAMOND and SQLite files |
| `PLUSPF_DB`, `REGIONAL_DB` | Kraken2 database directories with Bracken 150-bp distribution files |
| `CHOCO_DB`, `UNIREF_DB`, `METAPHLAN_DB` | HUMAnN nucleotide/protein reference directories and the MetaPhlAn marker directory |

Fresh revision scripts also look for `$P/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd`, `$P/db/eggnog` and `$P/db/gtdbtk_r226`. Make those resources accessible at the staged paths. The Qt-free figure route uses Matplotlib and R graphics. Install `Maaslin2` in the R environment for fresh genus modelling; the frozen figure and phylogenetic routes can run independently. The Nonpareil executable and R library are needed to refit raw coverage curves.

## Executing the raw-data branches

Define `P`, `W`, `BOOKS` and the manifest as in the README. Use a new fresh-mode project and a separate work root, with sufficient storage for FASTQs, indexes, BAMs and assemblies. For native shell templates, choose Linux project and Conda paths without spaces. Filesystem symlinks, hard links within the work volume and node-local temporary storage are used. `run.py` executes one command at a time inside the allocation. Calling the underlying array templates through your site's scheduler is also possible. Some native templates request or use 32 threads explicitly. The historical jobs' resource requests are visible in their template headers.

Install the two fastp versions separately, or use equivalent site installations. `reads` uses the fastp on `PATH`. `--historical-fastp` selects the executable for the four original assemblies. Both steps check the executable version.

```bash
QC_ENVS="$(pwd)/../qc-envs"
conda create --prefix "$QC_ENVS/fastp-0.22.0" --strict-channel-priority -c conda-forge -c bioconda fastp=0.22.0
conda create --prefix "$QC_ENVS/fastp-1.3.2" --strict-channel-priority -c conda-forge -c bioconda fastp=1.3.2
FASTP_0220="$QC_ENVS/fastp-0.22.0/bin/fastp"
export PATH="$QC_ENVS/fastp-1.3.2/bin:$PATH"
```

```bash
common=(--samples ../samples.tsv --project "$P" --work "$W")
python workflow/run.py reads "${common[@]}" --threads 8 --execute
python workflow/run.py historical "${common[@]}" --historical-fastp "$FASTP_0220" --assembly-threads 48 --threads 32 --execute
python workflow/run.py primary-traits "${common[@]}" --execute
python workflow/run.py sensitivity "${common[@]}" --execute
python workflow/run.py trees "${common[@]}" --execute
```

The historical command above starts from the raw pairs in `$W/fetch`. If the original raw archive is available, add `--historical-raw DIR` with files named `<sample_id>.R1.raw.fastq.gz` and `<sample_id>.R2.raw.fastq.gz`. Alternatively, replace `--historical-fastp "$FASTP_0220"` with `--assembly-reads DIR` to use archived `00_qc_raw` files named `<sample_id>.R1.clean.fastq.gz` and `<sample_id>.R2.clean.fastq.gz` directly. The `--provider` spelling remains a legacy alias for this directory. Run `reads` before either route to supply the later depth-mapping reads. Use a separate work directory when comparing different sources.

`manifest.py` takes `historical_input_reads` from each original MEGAHIT log in S5. This is the number of individual reads across both mates. For reconstructed pairs the workflow compares this value to `summary.after_filtering.total_reads` in the fastp JSON and stops before coassembly on a mismatch. Regenerate a v1.0.0 five-column manifest with the README command to enable this check. The original-archive replay matched all 35 counts. The FASTQs exported from SRA were outside that validation, and no historical FASTQ digest was available for comparison.

1. **Read processing.** Each SRA run supplies a pair in `$W/fetch/<sample_id>/`. The `reads` branch runs fastp 1.3.2 with paired-end adapter detection, a Phred threshold of 20 for qualified bases and a minimum read length of 50 bp. These reads, filtered without `--correction`, go to `$W/trim/` for profiling, recruitment, binning depth and the new 14-library assemblies.
2. **Four historical grouped assemblies.** The `historical-trim` step separately runs fastp 0.22.0 with `--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --correction`. The original run used eight threads. The wrapper uses eight when allocated at least that many. Filtered pairs and reports go to `$W/historical-trim/<sample_id>/`. MEGAHIT receives the pairs in the order recorded in S5 and used 48 threads in the original assemblies. MetaBAT2 used 32 threads. Historical MaxBin2 input uses `-abund` with total mean depth. The new independent/pooled-14 assemblies use per-library `-abund_list`. DAS Tool uses a 0.5 score threshold.
3. **MAG-linked response matrices.** Following quality assessment and dereplication, the primary branch calls proteins, assigns orthologues, applies the corrected empty-member scoring rules and recruits paired reads competitively. At mean depth ≥1 and breadth ≥0.50, abundance-weighted GMP uses graded module scores, SCA weights strict carriers by abundance, and detected-carrier prevalence gives each detected genome equal weight. Alternatives at breadth 0.10 and 0.25 are in the S5 catalogue branch. Protein pI is calculated after removing nonstandard residues and excluding sequences shorter than 20 residues. Historical acidic/basic “fraction” columns store percentages.
4. **Sensitivity.** The 29 tasks are 14 independent MEGAHIT, 14 independent metaSPAdes and one pooled-14 MEGAHIT job on the same fastp libraries. The sensitivity branch also re-bins original contigs, compares reciprocal ANI and aligned fractions, and includes both the historical all-reference recruitment and the restricted remap of the same 269 genomes. It requires the historical contigs, original bins/BAMs and analytical FASTAs at the paths created by `joins.py reference`. An assembly reconstructed with different representatives cannot be silently substituted for the archived 269. A stand-alone task can be run from its printed command with a valid node allocation.
5. **Trees and marker checks.** The original 143-tip r220 archaeal alignment came from an earlier collection. The archived tree in S5 supplies the pruned within-lineage baseline for replay. To rerun its inference from the original masked alignment, use the command below. The fresh 84-genome workflow compares the masked protein tree, a matched amino-acid tree and a codon-partitioned nucleotide tree to that baseline. `markers.py` handles the distinct potassium-marker gene tree. The study selected LG+F+G4 for the 84-tip masked tree and LG+F+R7 for the 187-tip marker tree. Protein sequence pI and the quality-adjusted PGLS refer to different response scales from gene activity or ion flux.

```bash
python workflow/reference_steps.py archival-tree --msa "$ARCHIVED_R220_AR53" --out ../archival-tree
python workflow/reference_steps.py kofam --query "$CANDIDATE_PROTEINS" \
  --profiles "$KOFAM/profiles" --ko-list "$KOFAM/ko_list" --out ../kofam
python workflow/markers.py --proteins "$W/proteins/project/proteins.faa" \
  --gff-dir "$W/proteins/project" --kofam-hits ../kofam/kofam_detail_output.tsv \
  --pfam-library "$PFAM_A" --pfam-candidates "$PFAM_CANDIDATES" --tree --out ../markers
```

Use candidate FASTAs with source-qualified protein IDs. The KOfam query and the Pfam assessment subset are distinct archival inputs. The full Pfam-A search uses the full-sequence reporting threshold. Preserve the reference release and identity crosswalk with the results.

## Read-based species and community genes

Prepare TPMC and QXLSG regional reference genomes and their original reference workbooks. The staged regional builder takes `TPMC_MAG_summary.xlsx` and `QXLSG_mmc2.xlsx`, with the matching FASTA directories. TPMC uses its published genome-species clusters; the regional QXLSG library uses the published species representatives. The external proteome comparison uses the larger QXLSG catalogue. [TPMC](https://doi.org/10.1038/s41467-024-45895-8) · [QXLSG](https://doi.org/10.1016/j.celrep.2025.116483)

```bash
T="$P/dual_species_reanalysis_20260819"
python "$T/code/00_prepare_regional_taxonomy.py" --input-dir "$REFERENCE_BOOKS" \
  --output-root "$T" --tpmc-genomes "$TPMC_FASTAS" --qxlsg-genomes "$QXLSG_FASTAS"
bash "$T/hpc/01_build_regional_species_db.slurm"
export REGIONAL_DB="$T/db/regional_species"
python workflow/run.py taxonomy "${common[@]}" --local-inputs "$LOCAL" --execute
python workflow/run.py community-ko "${common[@]}" --local-inputs "$LOCAL" --out ../fresh-ko --execute
```

The taxonomy branch applies PlusPF first, passes residual reads to the regional library and estimates both ranks with Bracken. Region-specific features keep their provenance in the merged matrix. The staged taxonomic R program includes the original 19-library descriptive display. Final 21-library display inputs are reconstructed from the paper workbooks by `paper_inputs.py`. `$LOCAL/taxonomy/manifest_with_site.tsv` supplies recorded site and sampling metadata. For community KOs, HUMAnN gene families are joined, regrouped to KOs and renormalised to copies per million. The optional model step takes the newly computed KO table and the design/QC controls in `$LOCAL/community_ko`. The latter must correspond to that mapping run.

## Statistical checkpoints and figures

```bash
python workflow/fresh_model_inputs.py --samples ../samples.tsv --work "$W" \
  --design "$LOCAL/ions" --s4 "$BOOKS/Supplementary_Data_S4.xlsx" --out ../fresh-model-inputs
Rscript --vanilla workflow/paired_models.R ../fresh-model-inputs ../paired
Rscript --vanilla workflow/ion_models.R ../fresh-model-inputs ../environment
Rscript --vanilla workflow/incidence.R ../fresh-model-inputs ../incidence
Rscript --vanilla workflow/crosslake.R ../fresh-model-inputs ../crosslake
python workflow/run.py study-models --project "$P" --books "$BOOKS" \
  --local-inputs "$LOCAL" --out ../study-fits --execute
```

`fresh_model_inputs.py` takes newly computed matrices from `work/responses` and joins them to the measured site chemistry, hydrology and sampling design. The study checkpoint route uses the author's archived `ions`, `community_ko`, `dom_engines`, `proteome`, `carrier` and `figure3/models` directories. These carry sample/covariate matrices, including historical alternatives used to assess software and annotation sensitivity. Keep the 15 DOM-matched libraries, the historical 19-library sensitivity and the seven independent station pairs separate. R 4.3.3/vegan 2.7-1 and R 4.5.1/vegan 2.7-3 are distinct DOM replay environments. Record the actual engine with each result.

`carrier_attribution.py` reconstructs the additive MAG and genus terms and checks their sum against the aggregate. It requires per-MAG depths, breadths, abundances, annotations and the reference contrasts named in the script. `figure3_tables.py` joins those outputs to the model tables and prepares all four panels. To draw the recomputed Figure 3, pass the parent of `figure3/source_data` as the local-inputs root to `figures.py --group carriers` (for the command above this root is `../study-fits`). Other plotted results use the figure-specific checkpoints in the [figure index](figures.md).

## External sediment comparison and curation

For table-level replay, `Rscript --vanilla workflow/external_models.R "$BOOKS/Supplementary_Data_S4.xlsx" ../lake-fits` refits the equal-weight lake OLS, the lake-intercept ML model in nlme and the crossed lake/phylum ML model in lme4, with a one-df likelihood-ratio test for the latter slope. Both OLS and lake-intercept leave-one-lake-out series are written to one labelled table. Lake means supply the ten environmental observations.

For sequence-based reconstruction, obtain the published QXLSG MAGs, predicted proteins and source workbook under the upstream terms. Stage the proteins as `$P/db/refdb/qxlsg_faa/*.faa`, the renamed MAGs as `$P/db/qxlsg/fasta_renamed/*.fa` and the original workbook as `$P/db/qxlsg/QXLSG_mmc2.xlsx`. Use GTDB-Tk 2.5.1 and r226 for this branch. The lake attribution step excludes cross-sample MAGs without a unique lake assignment. [NGDC project PRJCA037687](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA037687) · [QXLSG source](https://doi.org/10.1016/j.celrep.2025.116483)

```bash
python workflow/reference_steps.py external-taxonomy --genomes "$P/db/qxlsg/fasta_renamed" \
  --out "$P/db/qxlsg/gtdbtk_out_v25"
python "$P/external/24_qxlsg_pI_calculation.py"
E="$P/packages/pkg_CC_protein_pI_qxlsg"
python workflow/summarise_proteins.py "$E/protein_pI_summary.tsv" "$E/MAG_level_pI_summary_qxlsg.tsv"
Rscript --vanilla "$P/external/25_qxlsg_pI_statistics.R"
Rscript --vanilla workflow/external_models.R "$E/pI_with_full_taxonomy_lake_v3.tsv" ../lake-fits-from-proteins
```

For FCS-GX screening, `workflow/deposition_qc.py --help` describes the analytical FASTA, digest and taxid manifest. This operation produces reports. Subsequent contig dispositions require the reviewed decisions and reference sequences. Supplying S6 to `stage.py` installs its original follow-up scripts and frozen sensitivity tables in the study layout. S27/S28 use that frozen stage. Deposition copies and archived analytical references must keep their identities throughout.

The S6 ecological sensitivity changes the retained MAGs or their functional annotations while holding the original mapping reference fixed. The public matrix script compares its baseline with all twelve archived response matrices. Include `--fit` to rerun the four catalogues' exact tests, dispersion checks and module analyses. The phylogenetic script separately refits each archived tree after pruning the flagged genome identified by the local action table. Both commands write to new directories.

```bash
REV="$P/Major_revised_FM/R1_20261110"
python workflow/fcs_matrices.py --project "$P" --out ../s6-models --fit
Rscript --vanilla native/revision/fcs_phylogenetic_sensitivity.R "$REV" ../s6-trees
Rscript --vanilla plots/fcs_inference_figures.R "$REV" ../s6-figure ../s6-models ../s6-trees
```
