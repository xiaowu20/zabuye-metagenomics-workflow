# Rebuild the Figure 4, S9 and S8 pI grids from the study's rounded protein pI.
# Each resampling draw in the equal-weight curves samples complete MAGs.
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript protein_density.R PROTEIN_PI_TSV CATALOGUE_TSV_OR_S4_XLSX OUTPUT_DIRECTORY")
}
proteins <- read_tsv(args[1], col_select = c("mag_id", "pI"),
                     col_types = cols(mag_id = col_character(), pI = col_double()))
if (tolower(tools::file_ext(args[2])) == "xlsx") {
  catalogue <- readxl::read_excel(args[2], sheet = "ZBY_proteome_pI")
} else {
  catalogue <- read_tsv(args[2], show_col_types = FALSE)
}
required <- c("mag_id", "gtdb_phylum", "quality_tier")
stopifnot(all(required %in% names(catalogue)), !anyDuplicated(catalogue$mag_id))
if ("source_assembly" %in% names(catalogue)) {
  catalogue <- catalogue[catalogue$source_assembly %in% "ZBY_water", ]
}
if (!nrow(catalogue) || anyNA(catalogue[required]) ||
    any(!nzchar(catalogue$gtdb_phylum)) ||
    !all(catalogue$quality_tier %in% c("HQ", "MQ"))) {
  stop("Check the source, phylum and quality labels for the Zabuye water catalogue.")
}
proteins <- proteins[proteins$mag_id %in% catalogue$mag_id, ]
if (!setequal(proteins$mag_id, catalogue$mag_id) || any(!is.finite(proteins$pI))) {
  stop("Each selected MAG must have finite protein pI records.")
}
output <- file.path(args[3], c("all_MAG_density.tsv", "HQ_density.tsv", "pooled_protein_density.tsv"))
if (any(file.exists(output))) stop("An output grid already exists. Select a new destination.")
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)

# Bytewise MAG-ID order reproduces the archived dplyr arrange(mag_id) draw order.
catalogue <- catalogue[order(catalogue$mag_id, method = "radix"), ]
group <- ifelse(catalogue$gtdb_phylum == "Halobacteriota", "Halobacteriota", "Other phyla")
groups <- c("Halobacteriota", "Other phyla")
if (!all(groups %in% group[catalogue$quality_tier == "HQ"])) stop("Both groups require HQ MAGs.")
per_mag <- split(proteins$pI, factor(proteins$mag_id, levels = catalogue$mag_id))
if (any(lengths(per_mag) < 2)) stop("Each per-MAG density requires at least two proteins.")

grid <- seq(2, 11, length.out = 512)
bandwidth <- bw.nrd0(proteins$pI)
density_args <- list(bw = bandwidth, from = min(grid), to = max(grid), n = length(grid))
# R 4.4 changed the FFT grid spacing. The archived curves use the earlier convention.
if (getRversion() >= "4.4.0") density_args$old.coords <- TRUE
on_grid <- function(values) {
  do.call(stats::density, c(list(x = values), density_args))$y
}
curves <- t(vapply(per_mag, on_grid, numeric(length(grid))))
group_curve <- function(label, keep) {
  subset <- curves[group == label & keep, , drop = FALSE]
  size <- nrow(subset)
  resampled <- replicate(1999, {
    sample_ids <- sample.int(size, size, replace = TRUE)
    colMeans(subset[sample_ids, , drop = FALSE])
  })
  data.frame(group = label, x = grid, y = colMeans(subset),
             lo = apply(resampled, 1, quantile, probs = 0.025, names = FALSE),
             hi = apply(resampled, 1, quantile, probs = 0.975, names = FALSE))
}
both_groups <- function(keep) {
  do.call(rbind, lapply(groups, group_curve, keep = keep))
}
# The archive drew the two complete-catalogue groups before the two HQ groups.
set.seed(20260724)
all_mag <- both_groups(rep(TRUE, nrow(catalogue)))
hq <- both_groups(catalogue$quality_tier == "HQ")
pooled <- do.call(rbind, lapply(groups, function(label) {
  pI <- proteins$pI[proteins$mag_id %in% catalogue$mag_id[group == label]]
  data.frame(group = label, x = grid, y = on_grid(pI))
}))
write_tsv(all_mag, output[1])
write_tsv(hq, output[2])
write_tsv(pooled, output[3])
message(sprintf("Wrote three pI grids for %d proteins in %d MAGs.", nrow(proteins), nrow(catalogue)))
