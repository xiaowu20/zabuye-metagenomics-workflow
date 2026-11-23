# Corrected HUMAnN3 sample-by-KO analysis for the ZBY water sampling design.
# Primary inference uses 14 paired spring-autumn station libraries. The n = 17
# sensitivity set adds the autumn-only N5, S4 and S5 libraries. Summer and
# winter basin composites are excluded from inferential models.

set.seed(20260814)

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2 || length(args)>3) stop("Usage: Rscript community_ko.R INPUTS OUTPUT [KO_CPM]")
input_dir <- args[1]; provenance_dir <- input_dir; humann_dir <- args[2]
result_dir <- file.path(humann_dir, "results")
figure_dir <- file.path(humann_dir, "figures")
log_dir <- file.path(humann_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  ko = file.path(input_dir, "HUMAnN3_sample_by_KO_CoPM_35samples.tsv"),
  modules = file.path(input_dir, "HUMAnN3_module_abundance_CoPM_35samples.tsv"),
  panel = file.path(input_dir, "KO_panel_members_v3.csv"),
  taxon = file.path(input_dir, "HUMAnN3_taxonomic_resolution_legacy_summary.tsv"),
  legacy_qc = file.path(input_dir, "MAG_catalogue_detection_and_CoverM_mapping_QC.tsv"),
  membership = file.path(provenance_dir, "ZBY_water_analysis_membership.tsv"),
  master = file.path(provenance_dir, "Metagenome_Sample_Master_35.tsv")
)
if (length(args)==3) files["ko"] <- args[3]
missing_files <- files[!file.exists(files)]
if (length(missing_files)) stop("Missing required files:\n", paste(missing_files, collapse = "\n"))

read_tsv_base <- function(path, ...) {
  read.delim(path, check.names = FALSE, quote = "", comment.char = "", ...)
}

ko <- read_tsv_base(files[["ko"]])
names(ko)[1] <- "feature"
sample_names <- sub("_Abundance-RPKs$", "", names(ko)[-1])
names(ko)[-1] <- sample_names
ko_values <- as.matrix(ko[-1])
storage.mode(ko_values) <- "double"
rownames(ko_values) <- ko$feature

panel <- read.csv(files[["panel"]], check.names = FALSE, stringsAsFactors = FALSE)
panel_map <- unique(panel[c("KO", "score_target")])
panel_kos <- unique(panel_map$KO)

community_mask <- !grepl("\\|", ko$feature, fixed = FALSE)
community_features <- ko$feature[community_mask]
community_values <- ko_values[community_mask, , drop = FALSE]
community_ko_mask <- grepl("^K[0-9]{5}$", community_features)
community_kos <- community_features[community_ko_mask]
community_ko_values <- community_values[community_ko_mask, , drop = FALSE]

present_panel <- intersect(panel_kos, rownames(community_ko_values))
missing_panel <- setdiff(panel_kos, present_panel)
panel_values <- matrix(
  0, nrow = length(panel_kos), ncol = length(sample_names),
  dimnames = list(panel_kos, sample_names)
)
panel_values[present_panel, ] <- community_ko_values[present_panel, , drop = FALSE]

# Recompute the sample-by-module matrix directly from the archived sample-by-KO
# CoPM table. This reproduces the original sum-over-member-KOs aggregation.
module_targets <- sort(unique(panel_map$score_target))
module_matrix <- matrix(
  0, nrow = length(sample_names), ncol = length(module_targets),
  dimnames = list(sample_names, module_targets)
)
for (target in module_targets) {
  target_kos <- unique(panel_map$KO[panel_map$score_target == target])
  module_matrix[, target] <- colSums(panel_values[target_kos, , drop = FALSE])
}
module_output <- data.frame(sample_id = rownames(module_matrix), module_matrix, check.names = FALSE)
write.table(
  module_output,
  file.path(result_dir, "HUMAnN3_module_abundance_recomputed_35samples.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

legacy_modules <- read_tsv_base(files[["modules"]])
legacy_ids <- as.character(legacy_modules[[1]])
legacy_matrix <- as.matrix(legacy_modules[-1])
storage.mode(legacy_matrix) <- "double"
rownames(legacy_matrix) <- legacy_ids
overlap_modules <- intersect(colnames(legacy_matrix), colnames(module_matrix))
overlap_samples <- intersect(rownames(legacy_matrix), rownames(module_matrix))
legacy_aligned <- legacy_matrix[overlap_samples, overlap_modules, drop = FALSE]
recomputed_aligned <- module_matrix[overlap_samples, overlap_modules, drop = FALSE]
validation <- data.frame(
  n_samples = length(overlap_samples),
  n_modules = length(overlap_modules),
  max_absolute_difference = max(abs(legacy_aligned - recomputed_aligned)),
  mean_absolute_difference = mean(abs(legacy_aligned - recomputed_aligned)),
  detected_panel_KOs = length(present_panel),
  missing_panel_KOs = length(missing_panel),
  missing_KO_list = paste(missing_panel, collapse = ";")
)
write.table(
  validation,
  file.path(result_dir, "HUMAnN3_module_recompute_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Sample-level HUMAnN3 quality metrics. CoPM totals are calculated only from
# unstratified rows to avoid double counting taxon-stratified subdivisions.
feature_index <- setNames(seq_along(community_features), community_features)
get_row <- function(feature) {
  if (!feature %in% names(feature_index)) return(rep(0, length(sample_names)))
  community_values[feature_index[[feature]], ]
}
unstratified_total <- colSums(community_values)
unmapped <- get_row("UNMAPPED")
ungrouped <- get_row("UNGROUPED")
grouped_ko_copm <- colSums(community_ko_values)
sample_qc <- data.frame(
  sample_id = sample_names,
  unstratified_total_CoPM = unstratified_total,
  unmapped_CoPM = unmapped,
  ungrouped_CoPM = ungrouped,
  grouped_KO_CoPM = grouped_ko_copm,
  mapped_fraction_pct = 100 * (1 - unmapped / unstratified_total),
  grouped_KO_fraction_pct = 100 * grouped_ko_copm / unstratified_total,
  detected_community_KOs = colSums(community_ko_values > 0),
  detected_panel_KOs = colSums(panel_values > 0),
  panel_KO_CoPM = colSums(panel_values),
  stringsAsFactors = FALSE
)

master <- read_tsv_base(files[["master"]], na.strings = c("", "NA"))
master$sample_id <- as.character(master$sample_id)
sample_qc <- merge(
  sample_qc,
  master[c("sample_id", "lake", "matrix", "basin", "season_or_campaign", "physical_sampling_unit")],
  by = "sample_id", all.x = TRUE, sort = FALSE
)
sample_qc <- sample_qc[match(sample_names, sample_qc$sample_id), ]
write.table(
  sample_qc,
  file.path(result_dir, "HUMAnN3_KO_QC_35samples.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

membership <- read_tsv_base(files[["membership"]], na.strings = c("", "NA"))
membership$sample_id <- as.character(membership$sample_id)
membership$basin <- factor(membership$basin, levels = c("north", "south"))
membership$season <- factor(membership$season, levels = c("autumn", "spring"))
membership$site <- factor(membership$site_corrected)
legacy_qc <- read_tsv_base(files[["legacy_qc"]])
legacy_qc$sample_id <- as.character(legacy_qc$sample_id)

ids14 <- membership$sample_id[membership$included_n14_primary]
ids17 <- membership$sample_id[membership$included_n17_sensitivity]
stopifnot(length(ids14) == 14L, length(ids17) == 17L)

make_meta <- function(ids) {
  m <- membership[match(ids, membership$sample_id), ]
  counts <- table(m$site)
  m$coverage_pattern <- ifelse(counts[as.character(m$site)] == 2L, "paired", "autumn_only")
  m$coverage_pattern <- factor(m$coverage_pattern, levels = c("paired", "autumn_only"))
  qc <- sample_qc[match(m$sample_id, sample_qc$sample_id), ]
  old_qc <- legacy_qc[match(m$sample_id, legacy_qc$sample_id), ]
  stopifnot(!anyNA(qc$mapped_fraction_pct), !anyNA(qc$detected_panel_KOs))
  stopifnot(!anyNA(old_qc$catalogue_detection_fraction), !anyNA(old_qc$mapping_rate_pct))
  m$detected_panel_KOs_z <- as.numeric(scale(qc$detected_panel_KOs))
  m$mapped_fraction_pct_z <- as.numeric(scale(qc$mapped_fraction_pct))
  m$MAG_catalogue_detection_fraction_z <- as.numeric(scale(old_qc$catalogue_detection_fraction))
  m$CoverM_mapping_rate_pct_z <- as.numeric(scale(old_qc$mapping_rate_pct))
  m
}

hat_matrix <- function(x) {
  q <- qr.Q(qr(x))
  r <- qr(x)$rank
  if (r == 0L) return(matrix(0, nrow(x), nrow(x)))
  tcrossprod(q[, seq_len(r), drop = FALSE])
}

gower_matrix <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  j <- diag(n) - matrix(1 / n, n, n)
  -0.5 * j %*% (distance_matrix ^ 2) %*% j
}

adjusted_r2 <- function(r2, n, predictor_df) {
  if (n - predictor_df - 1 <= 0) return(NA_real_)
  1 - (1 - r2) * (n - 1) / (n - predictor_df - 1)
}

term_stat_from_g <- function(g, meta, term, tech_terms = character()) {
  predictors <- c(tech_terms, "season", "basin")
  reduced_predictors <- setdiff(predictors, term)
  full_formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  reduced_formula <- if (length(reduced_predictors)) {
    as.formula(paste("~", paste(reduced_predictors, collapse = " + ")))
  } else {
    ~1
  }
  x_full <- model.matrix(full_formula, data = meta)
  x_reduced <- model.matrix(reduced_formula, data = meta)
  h_full <- hat_matrix(x_full)
  h_reduced <- hat_matrix(x_reduced)
  n <- nrow(g)
  identity_matrix <- diag(n)
  ss_total <- sum(diag(g))
  ss_term <- sum(diag((h_full - h_reduced) %*% g))
  ss_residual <- sum(diag((identity_matrix - h_full) %*% g))
  rank_full <- qr(x_full)$rank
  rank_reduced <- qr(x_reduced)$rank
  df_term <- rank_full - rank_reduced
  df_residual <- n - rank_full
  pseudo_f <- (ss_term / df_term) / (ss_residual / df_residual)
  r2_term <- ss_term / ss_total
  r2_full <- sum(diag(h_full %*% g)) / ss_total
  r2_reduced <- sum(diag(h_reduced %*% g)) / ss_total
  adj_full <- adjusted_r2(r2_full, n, rank_full - 1L)
  adj_reduced <- adjusted_r2(r2_reduced, n, rank_reduced - 1L)
  c(
    pseudo_F = pseudo_f,
    R2 = r2_term,
    adjusted_R2 = adj_full - adj_reduced,
    df_term = df_term,
    df_residual = df_residual,
    full_model_R2 = r2_full,
    full_model_adjusted_R2 = adj_full
  )
}

season_label_permutations <- function(meta) {
  paired_sites <- unique(as.character(meta$site[meta$coverage_pattern == "paired"]))
  combinations <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(paired_sites))))
  out <- vector("list", nrow(combinations))
  for (i in seq_len(nrow(combinations))) {
    permuted <- meta
    permuted$season <- as.character(permuted$season)
    for (j in seq_along(paired_sites)) {
      if (combinations[i, j]) {
        index <- which(as.character(permuted$site) == paired_sites[j])
        permuted$season[index] <- rev(permuted$season[index])
      }
    }
    permuted$season <- factor(permuted$season, levels = c("autumn", "spring"))
    out[[i]] <- permuted
  }
  out
}

basin_label_permutations <- function(meta) {
  site_frame <- unique(meta[, c("site", "basin", "coverage_pattern")])
  site_frame[] <- lapply(site_frame, as.character)
  patterns <- unique(site_frame$coverage_pattern)
  choices <- lapply(patterns, function(pattern) {
    subset <- site_frame[site_frame$coverage_pattern == pattern, ]
    n_north <- sum(subset$basin == "north")
    combn(subset$site, n_north, simplify = FALSE)
  })
  grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    north_sites <- unlist(mapply(
      function(options, index) options[[index]], choices,
      as.integer(grid[i, ]), SIMPLIFY = FALSE
    ))
    permuted <- meta
    permuted$basin <- factor(
      ifelse(as.character(permuted$site) %in% north_sites, "north", "south"),
      levels = c("north", "south")
    )
    out[[i]] <- permuted
  }
  out
}

exact_term_test <- function(g, meta, term, tech_terms = character()) {
  observed <- term_stat_from_g(g, meta, term, tech_terms)
  permutations <- if (term == "season") season_label_permutations(meta) else basin_label_permutations(meta)
  permuted_f <- vapply(
    permutations,
    function(permuted) term_stat_from_g(g, permuted, term, tech_terms)[["pseudo_F"]],
    numeric(1)
  )
  p_exact <- sum(permuted_f >= observed[["pseudo_F"]] - sqrt(.Machine$double.eps)) / length(permuted_f)
  list(observed = observed, p_exact = p_exact, n_permutations = length(permuted_f))
}

bootstrap_indices <- function(meta) {
  site_frame <- unique(meta[, c("site", "basin", "coverage_pattern")])
  site_frame$site <- as.character(site_frame$site)
  strata <- interaction(site_frame$basin, site_frame$coverage_pattern, drop = TRUE)
  sampled_sites <- unlist(lapply(split(site_frame$site, strata), function(z) sample(z, length(z), replace = TRUE)))
  index <- integer(0)
  boot_site <- character(0)
  counter <- 0L
  for (site_id in sampled_sites) {
    counter <- counter + 1L
    rows <- which(as.character(meta$site) == site_id)
    index <- c(index, rows)
    boot_site <- c(boot_site, rep(paste0("boot_site_", counter), length(rows)))
  }
  list(index = index, site = boot_site)
}

bootstrap_term <- function(distance_matrix, meta, term, tech_terms = character(), n_boot = 999L) {
  estimates <- rep(NA_real_, n_boot)
  for (iteration in seq_len(n_boot)) {
    draw <- bootstrap_indices(meta)
    resampled_meta <- meta[draw$index, , drop = FALSE]
    resampled_meta$site <- factor(draw$site)
    resampled_distance <- distance_matrix[draw$index, draw$index, drop = FALSE]
    stat <- tryCatch(
      term_stat_from_g(gower_matrix(resampled_distance), resampled_meta, term, tech_terms),
      error = function(error) NULL
    )
    if (!is.null(stat)) estimates[iteration] <- stat[["R2"]]
  }
  estimates[is.finite(estimates)]
}

dispersion_exact_test <- function(distance_matrix, meta, term) {
  get_f <- function(current_meta) {
    group <- if (term == "season") current_meta$season else current_meta$basin
    fit <- suppressWarnings(betadisper(as.dist(distance_matrix), group))
    unname(anova(fit)[1, "F value"])
  }
  observed <- get_f(meta)
  permutations <- if (term == "season") season_label_permutations(meta) else basin_label_permutations(meta)
  permuted <- vapply(permutations, get_f, numeric(1))
  data.frame(
    term = term,
    dispersion_F = observed,
    dispersion_P_exact = sum(permuted >= observed - sqrt(.Machine$double.eps)) / length(permuted),
    exact_permutations = length(permuted)
  )
}

model_specs <- list(
  M0_unadjusted = character(),
  M1_MAG_catalogue_detection = "MAG_catalogue_detection_fraction_z",
  M2_CoverM_mapping_rate = "CoverM_mapping_rate_pct_z",
  M3_HUMAnN_detected_panel_KOs = "detected_panel_KOs_z",
  M4_HUMAnN_mapped_fraction = "mapped_fraction_pct_z",
  M5_HUMAnN_both = c("detected_panel_KOs_z", "mapped_fraction_pct_z")
)
scopes <- list(paired14_primary = ids14, individual17_sensitivity = ids17)
results <- list()
dispersion_results <- list()
counter <- 0L
dispersion_counter <- 0L

for (scope_name in names(scopes)) {
  ids <- scopes[[scope_name]]
  meta <- make_meta(ids)
  response <- module_matrix[ids, , drop = FALSE]
  variable_columns <- apply(response, 2, function(x) is.finite(sd(x)) && sd(x) > 0)
  response <- response[, variable_columns, drop = FALSE]
  distance_matrix <- as.matrix(vegdist(response, method = "bray"))
  gower <- gower_matrix(distance_matrix)

  for (model_name in names(model_specs)) {
    tech_terms <- model_specs[[model_name]]
    for (term in c("season", "basin")) {
      test <- exact_term_test(gower, meta, term, tech_terms)
      boot <- bootstrap_term(distance_matrix, meta, term, tech_terms, n_boot = 999L)
      counter <- counter + 1L
      results[[counter]] <- data.frame(
        response = "HUMAnN3 community KO module abundance",
        sample_scope = scope_name,
        n_samples = nrow(meta),
        n_sites = length(unique(meta$site)),
        n_variable_modules = ncol(response),
        model = model_name,
        technical_covariates = if (length(tech_terms)) paste(tech_terms, collapse = ";") else "none",
        term = term,
        pseudo_F = test$observed[["pseudo_F"]],
        R2 = test$observed[["R2"]],
        adjusted_R2 = test$observed[["adjusted_R2"]],
        CI95_low = unname(quantile(boot, 0.025, na.rm = TRUE)),
        CI95_high = unname(quantile(boot, 0.975, na.rm = TRUE)),
        P_exact = test$p_exact,
        exact_permutations = test$n_permutations,
        bootstrap_replicates = length(boot),
        permutation_unit = if (term == "season") {
          "within-site spring-autumn label swaps; autumn-only observations fixed"
        } else {
          "whole-site basin labels; paired/autumn-only coverage preserved"
        },
        stringsAsFactors = FALSE
      )
    }
  }

  for (term in c("season", "basin")) {
    dispersion_counter <- dispersion_counter + 1L
    dispersion_results[[dispersion_counter]] <- cbind(
      data.frame(sample_scope = scope_name, n_samples = nrow(meta), stringsAsFactors = FALSE),
      dispersion_exact_test(distance_matrix, meta, term)
    )
  }
}

results <- do.call(rbind, results)
dispersion_results <- do.call(rbind, dispersion_results)
write.table(
  results,
  file.path(result_dir, "HUMAnN3_exact_PERMANOVA_n14_n17.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  dispersion_results,
  file.path(result_dir, "HUMAnN3_dispersion_n14_n17.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# PCoA plotting data for the n = 17 sensitivity set.
meta17 <- make_meta(ids17)
response17 <- module_matrix[ids17, , drop = FALSE]
response17 <- response17[, apply(response17, 2, sd) > 0, drop = FALSE]
ordination <- cmdscale(vegdist(response17, method = "bray"), k = 2, eig = TRUE, add = TRUE)
explained <- 100 * ordination$eig[1:2] / sum(ordination$eig[ordination$eig > 0])
pcoa <- data.frame(
  sample_id = rownames(ordination$points),
  PCoA1 = ordination$points[, 1],
  PCoA2 = ordination$points[, 2],
  stringsAsFactors = FALSE
)
pcoa <- merge(pcoa, meta17[c("sample_id", "basin", "season", "site_corrected", "coverage_pattern")], by = "sample_id", sort = FALSE)
pcoa <- pcoa[match(ids17, pcoa$sample_id), ]
write.table(pcoa,file.path(result_dir,"HUMAnN3_PCoA_n17.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(log_dir,"sessionInfo.txt"))
