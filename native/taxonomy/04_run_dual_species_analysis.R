# Reanalyse the read-taxonomy layer after sequential PlusPF + regional-MAG classification.
# Primary inference: seven spring-autumn station pairs (n = 14).
# Annual PCoA/heatmap: 19 libraries after excluding auxiliary N5/S5.

set.seed(20260819)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/")) else getwd()
root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
result_dir <- file.path(root, "results")
log_dir <- file.path(root, "logs")
maaslin_dir <- file.path(result_dir, "MaAsLin2_paired14")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(root, "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(vegan)
  library(Maaslin2)
})

required <- file.path(data_dir, c("dual_species_counts.tsv", "dual_genus_counts.tsv", "manifest_with_site.tsv"))
if (length(missing <- required[!file.exists(required)])) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

read_count_matrix <- function(path) {
  x <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  rownames(x) <- make.unique(as.character(x[[1]]))
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "double"
  m
}

species_all <- read_count_matrix(file.path(data_dir, "dual_species_counts.tsv"))
genus_all <- read_count_matrix(file.path(data_dir, "dual_genus_counts.tsv"))
manifest <- read.delim(file.path(data_dir, "manifest_with_site.tsv"), check.names = FALSE,
                       na.strings = c("", "NA"))
manifest <- manifest[manifest$lake == "ZBY" & manifest$matrix == "water", ]
manifest$sample_id <- as.character(manifest$sample_id)
manifest$site_corrected <- as.character(manifest$site)
manifest$site_corrected[manifest$sample_id == "Z1N5"] <- "N5"
manifest$site_corrected[manifest$sample_id == "Z1S5"] <- "S5"
manifest$site_corrected[manifest$sample_id %in% c("Z3N", "Z3S", "Z4N", "Z4S")] <- NA_character_
manifest$sampling_unit <- ifelse(manifest$sample_id %in% c("Z3N", "Z3S", "Z4N", "Z4S"),
                                 "basin_composite", "individual_station")
stopifnot(all(manifest$sample_id %in% colnames(species_all)), all(manifest$sample_id %in% colnames(genus_all)))

paired_sites <- c("N1", "N2", "N3", "N4", "S1", "S2", "S3")
ids14 <- manifest$sample_id[manifest$site_corrected %in% paired_sites & manifest$season %in% c("spring", "autumn")]
ids17 <- manifest$sample_id[manifest$sampling_unit == "individual_station" & manifest$season %in% c("spring", "autumn")]
ids19 <- manifest$sample_id[!manifest$sample_id %in% c("Z1N5", "Z1S5")]
stopifnot(length(ids14) == 14L, length(ids17) == 17L, length(ids19) == 19L)

make_meta <- function(ids) {
  m <- manifest[match(ids, manifest$sample_id), ]
  m$basin <- factor(m$basin, levels = c("north", "south"))
  m$season <- factor(m$season, levels = c("autumn", "spring"))
  m$site <- factor(m$site_corrected)
  site_n <- table(m$site)
  m$coverage_pattern <- factor(ifelse(site_n[as.character(m$site)] == 2L, "paired", "autumn_only"),
                               levels = c("paired", "autumn_only"))
  m
}
meta14 <- make_meta(ids14)
meta17 <- make_meta(ids17)

hat_matrix <- function(x) {
  qx <- qr(x); r <- qx$rank; q <- qr.Q(qx, complete = FALSE)
  if (r == 0L) return(matrix(0, nrow(x), nrow(x)))
  tcrossprod(q[, seq_len(r), drop = FALSE])
}
gower_matrix <- function(distance_matrix) {
  n <- nrow(distance_matrix); j <- diag(n) - matrix(1 / n, n, n)
  -0.5 * j %*% (distance_matrix ^ 2) %*% j
}
term_stat <- function(g, meta, term) {
  predictors <- c("season", "basin")
  reduced <- setdiff(predictors, term)
  xf <- model.matrix(as.formula(paste("~", paste(predictors, collapse = "+"))), meta)
  xr <- model.matrix(if (length(reduced)) as.formula(paste("~", paste(reduced, collapse = "+"))) else ~1, meta)
  hf <- hat_matrix(xf); hr <- hat_matrix(xr); n <- nrow(g); ident <- diag(n)
  ss_total <- sum(diag(g)); ss_term <- sum(diag((hf - hr) %*% g)); ss_resid <- sum(diag((ident - hf) %*% g))
  df_term <- qr(xf)$rank - qr(xr)$rank; df_resid <- n - qr(xf)$rank
  c(R2 = ss_term / ss_total, pseudo_F = (ss_term / df_term) / (ss_resid / df_resid),
    df_term = df_term, df_residual = df_resid)
}
season_permutations <- function(meta) {
  sites <- unique(as.character(meta$site[meta$coverage_pattern == "paired"]))
  combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(sites))))
  lapply(seq_len(nrow(combos)), function(i) {
    p <- meta; p$season <- as.character(p$season)
    for (j in seq_along(sites)) if (combos[i, j]) {
      idx <- which(as.character(p$site) == sites[j]); p$season[idx] <- rev(p$season[idx])
    }
    p$season <- factor(p$season, levels = c("autumn", "spring")); p
  })
}
basin_permutations <- function(meta) {
  sf <- unique(meta[, c("site", "basin", "coverage_pattern")])
  sf[] <- lapply(sf, as.character); patterns <- unique(sf$coverage_pattern)
  choices <- lapply(patterns, function(pattern) {
    z <- sf[sf$coverage_pattern == pattern, ]; combn(z$site, sum(z$basin == "north"), simplify = FALSE)
  })
  grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    north <- unlist(mapply(function(x, j) x[[j]], choices, as.integer(grid[i, ]), SIMPLIFY = FALSE))
    p <- meta; p$basin <- factor(ifelse(as.character(p$site) %in% north, "north", "south"),
                                 levels = c("north", "south")); p
  })
}
exact_test <- function(g, meta, term) {
  obs <- term_stat(g, meta, term)
  perms <- if (term == "season") season_permutations(meta) else basin_permutations(meta)
  f <- vapply(perms, function(p) term_stat(g, p, term)[["pseudo_F"]], numeric(1))
  c(obs, P_exact = sum(f >= obs[["pseudo_F"]] - sqrt(.Machine$double.eps)) / length(f),
    n_permutations = length(f))
}
bootstrap_indices <- function(meta) {
  sf <- unique(meta[, c("site", "basin", "coverage_pattern")]); sf[] <- lapply(sf, as.character)
  strata <- interaction(sf$basin, sf$coverage_pattern, drop = TRUE)
  sampled <- unlist(lapply(split(sf$site, strata), function(z) sample(z, length(z), replace = TRUE)), use.names = FALSE)
  idx <- integer(); boot_site <- character()
  for (i in seq_along(sampled)) {
    rows <- which(as.character(meta$site) == sampled[i]); idx <- c(idx, rows)
    boot_site <- c(boot_site, rep(paste0("boot_", i), length(rows)))
  }
  list(index = idx, site = boot_site)
}
bootstrap_r2 <- function(distance_matrix, meta, term, n_boot = 999L) {
  out <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    draw <- bootstrap_indices(meta); m <- meta[draw$index, ]; m$site <- factor(draw$site)
    d <- distance_matrix[draw$index, draw$index, drop = FALSE]
    out[b] <- tryCatch(term_stat(gower_matrix(d), m, term)[["R2"]], error = function(e) NA_real_)
  }
  out[is.finite(out)]
}
dispersion_test <- function(distance_matrix, meta, term) {
  get_f <- function(m) unname(anova(suppressWarnings(betadisper(as.dist(distance_matrix),
    if (term == "season") m$season else m$basin)))[1, "F value"])
  obs <- get_f(meta); perms <- if (term == "season") season_permutations(meta) else basin_permutations(meta)
  f <- vapply(perms, get_f, numeric(1))
  c(F = obs, P_exact = sum(f >= obs - sqrt(.Machine$double.eps)) / length(f), n_permutations = length(f))
}

# Descriptive annual composition for 19 libraries.
sp19 <- species_all[, ids19, drop = FALSE]
sp19 <- sp19[rowSums(sp19 > 0) >= 2, , drop = FALSE]
annual_species_ge2 <- nrow(sp19)
sp19_rel <- sweep(sp19, 2, colSums(sp19), "/")
bray19 <- vegdist(t(sp19_rel), method = "bray")
pcoa19 <- wcmdscale(bray19, k = 2, eig = TRUE)
eig <- pcoa19$eig[pcoa19$eig > 0]; pct <- 100 * eig[1:2] / sum(eig)
pcoa_out <- data.frame(sample_id = rownames(pcoa19$points), PCoA1 = pcoa19$points[, 1],
                       PCoA2 = pcoa19$points[, 2], PCoA1_percent = pct[1], PCoA2_percent = pct[2])
pcoa_out <- merge(pcoa_out, manifest[, c("sample_id", "season", "basin", "site_corrected", "sampling_unit")],
                  by = "sample_id", sort = FALSE)
write.table(pcoa_out, file.path(result_dir, "Figure_2b_dual_species_PCoA_19.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Genus heatmap data for the same descriptive set.
low_conf <- grepl("unclassified|uncultured|unknown|environmental sample|metagenome|^Homo$|^Mus$|^Rattus$",
                  rownames(genus_all), ignore.case = TRUE)
gen19 <- genus_all[!low_conf, ids19, drop = FALSE]
gen19_rel <- sweep(gen19, 2, colSums(gen19), "/")
top15 <- names(head(sort(rowMeans(gen19_rel), decreasing = TRUE), 15))
season_order <- c("spring", "summer", "autumn", "winter")
m19 <- manifest[match(ids19, manifest$sample_id), ]; m19$order <- match(m19$season, season_order) * 100 + match(m19$basin, c("north", "south")) * 10 + m19$replicate
ordered <- m19$sample_id[order(m19$order, m19$sample_id)]
heat <- do.call(rbind, lapply(top15, function(g) {
  x <- log10(gen19_rel[g, ordered] * 100 + min(gen19_rel[g, ordered][gen19_rel[g, ordered] > 0], na.rm = TRUE) * 50)
  z <- if (sd(x) > 0) as.numeric(scale(x)) else rep(0, length(x))
  data.frame(genus = g, sample_id = ordered, relative_abundance_pct = gen19_rel[g, ordered] * 100, z = z)
}))
heat <- merge(heat, manifest[, c("sample_id", "season", "basin", "site_corrected", "sampling_unit")], by = "sample_id", sort = FALSE)
heat$sample_order <- match(heat$sample_id, ordered); heat$genus_order <- match(heat$genus, top15)
write.table(heat, file.path(result_dir, "Figure_2a_dual_top15_genera.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(excluded_genus = rownames(genus_all)[low_conf]),
            file.path(result_dir, "excluded_low_confidence_genus_labels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Exact paired tests across primary and sensitivity occurrence filters.
primary_keep <- rowSums(species_all[, ids14, drop = FALSE] > 0)
rows <- list(); dispersion_rows <- list(); counter <- 0L
for (k in c(2L, 3L, 5L)) {
  features <- names(primary_keep)[primary_keep >= k]
  for (scope in c("paired14_primary", "individual17_sensitivity")) {
    ids <- if (scope == "paired14_primary") ids14 else ids17
    meta <- if (scope == "paired14_primary") meta14 else meta17
    counts <- species_all[features, ids, drop = FALSE]
    counts <- counts[apply(counts, 1, sd) > 0, , drop = FALSE]
    rel <- sweep(counts, 2, colSums(counts), "/")
    d <- as.matrix(vegdist(t(rel), method = "bray")); g <- gower_matrix(d)
    for (term in c("season", "basin")) {
      test <- exact_test(g, meta, term); boot <- bootstrap_r2(d, meta, term)
      counter <- counter + 1L
      rows[[counter]] <- data.frame(occurrence_filter = k, sample_scope = scope, n_samples = length(ids),
        n_sites = length(unique(meta$site)), n_species = nrow(counts), term = term,
        R2 = test[["R2"]], R2_CI_low = quantile(boot, 0.025), R2_CI_high = quantile(boot, 0.975),
        pseudo_F = test[["pseudo_F"]], P_exact = test[["P_exact"]], exact_permutations = test[["n_permutations"]])
      disp <- dispersion_test(d, meta, term)
      dispersion_rows[[counter]] <- data.frame(occurrence_filter = k, sample_scope = scope, term = term,
        dispersion_F = disp[["F"]], dispersion_P_exact = disp[["P_exact"]], permutations = disp[["n_permutations"]])
    }
  }
}
permanova <- do.call(rbind, rows); dispersion <- do.call(rbind, dispersion_rows)
write.table(permanova, file.path(result_dir, "dual_species_exact_PERMANOVA_n14_n17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dispersion, file.path(result_dir, "dual_species_dispersion_n14_n17.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Paired PCoA trajectories and leave-one-site-out stability.
features2 <- names(primary_keep)[primary_keep >= 2]
sp14 <- species_all[features2, ids14, drop = FALSE]; sp14 <- sweep(sp14, 2, colSums(sp14), "/")
pc14 <- wcmdscale(vegdist(t(sp14), method = "bray"), k = 2, eig = TRUE)
pc14_out <- cbind(meta14[, c("sample_id", "site", "season", "basin")], PCoA1 = pc14$points[, 1], PCoA2 = pc14$points[, 2])
write.table(pc14_out, file.path(result_dir, "Figure_S16a_dual_paired_trajectories.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
loso <- do.call(rbind, lapply(levels(meta14$site), function(drop_site) {
  keep <- meta14$site != drop_site; m <- droplevels(meta14[keep, ]); x <- sp14[, keep, drop = FALSE]
  d <- as.matrix(vegdist(t(x), method = "bray")); test <- exact_test(gower_matrix(d), m, "season")
  data.frame(omitted_site = drop_site, n = sum(keep), R2 = test[["R2"]], P_exact = test[["P_exact"]], permutations = test[["n_permutations"]])
}))
write.table(loso, file.path(result_dir, "Figure_S16b_dual_leave_one_site_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Genus-level paired differential abundance.
gen14 <- genus_all[!low_conf, ids14, drop = FALSE]
input_data <- as.data.frame(t(gen14), check.names = FALSE)
feature_original <- rownames(gen14)
feature_maaslin <- make.names(feature_original, unique = TRUE)
if (anyDuplicated(feature_maaslin)) stop("Non-unique MaAsLin2 feature identifiers after name sanitisation")
feature_lookup <- setNames(feature_original, feature_maaslin)
feature_lookup[feature_original] <- feature_original
metadata <- meta14[, c("sample_id", "season", "basin", "site")]
rownames(metadata) <- metadata$sample_id; metadata$sample_id <- NULL
metadata$season <- relevel(factor(metadata$season), ref = "autumn")
metadata$basin <- relevel(factor(metadata$basin), ref = "north")
if (dir.exists(maaslin_dir)) unlink(maaslin_dir, recursive = TRUE)
fit <- Maaslin2(input_data = input_data, input_metadata = metadata, output = maaslin_dir,
  fixed_effects = c("season", "basin"), random_effects = c("site"), normalization = "TSS",
  transform = "LOG", analysis_method = "LM", correction = "BH", standardize = FALSE,
  min_abundance = 0.0001, min_prevalence = 0.20, max_significance = 0.25,
  plot_heatmap = FALSE, plot_scatter = FALSE)
maaslin <- read.delim(file.path(maaslin_dir, "all_results.tsv"), check.names = FALSE)
mapped_feature <- unname(feature_lookup[as.character(maaslin$feature)])
if (anyNA(mapped_feature)) {
  stop("MaAsLin2 returned unmapped feature names: ",
       paste(unique(maaslin$feature[is.na(mapped_feature)]), collapse = ", "))
}
maaslin$feature <- mapped_feature
maaslin$conf_low <- maaslin$coef - 1.96 * maaslin$stderr
maaslin$conf_high <- maaslin$coef + 1.96 * maaslin$stderr
write.table(maaslin, file.path(result_dir, "dual_genus_MaAsLin2_all_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- maaslin[maaslin$metadata == "season" & maaslin$qval < 0.05, ]
sig <- sig[order(sig$qval, -abs(sig$coef)), ]
write.table(sig, file.path(result_dir, "dual_genus_MaAsLin2_season_q05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Abundance scale and paired-direction audit for any reported genus. This is
# descriptive support for the fitted coefficient, not an additional test.
gen14_rel <- sweep(gen14, 2, colSums(gen14), "/")
direction_audit <- do.call(rbind, lapply(as.character(sig$feature), function(feature) {
  values <- gen14_rel[feature, ]
  deltas <- vapply(levels(meta14$site), function(site_id) {
    idx <- which(meta14$site == site_id)
    values[idx[meta14$season[idx] == "spring"]] - values[idx[meta14$season[idx] == "autumn"]]
  }, numeric(1))
  data.frame(feature = feature, mean_relative_abundance = mean(values),
             median_relative_abundance = median(values),
             positive_site_pairs = sum(deltas > 0), negative_site_pairs = sum(deltas < 0),
             tied_site_pairs = sum(deltas == 0), n_site_pairs = length(deltas))
}))
if (is.null(direction_audit)) {
  direction_audit <- data.frame(feature = character(), mean_relative_abundance = numeric(),
                                median_relative_abundance = numeric(), positive_site_pairs = integer(),
                                negative_site_pairs = integer(), tied_site_pairs = integer(),
                                n_site_pairs = integer())
}
write.table(direction_audit, file.path(result_dir, "dual_genus_MaAsLin2_significant_abundance_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

membership <- manifest[, c("sample_id", "season", "basin", "site_corrected", "sampling_unit")]
membership$annual_descriptive_n19 <- membership$sample_id %in% ids19
membership$paired_primary_n14 <- membership$sample_id %in% ids14
membership$individual_sensitivity_n17 <- membership$sample_id %in% ids17
write.table(membership, file.path(result_dir, "dual_species_sample_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
denominators <- data.frame(
  metric = c("all_classified_libraries", "ZBY_annual_descriptive_libraries",
             "ZBY_paired_primary_libraries", "ZBY_individual_sensitivity_libraries",
             "merged_species_total", "merged_genera_total", "annual_species_occurrence_ge2",
             "paired_species_occurrence_ge2", "heatmap_genera"),
  value = c(ncol(species_all), length(ids19), length(ids14), length(ids17), nrow(species_all),
            nrow(genus_all), annual_species_ge2, length(features2), length(top15))
)
write.table(denominators, file.path(result_dir, "dual_species_analysis_denominators.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "R_sessionInfo_species_analysis.txt"))
print(permanova[permanova$occurrence_filter == 2, ])
cat("Significant season-associated genera at BH q < 0.05:", nrow(sig), "\n")
