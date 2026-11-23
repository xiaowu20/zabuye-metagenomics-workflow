# Reanalysis of ZBY basin and season effects using the corrected sampling units.
# Primary inference: 14 individually sequenced spring-autumn samples from seven
# recurrent stations. Sensitivity analysis: those 14 samples plus three
# autumn-only stations (N5, S4 and S5; n = 17).

set.seed(20260813)

suppressPackageStartupMessages({
  library(vegan)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript paired_models.R INPUTS OUTPUT")
processed_dir <- args[1]; results_dir <- file.path(args[2],"results"); logs_dir <- file.path(args[2],"logs")
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE); dir.create(logs_dir,showWarnings=FALSE)
fig2_input_dir <- processed_dir
mag_abundance_file <- file.path(processed_dir,"within_zby_relative_abundance.tsv")
required_files <- c(
  file.path(processed_dir, "ZBYinternal_manifest_with_site.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_sample_qc.tsv"),
  mag_abundance_file,
  file.path(processed_dir, "ZBYinternal_response_graded_weighted.tsv"),
  file.path(processed_dir, "ZBYinternal_response_strict_carrier_abundance.tsv"),
  file.path(processed_dir, "ZBYinternal_response_detected_carrier_fraction.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_graded_weighted.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_strict_carrier_abundance.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_detected_carrier_fraction.tsv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required files:\n", paste(missing_files, collapse = "\n"))

read_response <- function(path) {
  x <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  names(x)[1] <- "sample_id"
  x$sample_id <- as.character(x$sample_id)
  value_cols <- setdiff(names(x), "sample_id")
  x[value_cols] <- lapply(x[value_cols], as.numeric)
  x
}

manifest <- read.delim(
  file.path(processed_dir, "ZBYinternal_manifest_with_site.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
manifest <- manifest[manifest$lake == "ZBY" & manifest$matrix == "water", ]
manifest$sample_id <- as.character(manifest$sample_id)

# Correct the old manifest: N5 and S5 are genuine autumn-only stations, while
# summer and winter libraries are basin composites and are not station samples.
manifest$site_corrected <- as.character(manifest$site)
manifest$site_corrected[manifest$sample_id == "Z1N5"] <- "N5"
manifest$site_corrected[manifest$sample_id == "Z1S5"] <- "S5"
manifest$site_corrected[manifest$sample_id %in% c("Z3N", "Z3S", "Z4N", "Z4S")] <- NA_character_
manifest$sampling_unit <- ifelse(
  manifest$sample_id %in% c("Z3N", "Z3S", "Z4N", "Z4S"),
  "basin_composite", "individual_station"
)

paired_sites <- c("N1", "N2", "N3", "N4", "S1", "S2", "S3")
ids14 <- manifest$sample_id[
  manifest$site_corrected %in% paired_sites & manifest$season %in% c("spring", "autumn")
]
ids17 <- manifest$sample_id[
  manifest$sampling_unit == "individual_station" & manifest$season %in% c("spring", "autumn")
]
stopifnot(length(ids14) == 14L, length(ids17) == 17L)

make_meta <- function(ids) {
  m <- manifest[match(ids, manifest$sample_id), ]
  stopifnot(all(m$sample_id == ids), !anyNA(m$site_corrected))
  m$basin <- factor(m$basin, levels = c("north", "south"))
  m$season <- factor(m$season, levels = c("autumn", "spring"))
  m$site <- factor(m$site_corrected)
  counts <- table(m$site)
  m$coverage_pattern <- ifelse(counts[as.character(m$site)] == 2L, "paired", "autumn_only")
  m$coverage_pattern <- factor(m$coverage_pattern, levels = c("paired", "autumn_only"))
  m
}

meta14 <- make_meta(ids14)
meta17 <- make_meta(ids17)
stopifnot(
  all(table(meta14$site) == 2L),
  sum(meta17$coverage_pattern == "paired") == 14L,
  sum(meta17$coverage_pattern == "autumn_only") == 3L
)

qc414 <- read.delim(
  file.path(processed_dir, "zbyonly21_cat414_sample_qc.tsv"),
  check.names = FALSE
)
qc414$sample_id <- as.character(qc414$sample_id)

mag_raw <- read.delim(mag_abundance_file, check.names = FALSE, quote = "", comment.char = "")
rownames(mag_raw) <- make.unique(as.character(mag_raw[[1]]))
mag_raw[[1]] <- NULL
mag_matrix <- t(as.matrix(mag_raw))
storage.mode(mag_matrix) <- "double"
mag_detected_269 <- rowSums(mag_matrix > 0)

response_list <- list()
response_list$mag_relative_abundance_269 <- list(
  matrix = mag_matrix,
  catalogue = "fixed 269 ZBY-water MAG catalogue",
  response_class = "primary",
  n_features = ncol(mag_matrix),
  tech = mag_detected_269
)

add_tabular_response <- function(name, file_name, catalogue, response_class, tech_vector) {
  x <- read_response(file.path(processed_dir, file_name))
  mat <- as.matrix(x[, setdiff(names(x), "sample_id"), drop = FALSE])
  rownames(mat) <- x$sample_id
  storage.mode(mat) <- "double"
  response_list[[name]] <<- list(
    matrix = mat,
    catalogue = catalogue,
    response_class = response_class,
    n_features = ncol(mat),
    tech = tech_vector
  )
}

add_tabular_response(
  "graded_weighted_modules_269", "ZBYinternal_response_graded_weighted.tsv",
  "fixed 269 ZBY-water MAG catalogue", "primary", mag_detected_269
)
add_tabular_response(
  "strict_carrier_abundance_269", "ZBYinternal_response_strict_carrier_abundance.tsv",
  "fixed 269 ZBY-water MAG catalogue", "primary", mag_detected_269
)
add_tabular_response(
  "detected_carrier_fraction_269", "ZBYinternal_response_detected_carrier_fraction.tsv",
  "fixed 269 ZBY-water MAG catalogue", "primary", mag_detected_269
)

tech414 <- setNames(qc414$n_detected_mags, qc414$sample_id)
add_tabular_response(
  "graded_weighted_modules_414", "zbyonly21_cat414_response_graded_weighted.tsv",
  "joint 414-MAG catalogue", "catalogue_sensitivity", tech414
)
add_tabular_response(
  "strict_carrier_abundance_414", "zbyonly21_cat414_response_strict_carrier_abundance.tsv",
  "joint 414-MAG catalogue", "catalogue_sensitivity", tech414
)
add_tabular_response(
  "detected_carrier_fraction_414", "zbyonly21_cat414_response_detected_carrier_fraction.tsv",
  "joint 414-MAG catalogue", "catalogue_sensitivity", tech414
)

hat_matrix <- function(x) {
  qx <- qr(x)
  r <- qx$rank
  q <- qr.Q(qx, complete = FALSE)
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

term_stat_from_g <- function(g, meta, term, use_tech) {
  predictors <- c(if (use_tech) "n_detected_mags", "season", "basin")
  reduced_predictors <- setdiff(predictors, term)
  f_full <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  f_reduced <- if (length(reduced_predictors)) {
    as.formula(paste("~", paste(reduced_predictors, collapse = " + ")))
  } else {
    ~1
  }
  x_full <- model.matrix(f_full, data = meta)
  x_reduced <- model.matrix(f_reduced, data = meta)
  h_full <- hat_matrix(x_full)
  h_reduced <- hat_matrix(x_reduced)
  n <- nrow(g)
  i <- diag(n)
  ss_total <- sum(diag(g))
  ss_term <- sum(diag((h_full - h_reduced) %*% g))
  ss_residual <- sum(diag((i - h_full) %*% g))
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
  site_levels <- unique(as.character(meta$site[meta$coverage_pattern == "paired"]))
  combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(site_levels))))
  out <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    p <- meta
    p$season <- as.character(p$season)
    for (j in seq_along(site_levels)) {
      if (combos[i, j]) {
        idx <- which(as.character(p$site) == site_levels[j])
        p$season[idx] <- rev(p$season[idx])
      }
    }
    p$season <- factor(p$season, levels = c("autumn", "spring"))
    out[[i]] <- p
  }
  out
}

basin_label_permutations <- function(meta) {
  site_frame <- unique(meta[, c("site", "basin", "coverage_pattern")])
  site_frame$site <- as.character(site_frame$site)
  site_frame$basin <- as.character(site_frame$basin)
  site_frame$coverage_pattern <- as.character(site_frame$coverage_pattern)
  patterns <- unique(site_frame$coverage_pattern)
  choice_lists <- lapply(patterns, function(pattern) {
    z <- site_frame[site_frame$coverage_pattern == pattern, ]
    n_north <- sum(z$basin == "north")
    combn(z$site, n_north, simplify = FALSE)
  })
  grid <- expand.grid(lapply(choice_lists, seq_along), KEEP.OUT.ATTRS = FALSE)
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    north_sites <- unlist(mapply(
      function(choices, idx) choices[[idx]], choice_lists,
      as.integer(grid[i, ]), SIMPLIFY = FALSE
    ))
    p <- meta
    p$basin <- factor(
      ifelse(as.character(p$site) %in% north_sites, "north", "south"),
      levels = c("north", "south")
    )
    out[[i]] <- p
  }
  out
}

exact_term_test <- function(g, meta, term, use_tech) {
  observed <- term_stat_from_g(g, meta, term, use_tech)
  permutations <- if (term == "season") {
    season_label_permutations(meta)
  } else {
    basin_label_permutations(meta)
  }
  perm_f <- vapply(
    permutations,
    function(p) term_stat_from_g(g, p, term, use_tech)[["pseudo_F"]],
    numeric(1)
  )
  p_exact <- sum(perm_f >= observed[["pseudo_F"]] - sqrt(.Machine$double.eps)) / length(perm_f)
  list(observed = observed, p_exact = p_exact, n_permutations = length(perm_f), perm_f = perm_f)
}

bootstrap_indices <- function(meta) {
  site_frame <- unique(meta[, c("site", "basin", "coverage_pattern")])
  site_frame$site <- as.character(site_frame$site)
  strata <- interaction(site_frame$basin, site_frame$coverage_pattern, drop = TRUE)
  sampled_sites <- unlist(lapply(split(site_frame$site, strata), function(z) {
    sample(z, length(z), replace = TRUE)
  }), use.names = FALSE)
  idx <- integer(0)
  boot_site <- character(0)
  counter <- 0L
  for (s in sampled_sites) {
    counter <- counter + 1L
    rows <- which(as.character(meta$site) == s)
    idx <- c(idx, rows)
    boot_site <- c(boot_site, rep(paste0("boot_site_", counter), length(rows)))
  }
  list(index = idx, site = boot_site)
}

bootstrap_term <- function(distance_matrix, meta, term, use_tech, n_boot = 999L) {
  r2 <- rep(NA_real_, n_boot)
  adjusted_r2_values <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    draw <- bootstrap_indices(meta)
    m <- meta[draw$index, , drop = FALSE]
    m$site <- factor(draw$site)
    d <- distance_matrix[draw$index, draw$index, drop = FALSE]
    g <- gower_matrix(d)
    stat <- tryCatch(term_stat_from_g(g, m, term, use_tech), error = function(e) NULL)
    if (!is.null(stat)) {
      r2[b] <- stat[["R2"]]
      adjusted_r2_values[b] <- stat[["adjusted_R2"]]
    }
  }
  list(
    R2 = r2[is.finite(r2)],
    adjusted_R2 = adjusted_r2_values[is.finite(adjusted_r2_values)]
  )
}

dispersion_exact_test <- function(distance_matrix, meta, term) {
  get_f <- function(m) {
    group <- if (term == "season") m$season else m$basin
    fit <- suppressWarnings(betadisper(as.dist(distance_matrix), group))
    unname(anova(fit)[1, "F value"])
  }
  observed_f <- get_f(meta)
  permutations <- if (term == "season") {
    season_label_permutations(meta)
  } else {
    basin_label_permutations(meta)
  }
  perm_f <- vapply(permutations, get_f, numeric(1))
  data.frame(
    term = term,
    dispersion_F = observed_f,
    dispersion_P_exact = sum(perm_f >= observed_f - sqrt(.Machine$double.eps)) / length(perm_f),
    n_permutations = length(perm_f)
  )
}

validate_against_adonis <- function(mat, meta, term, use_tech) {
  dat <- meta
  if (use_tech) {
    fit <- adonis2(
      mat ~ n_detected_mags + season + basin,
      data = dat, method = "bray", by = "margin", permutations = 9
    )
  } else {
    fit <- adonis2(
      mat ~ season + basin,
      data = dat, method = "bray", by = "margin", permutations = 9
    )
  }
  c(
    vegan_R2 = unname(fit[term, "R2"]),
    vegan_F = unname(fit[term, "F"])
  )
}

result_rows <- list()
validation_rows <- list()
dispersion_rows <- list()
row_counter <- 0L
validation_counter <- 0L
dispersion_counter <- 0L

for (response_name in names(response_list)) {
  response <- response_list[[response_name]]
  for (scope_name in c("paired14_primary", "individual17_sensitivity")) {
    meta <- if (scope_name == "paired14_primary") meta14 else meta17
    ids <- meta$sample_id
    stopifnot(all(ids %in% rownames(response$matrix)))
    mat <- response$matrix[ids, , drop = FALSE]
    keep_columns <- apply(mat, 2, function(z) any(is.finite(z)) && stats::sd(z, na.rm = TRUE) > 0)
    mat <- mat[, keep_columns, drop = FALSE]
    if (!ncol(mat)) stop("No variable response columns for ", response_name, " / ", scope_name)
    distance_matrix <- as.matrix(vegdist(mat, method = "bray"))
    g <- gower_matrix(distance_matrix)

    if (!is.null(response$tech)) {
      tech <- response$tech[ids]
      if (anyNA(tech)) stop("Missing technical covariate for ", response_name, " / ", scope_name)
      meta$n_detected_mags <- as.numeric(scale(as.numeric(tech)))
    }

    adjustment_states <- if (is.null(response$tech)) FALSE else c(FALSE, TRUE)
    for (use_tech in adjustment_states) {
      for (term in c("season", "basin")) {
        test <- exact_term_test(g, meta, term, use_tech)
        n_boot <- if (response$response_class == "filter_sensitivity") 199L else 399L
        boot <- bootstrap_term(distance_matrix, meta, term, use_tech, n_boot = n_boot)
        stat <- test$observed
        row_counter <- row_counter + 1L
        result_rows[[row_counter]] <- data.frame(
          response = response_name,
          catalogue = response$catalogue,
          response_class = response$response_class,
          sample_scope = scope_name,
          n_samples = nrow(meta),
          n_sites = length(unique(meta$site)),
          n_features_input = response$n_features,
          n_features_variable = ncol(mat),
          technical_adjustment = if (use_tech) "n_detected_mags" else "none",
          term = term,
          R2 = stat[["R2"]],
          R2_CI_low = unname(quantile(boot$R2, 0.025, na.rm = TRUE)),
          R2_CI_high = unname(quantile(boot$R2, 0.975, na.rm = TRUE)),
          adjusted_R2 = stat[["adjusted_R2"]],
          adjusted_R2_CI_low = unname(quantile(boot$adjusted_R2, 0.025, na.rm = TRUE)),
          adjusted_R2_CI_high = unname(quantile(boot$adjusted_R2, 0.975, na.rm = TRUE)),
          pseudo_F = stat[["pseudo_F"]],
          df_term = stat[["df_term"]],
          df_residual = stat[["df_residual"]],
          P_exact = test$p_exact,
          exact_permutations = test$n_permutations,
          bootstrap_replicates = length(boot$R2),
          permutation_unit = if (term == "season") {
            "within-site spring-autumn label swaps; autumn-only sites fixed"
          } else {
            "whole-site basin labels; paired/autumn-only coverage pattern preserved"
          },
          stringsAsFactors = FALSE
        )

        vegan_stat <- validate_against_adonis(mat, meta, term, use_tech)
        validation_counter <- validation_counter + 1L
        validation_rows[[validation_counter]] <- data.frame(
          response = response_name,
          sample_scope = scope_name,
          technical_adjustment = if (use_tech) "n_detected_mags" else "none",
          term = term,
          custom_R2 = stat[["R2"]],
          vegan_R2 = vegan_stat[["vegan_R2"]],
          R2_absolute_difference = abs(stat[["R2"]] - vegan_stat[["vegan_R2"]]),
          custom_F = stat[["pseudo_F"]],
          vegan_F = vegan_stat[["vegan_F"]],
          F_absolute_difference = abs(stat[["pseudo_F"]] - vegan_stat[["vegan_F"]]),
          stringsAsFactors = FALSE
        )
      }
    }

    for (term in c("season", "basin")) {
      dispersion_counter <- dispersion_counter + 1L
      drow <- dispersion_exact_test(distance_matrix, meta, term)
      dispersion_rows[[dispersion_counter]] <- cbind(
        data.frame(
          response = response_name,
          catalogue = response$catalogue,
          response_class = response$response_class,
          sample_scope = scope_name,
          n_samples = nrow(meta),
          n_sites = length(unique(meta$site)),
          stringsAsFactors = FALSE
        ),
        drow
      )
    }
  }
}

results <- do.call(rbind, result_rows)
validation <- do.call(rbind, validation_rows)
dispersion <- do.call(rbind, dispersion_rows)

write.table(
  results,
  file.path(results_dir, "permanova_n14_n17_all_responses.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  validation,
  file.path(results_dir, "custom_statistics_validation_against_vegan.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  dispersion,
  file.path(results_dir, "dispersion_n14_n17_all_responses.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

membership <- manifest[, c(
  "sample_id", "basin", "season", "site_corrected", "sampling_unit"
)]
membership$included_n14_primary <- membership$sample_id %in% ids14
membership$included_n17_sensitivity <- membership$sample_id %in% ids17
membership$analysis_role <- ifelse(
  membership$included_n14_primary, "paired_station_primary",
  ifelse(membership$included_n17_sensitivity, "autumn_only_sensitivity",
         ifelse(membership$sampling_unit == "basin_composite", "annual_descriptive_only", "excluded"))
)
write.table(
  membership,
  file.path(results_dir, "sample_membership_corrected.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_columns <- c(
  "response", "catalogue", "sample_scope", "technical_adjustment", "term",
  "R2", "adjusted_R2", "P_exact", "exact_permutations"
)
print(results[results$response_class != "filter_sensitivity", summary_columns], row.names = FALSE)

capture.output(sessionInfo(), file = file.path(logs_dir, "sessionInfo.txt"))
