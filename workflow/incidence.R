# Reanalysis of the post hoc binary module-incidence layer using corrected
# spring-autumn station units. This script intentionally omits the invalid
# historical n = 19 set.

set.seed(20260813)
suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript incidence.R INPUTS OUTPUT")
processed_dir <- args[1]; results_dir <- file.path(args[2],"results"); logs_dir <- file.path(args[2],"logs")
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE); dir.create(logs_dir,showWarnings=FALSE)
read_response <- function(path) {
  x <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  names(x)[1] <- "sample_id"
  x$sample_id <- as.character(x$sample_id)
  x
}

manifest <- read.delim(
  file.path(processed_dir, "ZBYinternal_manifest_with_site.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
manifest <- manifest[manifest$lake == "ZBY" & manifest$matrix == "water", ]
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

make_meta <- function(ids) {
  m <- manifest[match(ids, manifest$sample_id), ]
  m$basin <- factor(m$basin, levels = c("north", "south"))
  m$season <- factor(m$season, levels = c("autumn", "spring"))
  m$site <- factor(m$site_corrected)
  counts <- table(m$site)
  m$coverage_pattern <- factor(
    ifelse(counts[as.character(m$site)] == 2L, "paired", "autumn_only"),
    levels = c("paired", "autumn_only")
  )
  m
}
meta14 <- make_meta(ids14)
meta17 <- make_meta(ids17)

hat_matrix <- function(x) {
  qx <- qr(x)
  q <- qr.Q(qx, complete = FALSE)
  tcrossprod(q[, seq_len(qx$rank), drop = FALSE])
}

gower_matrix <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  j <- diag(n) - matrix(1 / n, n, n)
  -0.5 * j %*% (distance_matrix ^ 2) %*% j
}

term_stat <- function(g, meta, term) {
  full <- model.matrix(~ season + basin, meta)
  reduced <- if (term == "season") model.matrix(~ basin, meta) else model.matrix(~ season, meta)
  h_full <- hat_matrix(full)
  h_reduced <- hat_matrix(reduced)
  ss_total <- sum(diag(g))
  ss_term <- sum(diag((h_full - h_reduced) %*% g))
  ss_res <- sum(diag((diag(nrow(g)) - h_full) %*% g))
  df_term <- qr(full)$rank - qr(reduced)$rank
  df_res <- nrow(g) - qr(full)$rank
  c(R2 = ss_term / ss_total, pseudo_F = (ss_term / df_term) / (ss_res / df_res))
}

season_permutations <- function(meta) {
  paired <- unique(as.character(meta$site[meta$coverage_pattern == "paired"]))
  combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(paired))))
  lapply(seq_len(nrow(combos)), function(i) {
    p <- meta
    p$season <- as.character(p$season)
    for (j in seq_along(paired)) {
      if (combos[i, j]) {
        idx <- which(as.character(p$site) == paired[j])
        p$season[idx] <- rev(p$season[idx])
      }
    }
    p$season <- factor(p$season, levels = c("autumn", "spring"))
    p
  })
}

basin_permutations <- function(meta) {
  sites <- unique(meta[, c("site", "basin", "coverage_pattern")])
  sites$site <- as.character(sites$site)
  sites$basin <- as.character(sites$basin)
  sites$coverage_pattern <- as.character(sites$coverage_pattern)
  patterns <- unique(sites$coverage_pattern)
  choices <- lapply(patterns, function(pattern) {
    z <- sites[sites$coverage_pattern == pattern, ]
    combn(z$site, sum(z$basin == "north"), simplify = FALSE)
  })
  grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    north <- unlist(mapply(
      function(z, idx) z[[idx]], choices, as.integer(grid[i, ]), SIMPLIFY = FALSE
    ))
    p <- meta
    p$basin <- factor(
      ifelse(as.character(p$site) %in% north, "north", "south"),
      levels = c("north", "south")
    )
    p
  })
}

exact_test <- function(g, meta, term) {
  observed <- term_stat(g, meta, term)
  permutations <- if (term == "season") season_permutations(meta) else basin_permutations(meta)
  perm_f <- vapply(permutations, function(p) term_stat(g, p, term)[["pseudo_F"]], numeric(1))
  c(
    R2 = observed[["R2"]],
    pseudo_F = observed[["pseudo_F"]],
    P_exact = sum(perm_f >= observed[["pseudo_F"]] - sqrt(.Machine$double.eps)) / length(perm_f),
    exact_permutations = length(perm_f)
  )
}

dispersion_test <- function(distance_matrix, meta, term) {
  get_f <- function(p) {
    group <- if (term == "season") p$season else p$basin
    unname(anova(suppressWarnings(betadisper(as.dist(distance_matrix), group)))[1, "F value"])
  }
  observed <- get_f(meta)
  permutations <- if (term == "season") season_permutations(meta) else basin_permutations(meta)
  perm_f <- vapply(permutations, get_f, numeric(1))
  c(
    dispersion_F = observed,
    dispersion_P_exact = sum(perm_f >= observed - sqrt(.Machine$double.eps)) / length(perm_f)
  )
}

carrier_files <- c(
  "269" = "ZBYinternal_response_detected_carrier_fraction.tsv",
  "414" = "zbyonly21_cat414_response_detected_carrier_fraction.tsv"
)
thresholds <- c(gt0 = 0, ge005 = 0.05, ge010 = 0.10)
rows <- list()
overlap_rows <- list()
counter <- 0L
overlap_counter <- 0L

for (catalogue in names(carrier_files)) {
  x <- read_response(file.path(processed_dir, carrier_files[[catalogue]]))
  values <- as.matrix(x[, setdiff(names(x), "sample_id"), drop = FALSE])
  rownames(values) <- x$sample_id
  storage.mode(values) <- "double"
  for (scope in c("paired14_primary", "individual17_sensitivity")) {
    meta <- if (scope == "paired14_primary") meta14 else meta17
    for (threshold_name in names(thresholds)) {
      threshold <- thresholds[[threshold_name]]
      binary <- if (threshold == 0) values[meta$sample_id, , drop = FALSE] > 0 else
        values[meta$sample_id, , drop = FALSE] >= threshold
      variable <- apply(binary, 2, function(z) length(unique(z)) > 1)
      binary_for_test <- binary[, variable, drop = FALSE]
      if (!ncol(binary_for_test)) next
      distance_matrix <- as.matrix(vegdist(binary_for_test * 1, method = "jaccard", binary = TRUE))
      g <- gower_matrix(distance_matrix)
      for (term in c("season", "basin")) {
        test <- exact_test(g, meta, term)
        dispersion <- dispersion_test(distance_matrix, meta, term)
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          catalogue = paste0(catalogue, "-MAG"),
          sample_scope = scope,
          n_samples = nrow(meta),
          n_sites = length(unique(meta$site)),
          threshold_name = threshold_name,
          carrier_fraction_threshold = threshold,
          n_variable_modules = ncol(binary_for_test),
          term = term,
          R2 = test[["R2"]],
          pseudo_F = test[["pseudo_F"]],
          P_exact = test[["P_exact"]],
          exact_permutations = test[["exact_permutations"]],
          dispersion_F = dispersion[["dispersion_F"]],
          dispersion_P_exact = dispersion[["dispersion_P_exact"]],
          stringsAsFactors = FALSE
        )
      }

      north_pool <- colSums(binary[meta$basin == "north", , drop = FALSE]) > 0
      south_pool <- colSums(binary[meta$basin == "south", , drop = FALSE]) > 0
      shared <- sum(north_pool & south_pool)
      union <- sum(north_pool | south_pool)
      overlap_counter <- overlap_counter + 1L
      overlap_rows[[overlap_counter]] <- data.frame(
        catalogue = paste0(catalogue, "-MAG"),
        sample_scope = scope,
        threshold_name = threshold_name,
        carrier_fraction_threshold = threshold,
        north_pool_modules = sum(north_pool),
        south_pool_modules = sum(south_pool),
        shared_modules = shared,
        union_modules = union,
        jaccard_similarity = if (union > 0) shared / union else NA_real_,
        north_only_modules = sum(north_pool & !south_pool),
        south_only_modules = sum(!north_pool & south_pool),
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, rows)
overlap <- do.call(rbind, overlap_rows)
write.table(
  results, file.path(results_dir, "module_incidence_n14_n17_exact.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  overlap, file.path(results_dir, "module_incidence_pooled_overlap_n14_n17.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(results, row.names = FALSE)
print(overlap, row.names = FALSE)
