# Tests whether site-level salinity and major-ion composition account for
# basin-associated genome-resolved responses in the paired spring-autumn set.

set.seed(20260814)

suppressPackageStartupMessages({
  library(readxl)
  library(vegan)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript ion_models.R INPUTS_DIR NEW_OUTPUT_DIR")
processed_dir <- normalizePath(args[1])
if(dir.exists(args[2])) stop("Choose a new analysis output directory")
analysis_root <- args[2]
results_dir <- file.path(analysis_root, "results")
data_dir <- file.path(analysis_root, "data")
figures_dir <- file.path(analysis_root, "figures")
logs_dir <- file.path(analysis_root, "logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

mag_abundance_file <- file.path(processed_dir, "within_zby_relative_abundance.tsv")
environment_file <- file.path(processed_dir, "site_chemistry.xlsx")
hydrology_file <- file.path(processed_dir, "hydrology.xlsx")
required_files <- c(
  file.path(processed_dir, "sample_membership_corrected.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_sample_qc.tsv"),
  mag_abundance_file,
  file.path(processed_dir, "ZBYinternal_response_graded_weighted.tsv"),
  file.path(processed_dir, "ZBYinternal_response_strict_carrier_abundance.tsv"),
  file.path(processed_dir, "ZBYinternal_response_detected_carrier_fraction.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_graded_weighted.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_strict_carrier_abundance.tsv"),
  file.path(processed_dir, "zbyonly21_cat414_response_detected_carrier_fraction.tsv"),
  environment_file,
  hydrology_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

read_response <- function(path) {
  x <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  names(x)[1] <- "sample_id"
  x$sample_id <- as.character(x$sample_id)
  value_cols <- setdiff(names(x), "sample_id")
  x[value_cols] <- lapply(x[value_cols], as.numeric)
  mat <- as.matrix(x[, value_cols, drop = FALSE])
  rownames(mat) <- x$sample_id
  storage.mode(mat) <- "double"
  mat
}

membership <- read.delim(
  file.path(processed_dir, "sample_membership_corrected.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
membership <- membership[membership$included_n14_primary, ]
stopifnot(nrow(membership) == 14L)

membership$site <- membership$site_corrected
membership$basin <- factor(membership$basin, levels = c("north", "south"))
membership$season <- factor(membership$season, levels = c("autumn", "spring"))

# Map metagenome IDs to the site-season identifiers in the hydrochemistry file.
season_suffix <- c(autumn = "Aut", spring = "Spr")
membership$environment_id <- paste0(
  membership$site, "_", season_suffix[as.character(membership$season)]
)

env_raw <- as.data.frame(read_excel(environment_file, sheet = "Dataset S1 complete"))
env <- env_raw[match(membership$environment_id, env_raw$`Sample ID`), ]
if (anyNA(env$`Sample ID`)) {
  stop("Unmatched environmental IDs: ", paste(membership$environment_id[is.na(env$`Sample ID`)], collapse = ", "))
}

meta <- membership[, c("sample_id", "site", "basin", "season", "environment_id")]
meta$salinity_g_L <- as.numeric(env$`Total salt (g L^-1)`)
meta$pH <- as.numeric(env$pH)

# Convert consistently positive major ions to mol L^-1 before CLR analysis.
# Nutrient ions and carbonate/bicarbonate are excluded because they include
# zeros or mixed chemical species and do not represent the dominant salt matrix.
ion_definition <- data.frame(
  ion = c("Na", "K", "Li", "Mg", "Ca", "Cl", "Br", "F", "SO4"),
  source_column = c(
    "Na+ (g L^-1)", "K+ (g L^-1)", "Li+ (g L^-1)",
    "Mg2+ (mg L^-1)", "Ca2+ (mg L^-1)", "Cl- (g L^-1)",
    "Br- (g L^-1)", "F- (mg L^-1)", "SO4^2- (g L^-1)"
  ),
  grams_per_mole = c(22.989769, 39.0983, 6.94, 24.305, 40.078, 35.45, 79.904, 18.998403, 96.06),
  source_unit = c("g/L", "g/L", "g/L", "mg/L", "mg/L", "g/L", "g/L", "mg/L", "g/L"),
  stringsAsFactors = FALSE
)

ion_molar <- sapply(seq_len(nrow(ion_definition)), function(i) {
  values <- as.numeric(env[[ion_definition$source_column[i]]])
  if (ion_definition$source_unit[i] == "mg/L") values <- values / 1000
  values / ion_definition$grams_per_mole[i]
})
colnames(ion_molar) <- ion_definition$ion
rownames(ion_molar) <- meta$sample_id
if (any(!is.finite(ion_molar)) || any(ion_molar <= 0)) {
  stop("Major-ion matrix contains non-positive or non-finite values; CLR is not valid.")
}

clr <- log(ion_molar) - rowMeans(log(ion_molar))
ion_pca <- prcomp(clr, center = TRUE, scale. = FALSE)
meta$ion_PC1 <- ion_pca$x[, 1]
meta$ion_PC2 <- ion_pca$x[, 2]
meta$salinity_z <- as.numeric(scale(meta$salinity_g_L))
meta$ion_PC1_z <- as.numeric(scale(meta$ion_PC1))
meta$ion_PC2_z <- as.numeric(scale(meta$ion_PC2))

write.table(meta, file.path(data_dir, "environment_n14.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  cbind(sample_id = rownames(ion_molar), as.data.frame(ion_molar)),
  file.path(data_dir, "major_ions_molar_n14.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cbind(sample_id = rownames(clr), as.data.frame(clr)),
  file.path(data_dir, "major_ions_CLR_n14.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
loadings <- data.frame(
  ion = rownames(ion_pca$rotation),
  PC1 = ion_pca$rotation[, 1],
  PC2 = ion_pca$rotation[, 2],
  stringsAsFactors = FALSE
)
write.table(loadings, file.path(results_dir, "ion_PCA_loadings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(
    axis = paste0("PC", seq_along(ion_pca$sdev)),
    variance_fraction = ion_pca$sdev^2 / sum(ion_pca$sdev^2)
  ),
  file.path(results_dir, "ion_PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(ion_definition, file.path(data_dir, "major_ion_definition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

group_summary <- aggregate(
  cbind(salinity_g_L, ion_PC1, ion_PC2) ~ basin + season,
  data = meta,
  FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x))
)
group_summary_flat <- data.frame(
  basin = group_summary$basin,
  season = group_summary$season,
  salinity_mean = group_summary$salinity_g_L[, "mean"],
  salinity_sd = group_summary$salinity_g_L[, "sd"],
  ion_PC1_mean = group_summary$ion_PC1[, "mean"],
  ion_PC1_sd = group_summary$ion_PC1[, "sd"],
  ion_PC2_mean = group_summary$ion_PC2[, "mean"],
  ion_PC2_sd = group_summary$ion_PC2[, "sd"],
  n = group_summary$salinity_g_L[, "n"]
)
write.table(
  group_summary_flat, file.path(results_dir, "environment_basin_season_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Hydrological metrics are retained as four basin-season contextual units.
hydro <- as.data.frame(read_excel(hydrology_file, sheet = "basin_year_season_metrics"))
hydro <- hydro[hydro$year == 2023 & hydro$season %in% c("Spring", "Autumn"), ]
hydro <- hydro[, c(
  "year", "season", "basin", "open_water_area_km2", "wet_playa_fraction",
  "active_surface_area_km2", "snow_ice_fraction", "mndwi_mean", "ndwi_mean", "scene_key"
)]
write.table(hydro, file.path(data_dir, "hydrology_2023_basin_season_context.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

qc414 <- read.delim(file.path(processed_dir, "zbyonly21_cat414_sample_qc.tsv"), check.names = FALSE)
tech414 <- setNames(qc414$n_detected_mags, qc414$sample_id)

mag_raw <- read.delim(mag_abundance_file, check.names = FALSE, quote = "", comment.char = "")
rownames(mag_raw) <- make.unique(as.character(mag_raw[[1]]))
mag_raw[[1]] <- NULL
mag_matrix <- t(as.matrix(mag_raw))
storage.mode(mag_matrix) <- "double"
tech269 <- rowSums(mag_matrix > 0)

responses <- list(
  MAG_relative_abundance_269 = list(matrix = mag_matrix, tech = tech269, class = "MAG composition"),
  graded_modules_269 = list(
    matrix = read_response(file.path(processed_dir, "ZBYinternal_response_graded_weighted.tsv")),
    tech = tech269, class = "abundance-weighted modules"
  ),
  strict_modules_269 = list(
    matrix = read_response(file.path(processed_dir, "ZBYinternal_response_strict_carrier_abundance.tsv")),
    tech = tech269, class = "abundance-weighted modules"
  ),
  carrier_fraction_269 = list(
    matrix = read_response(file.path(processed_dir, "ZBYinternal_response_detected_carrier_fraction.tsv")),
    tech = tech269, class = "carrier prevalence"
  ),
  graded_modules_414 = list(
    matrix = read_response(file.path(processed_dir, "zbyonly21_cat414_response_graded_weighted.tsv")),
    tech = tech414, class = "catalogue sensitivity"
  ),
  strict_modules_414 = list(
    matrix = read_response(file.path(processed_dir, "zbyonly21_cat414_response_strict_carrier_abundance.tsv")),
    tech = tech414, class = "catalogue sensitivity"
  ),
  carrier_fraction_414 = list(
    matrix = read_response(file.path(processed_dir, "zbyonly21_cat414_response_detected_carrier_fraction.tsv")),
    tech = tech414, class = "catalogue sensitivity"
  )
)

hat_matrix <- function(x) {
  qx <- qr(x)
  if (!qx$rank) return(matrix(0, nrow(x), nrow(x)))
  q <- qr.Q(qx, complete = FALSE)
  tcrossprod(q[, seq_len(qx$rank), drop = FALSE])
}

gower_matrix <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  j <- diag(n) - matrix(1 / n, n, n)
  -0.5 * j %*% (distance_matrix^2) %*% j
}

term_stat <- function(g, dat, predictors, term) {
  reduced <- setdiff(predictors, term)
  full_formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  reduced_formula <- if (length(reduced)) {
    as.formula(paste("~", paste(reduced, collapse = " + ")))
  } else {
    ~1
  }
  x_full <- model.matrix(full_formula, dat)
  x_reduced <- model.matrix(reduced_formula, dat)
  h_full <- hat_matrix(x_full)
  h_reduced <- hat_matrix(x_reduced)
  ss_total <- sum(diag(g))
  ss_term <- sum(diag((h_full - h_reduced) %*% g))
  ss_res <- sum(diag((diag(nrow(g)) - h_full) %*% g))
  df_term <- qr(x_full)$rank - qr(x_reduced)$rank
  df_res <- nrow(g) - qr(x_full)$rank
  if (df_term < 1L || df_res < 1L || ss_res <= 0) {
    return(c(R2 = NA_real_, pseudo_F = NA_real_, df_term = df_term, df_residual = df_res))
  }
  c(
    R2 = ss_term / ss_total,
    pseudo_F = (ss_term / df_term) / (ss_res / df_res),
    df_term = df_term,
    df_residual = df_res
  )
}

basin_permutations <- function(dat) {
  sites <- unique(dat[, c("site", "basin")])
  sites$site <- as.character(sites$site)
  n_north <- sum(sites$basin == "north")
  north_sets <- combn(sites$site, n_north, simplify = FALSE)
  lapply(north_sets, function(north) {
    p <- dat
    p$basin <- factor(ifelse(p$site %in% north, "north", "south"), levels = c("north", "south"))
    p
  })
}

paired_response_permutations <- function(dat) {
  sites <- unique(as.character(dat$site))
  grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(sites))))
  lapply(seq_len(nrow(grid)), function(i) {
    idx <- seq_len(nrow(dat))
    for (j in seq_along(sites)) {
      if (grid[i, j]) {
        rows <- which(dat$site == sites[j])
        idx[rows] <- rev(rows)
      }
    }
    idx
  })
}

permutation_test <- function(g, dat, predictors, term, scheme) {
  observed <- term_stat(g, dat, predictors, term)
  if (!is.finite(observed[["pseudo_F"]])) return(c(observed, P = NA_real_, permutations = NA_real_))
  if (scheme == "whole_site_basin_labels") {
    perms <- basin_permutations(dat)
    perm_f <- vapply(perms, function(p) term_stat(g, p, predictors, term)[["pseudo_F"]], numeric(1))
  } else {
    perms <- paired_response_permutations(dat)
    perm_f <- vapply(perms, function(idx) {
      term_stat(g[idx, idx, drop = FALSE], dat, predictors, term)[["pseudo_F"]]
    }, numeric(1))
  }
  c(
    observed,
    P = sum(perm_f >= observed[["pseudo_F"]] - sqrt(.Machine$double.eps), na.rm = TRUE) / sum(is.finite(perm_f)),
    permutations = sum(is.finite(perm_f))
  )
}

model_specs <- list(
  baseline = c("n_detected_mags", "season", "basin"),
  salinity_conditioned = c("n_detected_mags", "season", "salinity_z", "basin"),
  ion_PC1_conditioned = c("n_detected_mags", "season", "ion_PC1_z", "basin"),
  salinity_ion_PC1_conditioned = c("n_detected_mags", "season", "salinity_z", "ion_PC1_z", "basin"),
  salinity_ion_PC1_PC2_sensitivity = c(
    "n_detected_mags", "season", "salinity_z", "ion_PC1_z", "ion_PC2_z", "basin"
  )
)

result_rows <- list()
row_i <- 0L
validation_rows <- list()
validation_i <- 0L

for (response_name in names(responses)) {
  response <- responses[[response_name]]
  ids <- meta$sample_id
  if (!all(ids %in% rownames(response$matrix))) stop("Response lacks n=14 IDs: ", response_name)
  mat <- response$matrix[ids, , drop = FALSE]
  keep <- apply(mat, 2, function(z) all(is.finite(z)) && stats::sd(z) > 0)
  mat <- mat[, keep, drop = FALSE]
  d <- as.matrix(vegdist(mat, method = "bray"))
  g <- gower_matrix(d)
  dat <- meta
  tech <- response$tech[ids]
  if (anyNA(tech)) stop("Missing technical covariate: ", response_name)
  dat$n_detected_mags <- as.numeric(scale(as.numeric(tech)))

  for (model_name in names(model_specs)) {
    predictors <- model_specs[[model_name]]
    terms_to_test <- intersect(c("basin", "salinity_z", "ion_PC1_z", "ion_PC2_z"), predictors)
    for (term in terms_to_test) {
      scheme <- if (term == "basin") "whole_site_basin_labels" else "within_site_response_swaps"
      test <- permutation_test(g, dat, predictors, term, scheme)
      row_i <- row_i + 1L
      result_rows[[row_i]] <- data.frame(
        response = response_name,
        response_class = response$class,
        model = model_name,
        predictors = paste(predictors, collapse = " + "),
        term = term,
        n_samples = nrow(dat),
        n_sites = length(unique(dat$site)),
        n_features = ncol(mat),
        R2 = unname(test[["R2"]]),
        pseudo_F = unname(test[["pseudo_F"]]),
        df_term = unname(test[["df_term"]]),
        df_residual = unname(test[["df_residual"]]),
        P_restricted = unname(test[["P"]]),
        permutations = unname(test[["permutations"]]),
        permutation_scheme = scheme,
        stringsAsFactors = FALSE
      )

      vegan_formula <- as.formula(paste("as.dist(d) ~", paste(predictors, collapse = " + ")))
      vegan_fit <- adonis2(
        vegan_formula, data = dat, by = "margin", permutations = 9
      )
      validation_i <- validation_i + 1L
      validation_rows[[validation_i]] <- data.frame(
        response = response_name,
        model = model_name,
        term = term,
        custom_R2 = unname(test[["R2"]]),
        vegan_R2 = unname(vegan_fit[term, "R2"]),
        R2_absolute_difference = abs(unname(test[["R2"]]) - unname(vegan_fit[term, "R2"])),
        custom_pseudo_F = unname(test[["pseudo_F"]]),
        vegan_pseudo_F = unname(vegan_fit[term, "F"]),
        F_absolute_difference = abs(unname(test[["pseudo_F"]]) - unname(vegan_fit[term, "F"])),
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, result_rows)
results$P_BH_within_model_term <- ave(
  results$P_restricted,
  interaction(results$model, results$term, drop = TRUE),
  FUN = function(x) p.adjust(x, method = "BH")
)
write.table(
  results, file.path(results_dir, "continuous_environment_PERMANOVA_n14.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
validation <- do.call(rbind, validation_rows)
write.table(
  validation, file.path(results_dir, "custom_statistics_validation_against_vegan.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

baseline <- results[results$model == "baseline" & results$term == "basin", c("response", "R2", "P_restricted")]
names(baseline)[2:3] <- c("basin_R2_baseline", "basin_P_baseline")
joint <- results[
  results$model == "salinity_ion_PC1_conditioned" & results$term == "basin",
  c("response", "R2", "P_restricted")
]
names(joint)[2:3] <- c("basin_R2_after_salinity_ion_PC1", "basin_P_after_salinity_ion_PC1")
attenuation <- merge(baseline, joint, by = "response", all = TRUE)
attenuation$absolute_R2_change <- attenuation$basin_R2_after_salinity_ion_PC1 - attenuation$basin_R2_baseline
attenuation$proportional_attenuation <- with(
  attenuation,
  1 - basin_R2_after_salinity_ion_PC1 / basin_R2_baseline
)
write.table(
  attenuation, file.path(results_dir, "basin_attenuation_after_environment_n14.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

basin_by_model <- results[results$term == "basin", c(
  "response", "response_class", "model", "R2", "P_restricted", "P_BH_within_model_term"
)]
basin_by_model <- merge(
  basin_by_model,
  baseline[, c("response", "basin_R2_baseline")],
  by = "response", all.x = TRUE
)
basin_by_model$proportional_attenuation_from_baseline <- with(
  basin_by_model,
  1 - R2 / basin_R2_baseline
)
write.table(
  basin_by_model, file.path(results_dir, "basin_effect_by_environment_model_n14.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cor_vars <- meta[, c("salinity_g_L", "pH", "ion_PC1", "ion_PC2")]
cor_matrix <- cor(cor_vars, method = "spearman")
write.table(
  cbind(variable = rownames(cor_matrix), as.data.frame(cor_matrix)),
  file.path(results_dir, "environment_spearman_correlations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

vif_for <- function(dat, predictors) {
  numeric_predictors <- predictors[predictors %in% c("n_detected_mags", "salinity_z", "ion_PC1_z", "ion_PC2_z")]
  out <- lapply(numeric_predictors, function(term) {
    others <- setdiff(predictors, term)
    fit <- lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = dat)
    r2 <- summary(fit)$r.squared
    data.frame(term = term, VIF = 1 / (1 - r2))
  })
  do.call(rbind, out)
}

tech_reference <- responses$graded_modules_269$tech[meta$sample_id]
vif_dat <- meta
vif_dat$n_detected_mags <- as.numeric(scale(as.numeric(tech_reference)))
vif_rows <- do.call(rbind, lapply(names(model_specs), function(model_name) {
  z <- vif_for(vif_dat, model_specs[[model_name]])
  cbind(model = model_name, z)
}))
write.table(vif_rows, file.path(results_dir, "predictor_VIF_n14.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Compact diagnostic figure: environmental space and basin-term attenuation.
plot_meta <- meta
plot_meta$season <- factor(plot_meta$season, levels = c("spring", "autumn"))
p1 <- ggplot(plot_meta, aes(ion_PC1_z, salinity_z, colour = basin, shape = season)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#B8BEC5") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#B8BEC5") +
  geom_path(aes(group = site), linewidth = 0.45, colour = "#A0A6AD") +
  geom_point(size = 2.8, stroke = 0.7) +
  scale_colour_manual(values = c(north = "#355F8A", south = "#E6A33D")) +
  scale_shape_manual(values = c(spring = 16, autumn = 15)) +
  labs(x = "Major-ion composition PC1 (standardised)", y = "Salinity (standardised)") +
  theme_bw(base_size = 9) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

plot_att <- attenuation[attenuation$response %in% c(
  "MAG_relative_abundance_269", "graded_modules_269", "strict_modules_269", "carrier_fraction_269"
), ]
plot_att <- rbind(
  data.frame(response = plot_att$response, state = "Detection + season", R2 = plot_att$basin_R2_baseline, P = plot_att$basin_P_baseline),
  data.frame(
    response = plot_att$response,
    state = "+ salinity",
    R2 = results$R2[match(
      paste(plot_att$response, "salinity_conditioned", "basin"),
      paste(results$response, results$model, results$term)
    )],
    P = results$P_restricted[match(
      paste(plot_att$response, "salinity_conditioned", "basin"),
      paste(results$response, results$model, results$term)
    )]
  )
)
plot_att$response <- factor(
  plot_att$response,
  levels = rev(c("MAG_relative_abundance_269", "graded_modules_269", "strict_modules_269", "carrier_fraction_269")),
  labels = rev(c("MAG abundance", "Graded modules", "Strict modules", "Carrier fraction"))
)
p2 <- ggplot(plot_att, aes(R2, response, colour = state, shape = P <= 0.05)) +
  geom_line(aes(group = response), colour = "#B8BEC5", linewidth = 0.5) +
  geom_point(size = 2.8, stroke = 0.8) +
  scale_colour_manual(values = c("Detection + season" = "#D18B2C", "+ salinity" = "#2A8C82")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "P <= 0.05", `FALSE` = "P > 0.05")) +
  labs(
    x = expression("Partial basin " * R^2), y = NULL, colour = NULL, shape = NULL,
    subtitle = "Filled: P <= 0.05; open: P > 0.05"
  ) +
  theme_bw(base_size = 9) +
  guides(colour = guide_legend(order = 1), shape = "none") +
  theme(panel.grid = element_blank(), legend.position = "bottom", plot.subtitle = element_text(size = 8))

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- p1 + p2 + patchwork::plot_annotation(tag_levels = "a")
  ggsave(file.path(figures_dir, "continuous_environment_MAG_diagnostic.pdf"), combined, width = 8.2, height = 4.2, device = cairo_pdf)
  ggsave(file.path(figures_dir, "continuous_environment_MAG_diagnostic.png"), combined, width = 8.2, height = 4.2, dpi = 600, bg = "white")
} else {
  ggsave(file.path(figures_dir, "environment_space.pdf"), p1, width = 3.6, height = 3.6, device = cairo_pdf)
  ggsave(file.path(figures_dir, "environment_space.png"), p1, width = 3.6, height = 3.6, dpi = 600, bg = "white")
  ggsave(file.path(figures_dir, "basin_attenuation.pdf"), p2, width = 3.6, height = 3.6, device = cairo_pdf)
  ggsave(file.path(figures_dir, "basin_attenuation.png"), p2, width = 3.6, height = 3.6, dpi = 600, bg = "white")
}

capture.output(sessionInfo(), file = file.path(logs_dir, "sessionInfo.txt"))
writeLines(capture.output(print(results, row.names = FALSE)), file.path(logs_dir, "analysis_console.txt"))
