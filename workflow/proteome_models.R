suppressPackageStartupMessages({library(tidyverse); library(lme4); library(nlme); library(ape); library(phangorn); library(phytools); library(car); library(lmtest); library(sandwich)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript proteome_models.R INPUTS OUTPUT")
data_dir <- args[1]; table_dir <- args[2]; dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
options(contrasts=c("contr.treatment","contr.poly"))
set.seed(20260729)
safe_ci <- function(model, term) {
  ci <- suppressMessages(confint(model, parm = term, method = "Wald"))
  c(lo = unname(ci[1]), hi = unname(ci[2]))
}

tidy_lm_term <- function(model, term, model_name, n_units, unit_label) {
  co <- summary(model)$coefficients[term, ]
  ci <- confint(model, term)[1, ]
  n_obs <- length(residuals(model))
  r2 <- summary(model)$r.squared
  tibble(
    model = model_name,
    term = term,
    estimate = unname(co[1]),
    se = unname(co[2]),
    conf_low = unname(ci[1]),
    conf_high = unname(ci[2]),
    statistic = unname(co[3]),
    p_value = unname(co[4]),
    n = n_obs,
    independent_units = n_units,
    unit = unit_label,
    fit_metric = r2,
    fit_metric_name = "R2"
  )
}

## Analysis 2: confirmed per-protein hits, with inference at MAG level
gene_path <- file.path(data_dir, "ZBY_salt_gene_pI_vs_background.tsv")
hits_path <- file.path(data_dir,"target_ko_all_hits.tsv")
gene_all <- read_tsv(gene_path, show_col_types = FALSE)
hits <- read_tsv(
  hits_path,
  col_names = c("star", "protein_id", "KO", "score", "threshold", "evalue", "description"),
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  filter(star == "*", KO %in% c("K03498", "K03499", "K01546", "K01547", "K01548")) %>%
  distinct(protein_id, KO, .keep_all = TRUE)

gene_confirmed <- gene_all %>%
  inner_join(hits %>% select(protein_id, KO, score, threshold, evalue),
             by = c("protein_id", "KO")) %>%
  mutate(
    gene = dplyr::recode(KO,
      K03498 = "K03498 (TrkH/KtrB-like)",
      K03499 = "K03499 (TrkA/KtrA-like)",
      K01546 = "K01546 (KdpA)",
      K01547 = "K01547 (KdpB)",
      K01548 = "K01548 (KdpC)"
    ),
    system = if_else(KO %in% c("K01546", "K01547", "K01548"), "Kdp components", gene)
  )

protein_counts <- gene_confirmed %>%
  group_by(KO, gene) %>%
  summarise(
    n_proteins = n_distinct(protein_id),
    n_mags = n_distinct(mag_id),
    median_length = median(length),
    median_pI = median(pI),
    median_delta = median(pI_delta),
    .groups = "drop"
  )

mag_delta_gene <- gene_confirmed %>%
  group_by(mag_id, KO, gene) %>%
  summarise(
    n_hits = n_distinct(protein_id),
    median_delta = median(pI_delta),
    median_hit_pI = median(pI),
    .groups = "drop"
  )

mag_delta_system <- gene_confirmed %>%
  mutate(system = if_else(KO %in% c("K01546", "K01547", "K01548"),
                          "Kdp components", gene)) %>%
  group_by(mag_id, system) %>%
  summarise(
    n_hits = n_distinct(protein_id),
    median_delta = median(pI_delta),
    median_hit_pI = median(pI),
    .groups = "drop"
  )

wilcox_mag <- bind_rows(
  mag_delta_gene %>% rename(group = gene) %>% select(mag_id, group, median_delta),
  mag_delta_system %>% filter(system == "Kdp components") %>%
    rename(group = system) %>% select(mag_id, group, median_delta)
) %>%
  group_by(group) %>%
  group_modify(~{
    wt <- wilcox.test(.x$median_delta, mu = 0, exact = FALSE, conf.int = TRUE)
    tibble(
      n_mags = nrow(.x),
      median_delta = median(.x$median_delta),
      mean_delta = mean(.x$median_delta),
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      conf_low = unname(wt$conf.int[1]),
      conf_high = unname(wt$conf.int[2])
    )
  }) %>%
  ungroup() %>%
  mutate(p_BH = p.adjust(p_value, method = "BH"))

write_tsv(gene_confirmed, file.path(table_dir, "Table_S22_confirmed_salt_protein_hits.tsv"))
write_tsv(protein_counts, file.path(table_dir, "Table_S23_confirmed_salt_protein_summary.tsv"))
write_tsv(mag_delta_gene, file.path(table_dir, "per_MAG_gene_delta.tsv"))
write_tsv(mag_delta_system, file.path(table_dir, "per_MAG_system_delta.tsv"))
write_tsv(wilcox_mag, file.path(table_dir, "Table_S24_per_MAG_pI_offset_tests.tsv"))

## Analysis 3: corrected module scores, continuous MAG-quality controls, and organism-tree sensitivity
mv_path <- file.path(data_dir, "ZBY_pI_multivar_input.tsv")
mag_path <- file.path(data_dir, "mag_catalogue_414_unified.tsv")
ko_path <- file.path(data_dir, "new_ko_completeness_long.tsv")
tree_path <- file.path(data_dir,"ar53_iqtree.treefile")

pI_mag <- read_tsv(mv_path, show_col_types = FALSE) %>%
  select(mag_id, n_proteins, median_pI, mean_pI, mean_acidic_fraction,
         mean_basic_fraction, acidic_minus_basic)
mag <- read_tsv(mag_path, show_col_types = FALSE) %>%
  filter(source_assembly == "ZBY_water") %>%
  select(mag_id, completeness, contamination, genome_size_bp, gc_content, n50,
         n_genes, quality_tier, gtdb_phylum, gtdb_class, gtdb_order, gtdb_family, gtdb_genus)
ko <- read_tsv(ko_path, show_col_types = FALSE)

salt_modules <- c(
  "Salt_Trk_Ktr_core", "Salt_Trk_Ktr_regulatory", "Salt_Mrp_complex",
  "Salt_betaine_transport", "Salt_betaine_synthesis", "Salt_ectoine_synthesis",
  "Salt_trehalose_synthesis", "Salt_ProU_transporter", "Salt_NaH_antiport_family",
  "Salt_hydroxyectoine_ext", "Salt_Kdp_complex"
)

salt_wide <- ko %>%
  filter(module_id %in% salt_modules) %>%
  mutate(score = case_when(
    potential_status == "complete" ~ 1,
    potential_status %in% c("near_complete", "present") ~ 0.5,
    TRUE ~ 0
  )) %>%
  select(mag_id, module_id, score) %>%
  distinct() %>%
  pivot_wider(names_from = module_id, values_from = score, values_fill = 0)

df <- pI_mag %>%
  inner_join(mag, by = "mag_id") %>%
  inner_join(salt_wide, by = "mag_id") %>%
  mutate(
    salt_composite = rowSums(across(all_of(salt_modules))),
    is_halo = factor(if_else(gtdb_phylum == "Halobacteriota", "Halobacteriota", "Other"),
                     levels = c("Other", "Halobacteriota")),
    gtdb_phylum = fct_explicit_na(factor(gtdb_phylum), na_level = "Unclassified"),
    log_n_proteins = log10(n_proteins),
    z_pI = as.numeric(scale(median_pI)),
    z_salt = as.numeric(scale(salt_composite)),
    z_completeness = as.numeric(scale(completeness)),
    z_contamination = as.numeric(scale(contamination)),
    z_log_n_proteins = as.numeric(scale(log_n_proteins)),
    z_n50 = as.numeric(scale(log10(n50)))
  )

m_unadjusted <- lm(z_pI ~ z_salt, data = df)
m_binary <- lm(z_pI ~ z_salt + is_halo + z_completeness + z_contamination + z_log_n_proteins,
               data = df)
m_phylum <- lmer(
  z_pI ~ z_salt + z_completeness + z_contamination + z_log_n_proteins + (1 | gtdb_phylum),
  data = df, REML = FALSE, control = lmerControl(optimizer = "bobyqa")
)
m_phylum0 <- update(m_phylum, . ~ . - z_salt)
m_phylum_n50 <- update(m_phylum, . ~ . + z_n50)
m_phylum_n50_0 <- update(m_phylum_n50, . ~ . - z_salt)

extract_robust <- function(model, term, name) {
  ct <- coeftest(model, vcov. = vcovHC(model, type = "HC3"))
  est <- ct[term, 1]
  se <- ct[term, 2]
  n_obs <- length(residuals(model))
  cook <- cooks.distance(model)
  df_res <- df.residual(model)
  adj_r2 <- summary(model)$adj.r.squared
  tibble(
    model = name,
    term = term,
    estimate = est,
    se = se,
    conf_low = est - qt(0.975, df_res) * se,
    conf_high = est + qt(0.975, df_res) * se,
    statistic = ct[term, 3],
    p_value = ct[term, 4],
    n = n_obs,
    adjusted_R2 = adj_r2,
    max_Cooks_D = max(cook),
    n_Cooks_gt_4_over_n = sum(cook > 4 / n_obs)
  )
}

extract_lmer_lrt <- function(model, null_model, term, name) {
  co <- fixef(model)[term]
  se <- sqrt(diag(vcov(model)))[term]
  ci <- safe_ci(model, term)
  lrt <- anova(null_model, model)
  n_obs <- nrow(model.frame(model))
  tibble(
    model = name,
    term = term,
    estimate = unname(co),
    se = unname(se),
    conf_low = unname(ci["lo"]),
    conf_high = unname(ci["hi"]),
    statistic = unname(lrt$Chisq[2]),
    p_value = unname(lrt$`Pr(>Chisq)`[2]),
    n = n_obs,
    adjusted_R2 = NA_real_,
    max_Cooks_D = NA_real_,
    n_Cooks_gt_4_over_n = NA_integer_
  )
}

multivar_core <- bind_rows(
  extract_robust(m_unadjusted, "z_salt", "Unadjusted"),
  extract_robust(m_binary, "z_salt", "Adjusted: Halobacteriota contrast + continuous MAG quality"),
  extract_lmer_lrt(m_phylum, m_phylum0, "z_salt", "Adjusted: phylum random intercept + continuous MAG quality"),
  extract_lmer_lrt(m_phylum_n50, m_phylum_n50_0, "z_salt", "Sensitivity: phylum random intercept + quality + N50")
)

quality_terms <- c("z_salt", "is_haloHalobacteriota", "z_completeness", "z_contamination", "z_log_n_proteins")
multivar_coefficients <- map_dfr(quality_terms, function(term) {
  extract_robust(m_binary, term, "Primary adjusted model")
}) %>%
  mutate(p_BH = p.adjust(p_value, method = "BH"))

module_effects <- map_dfr(salt_modules, function(module) {
  if (!module %in% names(df) || sd(df[[module]]) == 0) return(NULL)
  tmp <- df %>% mutate(z_module = as.numeric(scale(.data[[module]])))
  model <- lmer(
    z_pI ~ z_module + z_completeness + z_contamination + z_log_n_proteins + (1 | gtdb_phylum),
    data = tmp, REML = FALSE, control = lmerControl(optimizer = "bobyqa")
  )
  null_model <- update(model, . ~ . - z_module)
  extract_lmer_lrt(model, null_model, "z_module", module) %>% mutate(module = module)
}) %>%
  mutate(p_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_BH)

vif_binary <- car::vif(m_binary)
vif_table <- tibble(term = names(vif_binary), VIF = as.numeric(vif_binary))

## Proper organism-tree sensitivity within Halobacteriota.
tree <- phangorn::midpoint(read.tree(tree_path))
tree$edge.length <- pmax(tree$edge.length, 1e-8)
halo <- df %>% filter(is_halo == "Halobacteriota", mag_id %in% tree$tip.label)
tree_halo <- drop.tip(tree, setdiff(tree$tip.label, halo$mag_id))
halo <- halo %>% arrange(match(mag_id, tree_halo$tip.label))

lambda_result <- tryCatch(
  phytools::phylosig(tree_halo, setNames(halo$median_pI, halo$mag_id), method = "lambda", test = TRUE),
  error = function(e) e
)

fit_manual_pgls <- function(formula, data, phy) {
  data <- as.data.frame(data)
  rownames(data) <- data$mag_id
  data <- data[phy$tip.label, , drop = FALSE]
  X <- model.matrix(formula, data = data)
  y <- model.response(model.frame(formula, data = data))
  C <- vcv.phylo(phy, model = "Brownian", corr = FALSE)
  C <- C[rownames(data), rownames(data)]
  n <- length(y)
  p <- ncol(X)
  evaluate <- function(lambda, return_fit = FALSE) {
    V <- lambda * C
    diag(V) <- diag(C)
    V <- V + diag(1e-8, n)
    Vi <- solve(V)
    XtViX <- crossprod(X, Vi %*% X)
    beta <- solve(XtViX, crossprod(X, Vi %*% y))
    resid <- y - as.vector(X %*% beta)
    rss <- as.numeric(crossprod(resid, Vi %*% resid))
    sigma2 <- rss / n
    logdet <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
    nll <- 0.5 * (n * log(2 * pi) + n * log(sigma2) + logdet + n)
    if (!return_fit) return(nll)
    vc_beta <- sigma2 * solve(XtViX)
    se <- sqrt(diag(vc_beta))
    t_value <- as.vector(beta) / se
    p_value <- 2 * pt(abs(t_value), df = n - p, lower.tail = FALSE)
    tibble(
      term = colnames(X), estimate = as.vector(beta), se = se,
      statistic = t_value, p_value = p_value,
      conf_low = as.vector(beta) - qt(0.975, n - p) * se,
      conf_high = as.vector(beta) + qt(0.975, n - p) * se,
      n = n, lambda = lambda, sigma2 = sigma2, logLik = -nll
    )
  }
  opt <- optimize(evaluate, interval = c(0, 1), tol = 1e-6)
  evaluate(opt$minimum, return_fit = TRUE)
}

pgls_result <- tryCatch(
  fit_manual_pgls(
    z_pI ~ z_salt + z_completeness + z_contamination + z_log_n_proteins,
    halo, tree_halo
  ),
  error = function(e) structure(list(error = conditionMessage(e)), class = "pgls_error")
)

if (!inherits(lambda_result, "error")) {
  lambda_table <- tibble(
    n_tips = length(tree_halo$tip.label),
    lambda = lambda_result$lambda,
    logL = lambda_result$logL,
    p_value = lambda_result$P,
    method = "Pagel lambda on ar53 organism tree"
  )
} else {
  lambda_table <- tibble(
    n_tips = length(tree_halo$tip.label), lambda = NA_real_, logL = NA_real_, p_value = NA_real_,
    method = paste("FAILED:", conditionMessage(lambda_result))
  )
}

if (!inherits(pgls_result, "pgls_error")) {
  pgls_table <- pgls_result
} else {
  pgls_table <- tibble(
    term = "MODEL_FAILED", estimate = NA_real_, se = NA_real_, statistic = NA_real_,
    p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, n = nrow(halo),
    lambda = NA_real_, sigma2 = NA_real_, logLik = NA_real_, error = pgls_result$error
  )
}

write_tsv(df, file.path(table_dir, "Table_S25_ZBY_pI_multivariable_input_corrected.tsv"))
write_tsv(multivar_core, file.path(table_dir, "Table_S26_pI_salt_score_model_comparison.tsv"))
write_tsv(multivar_coefficients, file.path(table_dir, "Table_S27_pI_primary_model_coefficients.tsv"))
write_tsv(module_effects, file.path(table_dir, "Table_S28_module_specific_pI_effects_BH.tsv"))
write_tsv(vif_table, file.path(table_dir, "Table_S29_multivariable_VIF.tsv"))
write_tsv(lambda_table, file.path(table_dir, "Table_S30_Halobacteriota_phylogenetic_signal.tsv"))
write_tsv(pgls_table, file.path(table_dir, "Table_S31_Halobacteriota_PGLS.tsv"))
