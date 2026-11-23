suppressPackageStartupMessages({
  library(vegan)
  library(tidyverse)
})

set.seed(20260720)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript dom_models.R INPUTS OUTPUT")
IN <- args[1]; IN2 <- IN; OUT <- args[2]
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ── 读取并合并DOM/离子轴、检出比例、样本manifest ──
manifest <- read_tsv(file.path(IN,"ZBYinternal_manifest_with_site.tsv"), show_col_types=FALSE)
qc_detfrac <- read_tsv(file.path(IN2,"sample_detection_qc.tsv"), show_col_types=FALSE)
dom_axes   <- read_tsv(file.path(IN2,"dom_environmental_axes_v2.tsv"), show_col_types=FALSE)

# match_id 已在dom_environmental_axes_v2.tsv中正确构造(season为Autumn/Spring/Summer/Winter全称大写)

meta_full <- manifest %>%
  left_join(qc_detfrac, by="sample_id") %>%
  mutate(basin=factor(basin, levels=c("north","south")),
         season=factor(season, levels=c("spring","summer","autumn","winter")),
         site=factor(site))

meta_full <- meta_full %>%
  left_join(dom_axes %>% select(match_id, logSalinity, IonComp_PC1,
                                logDOC, DOM_quality_PC1),
            by=c("matched_environment_id"="match_id"))

cat("=== 合并后DOM变量缺失情况 ===\n")
print(meta_full %>% filter(in_main_model==TRUE | in_main_model=="True") %>%
      select(sample_id, basin, season, matched_environment_id,
             logDOC, DOM_quality_PC1) %>%
      mutate(has_dom = !is.na(logDOC)))

# ── 读取响应矩阵 ──
graded <- read_tsv(file.path(IN,"ZBYinternal_response_graded_weighted.tsv"), show_col_types=FALSE)
colnames(graded)[1] <- "sample_id"

run_partial_dbrda <- function(sample_ids, label, use_strata_site=FALSE) {
  m <- meta_full %>% filter(sample_id %in% sample_ids, !is.na(logDOC)) %>%
    droplevels()
  resp <- graded %>% filter(sample_id %in% m$sample_id) %>%
    arrange(match(sample_id, m$sample_id))
  m <- m %>% arrange(match(sample_id, resp$sample_id))
  mat <- as.matrix(resp %>% select(-sample_id))
  rownames(mat) <- resp$sample_id
  mat <- mat[, colSums(mat) > 0]  # 删除全零模块

  cat(sprintf("\n\n========== %s (n=%d) ==========\n", label, nrow(m)))
  cat("零方差模块过滤后剩余模块数:", ncol(mat), "\n")

  if (nrow(m) < 8) {
    cat("样本量过小, 跳过\n"); return(NULL)
  }

  # 置换设计
  if (use_strata_site && length(unique(m$site)) > 1) {
    perm_ctrl <- how(blocks = m$site, within=Within(type="free"), nperm=9999)
  } else {
    perm_ctrl <- 9999
  }

  # ── Bray-Curtis partial dbRDA (主分析) ──
  mod_full <- capscale(
    mat ~ logDOC + DOM_quality_PC1 +
      Condition(logSalinity + IonComp_PC1 + basin + season + catalogue_detection_fraction_z),
    data = m, distance = "bray", sqrt.dist = TRUE)

  cat("\n--- Bray-Curtis partial dbRDA: DOM组整体检验 ---\n")
  anova_full <- tryCatch(anova(mod_full, permutations=perm_ctrl),
                         error=function(e){cat("检验失败:",conditionMessage(e),"\n"); NULL})
  if (!is.null(anova_full)) print(anova_full)

  r2adj_full <- tryCatch(RsquareAdj(mod_full)$adj.r.squared, error=function(e) NA)
  cat("DOM组独立(条件化后)adjusted R2:", r2adj_full, "\n")

  # DOM总效应(不条件化, 用于估算共享解释量)
  mod_dom_only <- capscale(mat ~ logDOC + DOM_quality_PC1, data=m,
                           distance="bray", sqrt.dist=TRUE)
  r2adj_total <- tryCatch(RsquareAdj(mod_dom_only)$adj.r.squared, error=function(e) NA)
  shared_approx <- r2adj_total - r2adj_full
  cat("DOM组总(未条件化)adjusted R2:", r2adj_total, "\n")
  cat("近似共享解释量(总-独立):", shared_approx, "(注: 近似值, 非严格variance partition)\n")

  # ── Hellinger变换 + RDA (敏感性) ──
  mat_hell <- decostand(mat, method="hellinger")
  mod_full_hell <- rda(
    mat_hell ~ logDOC + DOM_quality_PC1 +
      Condition(logSalinity + IonComp_PC1 + basin + season + catalogue_detection_fraction_z),
    data = m)
  cat("\n--- Hellinger + RDA: DOM组整体检验(敏感性) ---\n")
  anova_hell <- tryCatch(anova(mod_full_hell, permutations=perm_ctrl),
                        error=function(e){cat("检验失败:",conditionMessage(e),"\n"); NULL})
  if (!is.null(anova_hell)) print(anova_hell)
  r2adj_hell <- tryCatch(RsquareAdj(mod_full_hell)$adj.r.squared, error=function(e) NA)
  cat("Hellinger敏感性 adjusted R2:", r2adj_hell, "\n")

  list(bray_anova=anova_full, bray_r2adj=r2adj_full, bray_r2adj_total=r2adj_total,
       hellinger_anova=anova_hell, hellinger_r2adj=r2adj_hell, n=nrow(m))
}

# ── 三层样本设计 ──
main <- meta_full %>% filter(in_main_model==TRUE | in_main_model=="True")

ids_15 <- main %>% filter(!is.na(site), season %in% c("spring","autumn")) %>% pull(sample_id)
ids_19 <- main %>% filter(!is.na(site)) %>% pull(sample_id)

cat("15样本(主分析, 环境匹配春秋):", length(ids_15), "\n")
cat("19样本(敏感性, 全季节环境匹配):", length(ids_19), "\n")

res_15 <- run_partial_dbrda(ids_15, "15样本主分析(春秋, basin_site限制置换)", use_strata_site=TRUE)
res_19 <- run_partial_dbrda(ids_19, "19样本敏感性(全季节, 不检验交互)", use_strata_site=FALSE)

saveRDS(list(res_15=res_15, res_19=res_19),
        file.path(OUT, "partial_dbrda_results.rds"))
cat("\n\n完成, 结果已保存\n")
