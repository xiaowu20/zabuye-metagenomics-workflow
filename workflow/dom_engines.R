args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1 || length(args) > 2) stop("Usage: Rscript dom_engines.R INPUTS [OUTPUT]")
inputs <- args[1]
suppressPackageStartupMessages({library(vegan);library(tidyverse)})
set.seed(20260720)
cat("R_VERSION",R.version.string,"\n");cat("VEGAN_VERSION",as.character(packageVersion("vegan")),"\n")
manifest <- read_tsv(file.path(inputs,"manifest.tsv"),show_col_types=FALSE)
qc_detfrac <- read_tsv(file.path(inputs,"detection.tsv"),show_col_types=FALSE)
dom_axes <- read_tsv(file.path(inputs,"axes.tsv"),show_col_types=FALSE)
meta_full <- manifest %>%
  left_join(qc_detfrac, by="sample_id") %>%
  mutate(basin=factor(basin, levels=c("north","south")),
         season=factor(season, levels=c("spring","summer","autumn","winter")),
         site=factor(site))

meta_full <- meta_full %>%
  left_join(dom_axes %>% select(match_id, logSalinity, IonComp_PC1,
                                logDOC, DOM_quality_PC1),
            by=c("matched_environment_id"="match_id"))

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

main<-meta_full %>% filter(in_main_model==TRUE | in_main_model=="True")
ids_15<-main %>% filter(!is.na(site),season %in% c("spring","autumn")) %>% pull(sample_id)
out<-list()
graded<-read_tsv(file.path(inputs,"preTrkFix.tsv"),show_col_types=FALSE);names(graded)[1]<-"sample_id"
res<-run_partial_dbrda(ids_15,"archived_preTrkFix",TRUE)
out[["archived_preTrkFix"]]<-data.frame(version="archived_preTrkFix",method=c("bray","hellinger"),adjusted_R2=c(res$bray_r2adj,res$hellinger_r2adj),P=c(res$bray_anova[1,"Pr(>F)"],res$hellinger_anova[1,"Pr(>F)"]))
graded<-read_tsv(file.path(inputs,"previous.tsv"),show_col_types=FALSE);names(graded)[1]<-"sample_id"
res<-run_partial_dbrda(ids_15,"legacy_current",TRUE)
out[["legacy_current"]]<-data.frame(version="legacy_current",method=c("bray","hellinger"),adjusted_R2=c(res$bray_r2adj,res$hellinger_r2adj),P=c(res$bray_anova[1,"Pr(>F)"],res$hellinger_anova[1,"Pr(>F)"]))
graded<-read_tsv(file.path(inputs,"annotation_only.tsv"),show_col_types=FALSE);names(graded)[1]<-"sample_id"
res<-run_partial_dbrda(ids_15,"new_annotation_only",TRUE)
out[["new_annotation_only"]]<-data.frame(version="new_annotation_only",method=c("bray","hellinger"),adjusted_R2=c(res$bray_r2adj,res$hellinger_r2adj),P=c(res$bray_anova[1,"Pr(>F)"],res$hellinger_anova[1,"Pr(>F)"]))
graded<-read_tsv(file.path(inputs,"corrected.tsv"),show_col_types=FALSE);names(graded)[1]<-"sample_id"
res<-run_partial_dbrda(ids_15,"new_empty_member_fix",TRUE)
out[["new_empty_member_fix"]]<-data.frame(version="new_empty_member_fix",method=c("bray","hellinger"),adjusted_R2=c(res$bray_r2adj,res$hellinger_r2adj),P=c(res$bray_anova[1,"Pr(>F)"],res$hellinger_anova[1,"Pr(>F)"]))
cat("RESULT_TSV_BEGIN\n");write.table(bind_rows(out),stdout(),sep="\t",row.names=FALSE,quote=FALSE);cat("RESULT_TSV_END\n")
if (length(args) == 2) {
  dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
  write_tsv(bind_rows(out),file.path(args[2],"DOM_engine_replay.tsv"))
  capture.output(sessionInfo(),file=file.path(args[2],"sessionInfo.txt"))
}
