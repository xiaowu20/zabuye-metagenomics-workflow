suppressPackageStartupMessages({library(ape); library(phangorn); library(phytools); library(tibble)})
root <- commandArgs(trailingOnly=TRUE)[1]
out <- file.path(root, "results/phylogenetic_sensitivity")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.delim(file.path(root, "inputs/Table_S25_ZBY_pI_multivariable_input_corrected.tsv"), check.names=FALSE)
d <- d[d$gtdb_phylum == "Halobacteriota", ]
stopifnot(nrow(d) == 84, !anyDuplicated(d$mag_id))
# Reuse the frozen GLS implementation without executing the archival replay script.
for (expr in parse(file.path(root, "code/replay_phylogenetic_models.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("fit_manual_pgls"))) eval(expr)
}
paths <- c(archival_143_pruned="inputs/ar53_iqtree.treefile",
  GTDB_masked_AA84="phylogeny/aa84/ar53_BIC.treefile",
  matched_AA84="phylogeny/codon84/AA_BIC.treefile",
  matched_codon_NT84="phylogeny/codon84/NT_BIC.treefile")
paths <- paths[file.exists(file.path(root, paths))]
trees <- list(); signal <- list(); coefficients <- list(); supports <- list()
set.seed(20261110)
for (id in names(paths)) {
  raw <- read.tree(file.path(root, paths[id]))
  tree <- phangorn::midpoint(raw)
  tree$edge.length <- pmax(tree$edge.length, 1e-8)
  tree <- keep.tip(tree, d$mag_id)
  dat <- d[match(tree$tip.label, d$mag_id), ]
  stopifnot(identical(as.character(dat$mag_id), tree$tip.label))
  trees[[id]] <- tree
  write.tree(tree, file.path(out, paste0(id, "_midpoint.nwk")))
  C <- vcv.phylo(tree, model="Brownian", corr=FALSE)
  write.table(C, file.path(out,paste0(id,"_Brownian_covariance.tsv")),sep="\t",quote=FALSE,col.names=NA)
  sig <- phylosig(tree, setNames(dat$median_pI,dat$mag_id), method="lambda", test=TRUE)
  signal[[id]] <- data.frame(tree=id,n=84,lambda=sig$lambda,logL=sig$logL,logL0=sig$logL0,
    LRT=2*(sig$logL-sig$logL0),p_chisq1=sig$P)
  fit <- as.data.frame(fit_manual_pgls(z_pI ~ z_salt + z_completeness + z_contamination + z_log_n_proteins,dat,tree))
  fit$tree <- id; fit$interval_method <- "archival_ML_scale_conditional_t79"
  finite <- fit
  finite$se <- finite$se*sqrt(84/79)
  finite$statistic <- finite$estimate/finite$se
  finite$p_value <- 2*pt(abs(finite$statistic),df=79,lower.tail=FALSE)
  finite$conf_low <- finite$estimate-qt(.975,79)*finite$se
  finite$conf_high <- finite$estimate+qt(.975,79)*finite$se
  finite$interval_method <- "residual_df_scale_conditional_t79"
  coefficients[[id]] <- rbind(fit,finite)
  # Branch supports are parsed on the actual inferred, unrooted tree.
  internal <- raw$edge[raw$edge[,2] > Ntip(raw),2]
  labels <- raw$node.label[internal-Ntip(raw)]
  parts <- strsplit(labels,"/",fixed=TRUE)
  supports[[id]] <- data.frame(tree=id,inferred_tips=Ntip(raw),node=internal,label=labels,
    SH_aLRT=vapply(parts,function(z) if(length(z)==2) as.numeric(z[1]) else NA_real_,0),
    UFBoot=vapply(parts,function(z) as.numeric(tail(z,1)),0))
}
write_tsv <- function(x,name) write.table(x,file.path(out,name),sep="\t",quote=FALSE,row.names=FALSE)
write_tsv(do.call(rbind,signal),"Pagel_lambda_by_tree.tsv")
write_tsv(do.call(rbind,coefficients),"PGLS_coefficients_by_tree.tsv")
supp <- do.call(rbind,supports)
write_tsv(supp,"inferred_branch_supports.tsv")
ss <- do.call(rbind,lapply(split(supp,supp$tree),function(x) data.frame(tree=x$tree[1],
  inferred_tips=x$inferred_tips[1],internal_branches=nrow(x),UFBoot_median=median(x$UFBoot,na.rm=TRUE),
  UFBoot_min=min(x$UFBoot,na.rm=TRUE),UFBoot_ge95=sum(x$UFBoot>=95,na.rm=TRUE),
  SH_aLRT_ge80_UFBoot_ge95=sum(x$SH_aLRT>=80 & x$UFBoot>=95,na.rm=TRUE))))
write_tsv(ss,"branch_support_summary.tsv")
pairs <- combn(names(trees),2,simplify=FALSE)
cmp <- do.call(rbind,lapply(pairs,function(ids) {
  a<-unroot(trees[[ids[1]]]); b<-unroot(trees[[ids[2]]])
  a<-keep.tip(a,sort(d$mag_id)); b<-keep.tip(b,sort(d$mag_id))
  ca<-cophenetic.phylo(a)[sort(d$mag_id),sort(d$mag_id)]
  cb<-cophenetic.phylo(b)[sort(d$mag_id),sort(d$mag_id)]
  data.frame(tree1=ids[1],tree2=ids[2],n_tips=84,
    RF=as.numeric(RF.dist(a,b,normalize=FALSE,rooted=FALSE,check.labels=TRUE)),
    normalized_RF=as.numeric(RF.dist(a,b,normalize=TRUE,rooted=FALSE,check.labels=TRUE)),
    patristic_Spearman=cor(ca[lower.tri(ca)],cb[lower.tri(cb)],method="spearman"))
}))
write_tsv(cmp,"topology_comparisons.tsv")
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
saveRDS(list(trees=trees,data=d,signal=signal,coefficients=coefficients),file.path(out,"comparison.rds"))
print(ss);print(cmp)
print(do.call(rbind,coefficients)[do.call(rbind,coefficients)$term=="z_salt",])
