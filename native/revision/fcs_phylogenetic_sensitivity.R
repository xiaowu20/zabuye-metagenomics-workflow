suppressPackageStartupMessages({library(ape);library(phangorn);library(phytools);library(tibble)})
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
follow<-file.path(root,"deposition_qc/followup_20261112")
out <- if (length(args) >= 2) args[2] else file.path(follow, "phylogenetic_sensitivity")
stopifnot(!dir.exists(out));dir.create(out)
d<-read.delim(file.path(root,"inputs/Table_S25_ZBY_pI_multivariable_input_corrected.tsv"),check.names=FALSE)
d<-d[d$gtdb_phylum=="Halobacteriota",];stopifnot(nrow(d)==84)
actions<-read.delim(file.path(follow,"tables/FCS_action_calls.tsv"))
excluded<-intersect(d$mag_id,actions$mag_id[actions$action=="EXCLUDE"])
flagged<-intersect(d$mag_id,unique(actions$mag_id))
rare <- intersect(d$mag_id, actions$mag_id[actions$action == "REVIEW_RARE"])
stopifnot(!length(excluded), length(flagged) == 1, setequal(flagged, rare))
for(expr in parse(file.path(root,"code/replay_phylogenetic_models.R"))) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&identical(expr[[2]],as.name("fit_manual_pgls")))eval(expr)
}
paths<-c(archival_143_pruned="inputs/ar53_iqtree.treefile",GTDB_masked_AA84="phylogeny/aa84/ar53_BIC.treefile",
         matched_AA84="phylogeny/codon84/AA_BIC.treefile",matched_codon_NT84="phylogeny/codon84/NT_BIC.treefile")
stopifnot(all(file.exists(file.path(root,paths))))
coefs<-list();signals<-list();set.seed(20261112)
for(id in names(paths)) {
  tree<-phangorn::midpoint(read.tree(file.path(root,paths[id])))
  tree$edge.length<-pmax(tree$edge.length,1e-8)
  # Match the archived rooting/pruning order; retain the original standardisation.
  tree<-keep.tip(tree,d$mag_id)
  for(scenario in c("baseline84","omit_REVIEW_RARE_MAG83")) {
    phy<-if(scenario=="baseline84")tree else drop.tip(tree,flagged)
    dat<-d[match(phy$tip.label,d$mag_id),];n<-nrow(dat);df<-n-5L
    stopifnot(identical(as.character(dat$mag_id),phy$tip.label))
    fit<-as.data.frame(fit_manual_pgls(z_pI~z_salt+z_completeness+z_contamination+z_log_n_proteins,dat,phy))
    fit$tree<-id;fit$scenario<-scenario;fit$interval_method<-paste0("archival_ML_scale_conditional_t",df)
    finite<-fit;finite$se<-finite$se*sqrt(n/df);finite$statistic<-finite$estimate/finite$se
    finite$p_value<-2*pt(abs(finite$statistic),df=df,lower.tail=FALSE)
    finite$conf_low<-finite$estimate-qt(.975,df)*finite$se;finite$conf_high<-finite$estimate+qt(.975,df)*finite$se
    finite$interval_method<-paste0("residual_df_scale_conditional_t",df)
    coefs[[length(coefs)+1]]<-rbind(fit,finite)
    sig<-phylosig(phy,setNames(dat$median_pI,dat$mag_id),method="lambda",test=TRUE)
    signals[[length(signals)+1]]<-data.frame(tree=id,scenario=scenario,n=n,lambda=sig$lambda,logL=sig$logL,logL0=sig$logL0,
      LRT=2*(sig$logL-sig$logL0),p_chisq1=sig$P)
    write.tree(phy,file.path(out,paste0(id,"_",scenario,".nwk")))
  }
}
all<-do.call(rbind,coefs)
old<-read.delim(file.path(root,"results/phylogenetic_sensitivity/PGLS_coefficients_by_tree.tsv"))
new<-all[all$scenario=="baseline84",]
idx<-match(paste(new$tree,new$term,new$interval_method),paste(old$tree,old$term,old$interval_method))
stopifnot(!anyNA(idx),max(abs(new$estimate-old$estimate[idx]))<1e-9,max(abs(new$p_value-old$p_value[idx]))<1e-8)
write.table(all,file.path(out,"PGLS_coefficients.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,signals),file.path(out,"Pagel_lambda.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(tree=new$tree,term=new$term,method=new$interval_method,
  estimate_difference=abs(new$estimate-old$estimate[idx]),P_difference=abs(new$p_value-old$p_value[idx])),
  file.path(out,"baseline_reproduction_QA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("The sole FCS action genome in the 84-MAG Halobacteriota analysis is ZBY_water__metabat2_bin.249_sub (REVIEW_RARE).",
  "Whole-MAG omission is a sensitivity analysis, not approval to delete its flagged contig or discard the MAG.",
  "Existing trees are pruned; no new phylogenetic reconstruction is claimed. Original 269-MAG z-score standardisation is preserved."),file.path(out,"scope.txt"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines("complete",file.path(out,"COMPLETE.txt"));print(all[all$term=="z_salt",])
