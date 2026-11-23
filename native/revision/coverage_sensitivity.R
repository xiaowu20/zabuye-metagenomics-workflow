suppressPackageStartupMessages({library(vegan);library(dplyr)})
args<-commandArgs(trailingOnly=TRUE);root<-args[1];group<-args[2]
out<-file.path(root,"results/coverage_sensitivity",group);dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(expr in parse(file.path(root,"code/paired_station_method_source.R"))) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&length(expr)==3&&
     is.call(expr[[3]])&&identical(expr[[3]][[1]],as.name("function")))eval(expr)
}
meta<-read.delim(file.path(root,"inputs/paired14.tsv"))
meta$sample_id<-meta$sample;meta$site<-factor(meta$station)
meta$season<-factor(meta$season,levels=c("autumn","spring"))
meta$basin<-factor(meta$basin,levels=c("north","south"));meta$coverage_pattern<-factor("paired")
np<-read.delim(file.path(root,"results/nonpareil/sample_coverage_summary.tsv"))
meta$Nonpareil_coverage<-as.numeric(scale(np$observed_average_coverage[match(meta$sample,np$sample)]))
stopifnot(!anyNA(meta$Nonpareil_coverage))
responses<-list()
if(group=="read_species") {
  counts<-as.matrix(read.delim(file.path(root,"inputs/dual_species_counts.tsv"),row.names=1,check.names=FALSE))
  counts<-counts[,meta$sample_id,drop=FALSE]
  counts<-counts[rowSums(counts>0)>=2,,drop=FALSE]
  counts<-counts[apply(counts,1,sd)>0,,drop=FALSE]
  stopifnot(nrow(counts)>0, all(colSums(counts)>0))
  responses$read_species<-t(sweep(counts,2,colSums(counts),"/"))
} else {
  for(r in c("MAG","GMP","SCA","DCP"))responses[[r]]<-as.matrix(read.delim(file.path(root,"catalogues",group,"responses",paste0(r,"_breadth50.tsv")),row.names=1,check.names=FALSE))
  qc<-read.delim(file.path(root,"catalogues",group,"responses/sample_QC.tsv"))
  qc<-qc[qc$breadth_cutoff==.5,]
  meta$n_detected_mags<-as.numeric(scale(qc$n_detected_mags[match(meta$sample_id,qc$sample_id)]))
}
# Same partial sums of squares and restricted label enumeration as the archived
# method. Only the explicitly listed technical covariates change here.
term_stat_from_g<-function(g,meta,term,use_tech=FALSE) {
  pred<-c(active_covariates,"season","basin");red<-setdiff(pred,term)
  xf<-model.matrix(reformulate(pred),meta);xr<-model.matrix(reformulate(red),meta)
  hf<-hat_matrix(xf);hr<-hat_matrix(xr);n<-nrow(g)
  ss<-sum(diag((hf-hr)%*%g));resid<-sum(diag((diag(n)-hf)%*%g))
  dt<-qr(xf)$rank-qr(xr)$rank;dr<-n-qr(xf)$rank
  c(pseudo_F=(ss/dt)/(resid/dr),R2=ss/sum(diag(g)),df_term=dt,df_residual=dr)
}
models<-list(unadjusted=character(),Nonpareil_coverage="Nonpareil_coverage")
if(group!="read_species")models<-c(models,list(detected_MAG_count="n_detected_mags",both=c("n_detected_mags","Nonpareil_coverage")))
res<-list();checks<-list();set.seed(20261110)
for(response in names(responses)) {
  mat<-responses[[response]][meta$sample_id,,drop=FALSE]
  mat<-mat[,apply(mat,2,sd)>0,drop=FALSE]
  d<-as.matrix(vegdist(mat,method="bray"))
  for(model in names(models)) {
    active_covariates<-models[[model]]
    for(scope in c("paired14","omit_station_N3")) {
      keep<-if(scope=="paired14")rep(TRUE,14)else meta$station!="N3"
      m<-droplevels(meta[keep,]);dm<-d[keep,keep];g<-gower_matrix(dm)
      for(term in c("season","basin")) {
        test<-exact_term_test(g,m,term,FALSE)
        b<-rep(NA_real_,999)
        for(i in seq_along(b)) {
          draw<-bootstrap_indices(m);mb<-m[draw$index,];mb$site<-factor(draw$site)
          b[i]<-term_stat_from_g(gower_matrix(dm[draw$index,draw$index]),mb,term)[['R2']]
        }
        ci<-quantile(b,c(.025,.975),na.rm=TRUE)
        res[[length(res)+1]]<-data.frame(catalogue=group,response=response,model=model,scope=scope,term=term,
          n_samples=nrow(m),n_stations=nlevels(m$site),variable_features=ncol(mat),R2=test$observed[['R2']],
          CI_low=ci[1],CI_high=ci[2],pseudo_F=test$observed[['pseudo_F']],P_exact=test$p_exact,
          permutations=test$n_permutations,bootstrap_valid=sum(is.finite(b)))
        z<-mat[keep,,drop=FALSE];f<-reformulate(c(active_covariates,"season","basin"),response="z")
        v<-adonis2(f,data=m,method="bray",by="margin",permutations=9)
        checks[[length(checks)+1]]<-data.frame(response=response,model=model,scope=scope,term=term,
          R2_abs_difference=abs(v[term,'R2']-test$observed[['R2']]),F_abs_difference=abs(v[term,'F']-test$observed[['pseudo_F']]))
      }
    }
  }
}
r<-bind_rows(res);qa<-bind_rows(checks)
stopifnot(max(qa$R2_abs_difference)<1e-9,max(qa$F_abs_difference)<1e-8)
if(group=="read_species") {
  baseline<-r[r$model=="unadjusted"&r$scope=="paired14"&r$term=="season",]
  stopifnot(abs(baseline$R2-.261804801327055)<1e-10,baseline$P_exact==.015625)
}
write.table(r,file.path(out,"PERMANOVA_Nonpareil_sensitivity.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(qa,file.path(out,"custom_vs_vegan_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(r);writeLines("complete",file.path(out,"COMPLETE.txt"))
