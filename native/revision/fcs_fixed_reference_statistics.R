suppressPackageStartupMessages({library(vegan);library(dplyr)})
args<-commandArgs(trailingOnly=TRUE);root<-args[1];original_root<-args[2]
for(expr in parse(file.path(root,"code/paired_station_method_source.R"))) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&length(expr)==3&&
     is.call(expr[[3]])&&identical(expr[[3]][[1]],as.name("function")))eval(expr)
}
source(file.path(root,"code/module_feature_tests.R"))
meta<-read.delim(file.path(root,"inputs/paired14.tsv"))
meta$sample_id<-meta$sample;meta$site<-factor(meta$station)
meta$season<-factor(meta$season,levels=c("autumn","spring"))
meta$basin<-factor(meta$basin,levels=c("north","south"));meta$coverage_pattern<-factor("paired")
stopifnot(nrow(meta)==14,nlevels(meta$site)==7)
np<-read.delim(file.path(original_root,"results/nonpareil/sample_coverage_summary.tsv"))
meta$Nonpareil_coverage<-as.numeric(scale(np$observed_average_coverage[match(meta$sample,np$sample)]))
stopifnot(!anyNA(meta$Nonpareil_coverage))
# Same partial sums of squares and exact restrictions; vary listed covariates only.
term_stat_from_g<-function(g,meta,term,use_tech=FALSE) {
  pred<-c(active_covariates,"season","basin");red<-setdiff(pred,term)
  xf<-model.matrix(reformulate(pred),meta);xr<-model.matrix(reformulate(red),meta)
  hf<-hat_matrix(xf);hr<-hat_matrix(xr);n<-nrow(g)
  ss<-sum(diag((hf-hr)%*%g));resid<-sum(diag((diag(n)-hf)%*%g))
  dt<-qr(xf)$rank-qr(xr)$rank;dr<-n-qr(xf)$rank
  total<-sum(diag(g))
  adjusted<-adjusted_r2(sum(diag(hf%*%g))/total,n,qr(xf)$rank-1L)-
    adjusted_r2(sum(diag(hr%*%g))/total,n,qr(xr)$rank-1L)
  c(pseudo_F=(ss/dt)/(resid/dr),R2=ss/total,adjusted_R2=adjusted,df_term=dt,df_residual=dr)
}
models<-list(detected_MAG_count="n_detected_mags",both=c("n_detected_mags","Nonpareil_coverage"))
groups<-c("baseline_original269","without_EXCLUDE_MAGs","without_all_action_MAGs","EXCLUDE_contigs_function_only")
for(group in groups) {
  message("Statistical scenario: ",group)
  out<-file.path(root,"catalogues",group,"statistics")
  qc<-read.delim(file.path(root,"catalogues",group,"responses/sample_QC.tsv"));qc<-qc[qc$breadth_cutoff==.5,]
  m<-meta;m$n_detected_mags<-as.numeric(scale(qc$n_detected_mags[match(m$sample,qc$sample_id)]))
  stopifnot(!anyNA(m$n_detected_mags))
  result<-list();dispersion<-list();qa<-list();loso<-list();set.seed(20261112)
  for(response in c("MAG","GMP","SCA","DCP")) {
    message("Response: ",response)
    mat<-as.matrix(read.delim(file.path(root,"catalogues",group,"responses",paste0(response,"_breadth50.tsv")),row.names=1,check.names=FALSE))
    mat<-mat[m$sample_id,,drop=FALSE];mat<-mat[,apply(mat,2,sd)>0,drop=FALSE]
    d<-as.matrix(vegdist(mat,method="bray"));g<-gower_matrix(d)
    for(term in c("season","basin")) {
      dispersion[[length(dispersion)+1]]<-cbind(catalogue=group,response=response,breadth=.5,dispersion_exact_test(d,m,term))
      for(model in names(models)) {
        active_covariates<-models[[model]]
        test<-exact_term_test(g,m,term,FALSE)
        stopifnot(all(c("R2","adjusted_R2")%in%names(test$observed)))
        boot<-bootstrap_term(d,m,term,FALSE,n_boot=999L)
        ci<-quantile(boot$R2,c(.025,.975),na.rm=TRUE)
        result[[length(result)+1]]<-data.frame(catalogue=group,response=response,breadth=.5,model=model,term=term,
          n_stations=7,n_samples=14,variable_features=ncol(mat),R2=test$observed[['R2']],CI_low=ci[1],CI_high=ci[2],
          pseudo_F=test$observed[['pseudo_F']],P_exact=test$p_exact,assignments=test$n_permutations,bootstrap_valid=length(boot$R2))
        z<-mat;f<-reformulate(c(active_covariates,"season","basin"),response="z")
        v<-adonis2(f,data=m,method="bray",by="margin",permutations=9)
        qa[[length(qa)+1]]<-data.frame(response=response,model=model,term=term,
          R2_difference=abs(v[term,'R2']-test$observed[['R2']]),F_difference=abs(v[term,'F']-test$observed[['pseudo_F']]))
        if(model=="detected_MAG_count")for(station in levels(m$site)) {
          keep<-m$site!=station;ml<-droplevels(m[keep,]);test_l<-exact_term_test(gower_matrix(d[keep,keep]),ml,term,FALSE)
          loso[[length(loso)+1]]<-data.frame(catalogue=group,response=response,term=term,omitted_station=station,
            R2=test_l$observed[['R2']],pseudo_F=test_l$observed[['pseudo_F']],P_exact=test_l$p_exact,assignments=test_l$n_permutations)
        }
      }
    }
  }
  result<-bind_rows(result);qa<-bind_rows(qa)
  stopifnot(max(qa$R2_difference)<1e-9,max(qa$F_difference)<1e-8)
  write_tsv<-function(x,name)write.table(x,file.path(out,name),sep="\t",row.names=FALSE,quote=FALSE)
  write_tsv(result,"PERMANOVA.tsv");write_tsv(qa,"custom_vs_vegan_QA.tsv")
  write_tsv(bind_rows(dispersion),"dispersion.tsv");write_tsv(bind_rows(loso),"leave_one_station_out.tsv")
  module_feature_tests(root,group)
  if(group=="baseline_original269") {
    old<-read.delim(file.path(original_root,"catalogues/archival_original269/statistics/PERMANOVA.tsv"))
    old<-old[old$breadth==.5&old$technical_adjustment=="detected_MAG_count",]
    new<-result[result$model=="detected_MAG_count",];idx<-match(paste(new$response,new$term),paste(old$response,old$term))
    stopifnot(!anyNA(idx),max(abs(new$R2-old$R2[idx]))<1e-10,max(abs(new$P_exact-old$P_exact[idx]))<1e-10)
    write_tsv(data.frame(response=new$response,term=new$term,R2_difference=abs(new$R2-old$R2[idx]),P_difference=abs(new$P_exact-old$P_exact[idx])),"original_PERMANOVA_reproduction.tsv")
  }
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  writeLines("complete",file.path(out,"COMPLETE.txt"));print(result)
}
