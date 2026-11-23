suppressPackageStartupMessages({library(vegan);library(dplyr);library(tidyr)})
args<-commandArgs(trailingOnly=TRUE);root<-args[1];group<-args[2]
out<-file.path(root,"catalogues",group,"statistics");dir.create(out,showWarnings=FALSE)
# Import only the frozen method function definitions, with no data reads or analyses.
for(expr in parse(file.path(root,"code/paired_station_method_source.R"))) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&length(expr)==3&&
     is.call(expr[[3]])&&identical(expr[[3]][[1]],as.name("function")))eval(expr)
}
meta<-read.delim(file.path(root,"inputs/paired14.tsv"))
meta$sample_id<-meta$sample;meta$site<-factor(meta$station)
meta$season<-factor(meta$season,levels=c("autumn","spring"))
meta$basin<-factor(ifelse(meta$basin %in% c("NZ","north"),"north","south"),levels=c("north","south"))
meta$coverage_pattern<-factor("paired")
stopifnot(nrow(meta)==14,all(table(meta$site)==2),sum(meta$basin=="north")==8)
qc<-read.delim(file.path(root,"catalogues",group,"responses/sample_QC.tsv"))
set.seed(20261110);res<-list();disp<-list();loso<-list();checks<-list()
for(cutoff in c(50,25,10)) {
  for(response in c("MAG","GMP","SCA","DCP")) {
    mat<-as.matrix(read.delim(file.path(root,"catalogues",group,"responses",paste0(response,"_breadth",cutoff,".tsv")),row.names=1,check.names=FALSE))
    mat<-mat[meta$sample_id,,drop=FALSE]
    mat_full<-mat
    mat<-mat[,apply(mat,2,sd)>0,drop=FALSE]
    stopifnot(ncol(mat)>0,all(rowSums(mat)>0))
    m<-meta
    q<-qc[qc$breadth_cutoff==cutoff/100,]
    tech<-q$n_detected_mags[match(meta$sample_id,q$sample_id)]
    m$n_detected_mags<-if(sd(tech)>0)as.numeric(scale(tech))else 0
    d<-as.matrix(vegdist(mat,method="bray"));g<-gower_matrix(d)
    for(term in c("season","basin")) {
      dp<-dispersion_exact_test(d,m,term)
      disp[[length(disp)+1]]<-cbind(catalogue=group,response=response,breadth=cutoff/100,dp)
      for(adjust in c(FALSE,TRUE)) {
        test<-exact_term_test(g,m,term,adjust)
        boot<-bootstrap_term(d,m,term,adjust,n_boot=999L)
        qci<-quantile(boot$R2,c(.025,.975),na.rm=TRUE)
        res[[length(res)+1]]<-data.frame(catalogue=group,response=response,breadth=cutoff/100,
          technical_adjustment=if(adjust)"detected_MAG_count"else"none",term=term,n_samples=14,n_stations=7,
          variable_features=ncol(mat),R2=test$observed[['R2']],R2_CI_low=qci[1],R2_CI_high=qci[2],
          pseudo_F=test$observed[['pseudo_F']],P_exact=test$p_exact,assignments=test$n_permutations,
          bootstrap_valid=length(boot$R2),bootstrap_requested=999)
        vg<-validate_against_adonis(mat,m,term,adjust)
        checks[[length(checks)+1]]<-data.frame(response=response,breadth=cutoff/100,term=term,adjust=adjust,
          R2_abs_difference=abs(vg['vegan_R2']-test$observed[['R2']]),F_abs_difference=abs(vg['vegan_F']-test$observed[['pseudo_F']]))
        if(cutoff==50&&adjust)for(site in levels(meta$site)) {
          keep<-m$site!=site;ml<-droplevels(m[keep,]);gl<-gower_matrix(d[keep,keep])
          te<-exact_term_test(gl,ml,term,adjust)
          loso[[length(loso)+1]]<-data.frame(catalogue=group,response=response,term=term,omitted_station=site,
            n_stations=6,R2=te$observed[['R2']],pseudo_F=te$observed[['pseudo_F']],P_exact=te$p_exact,assignments=te$n_permutations)
        }
      }
    }
  }
}
write_result<-function(x,name)write.table(bind_rows(x),file.path(out,name),sep="\t",row.names=FALSE,quote=FALSE)
write_result(res,"PERMANOVA.tsv");write_result(disp,"dispersion.tsv")
write_result(loso,"leave_one_station_out.tsv")
source(file.path(root,"code/module_feature_tests.R"))
module_feature_tests(root,group)
write_result(checks,"custom_vs_vegan_statistic_QA.tsv")
stopifnot(max(bind_rows(checks)$R2_abs_difference)<1e-9,max(bind_rows(checks)$F_abs_difference)<1e-8)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(bind_rows(res)%>%filter(breadth==.5,technical_adjustment=="detected_MAG_count"))
writeLines("complete",file.path(out,"COMPLETE.txt"))
