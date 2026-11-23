args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop("Usage: Rscript crosslake.R INPUTS OUTPUT")
inputs <- args[1]; target <- args[2]; version <- "new_empty_member_fix"
dir.create(file.path(target,"results"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(target,"logs"),recursive=TRUE,showWarnings=FALSE)
# Same cross-lake RNG seed, response sequence, filtering, formula sequence and
# 999 free permutations as final script 18_permanova_crosslake_modules.slurm.
suppressPackageStartupMessages({library(vegan);library(tidyverse)})
set.seed(20260724)
qc <- read_tsv(file.path(inputs,'crosslake26_cat414_sample_qc.tsv'),show_col_types=FALSE)
mani <- read_tsv(file.path(inputs,'ZBYinternal_manifest_with_site.tsv'),show_col_types=FALSE)
crossrows <- list(); dispersion_rows <- list()
run <- function(tag,form_list) {
  for(k in c('graded_weighted','strict_carrier_abundance','detected_carrier_fraction')) {
    d <- read_tsv(file.path(inputs,sprintf('%s_response_%s.tsv',tag,k)),show_col_types=FALSE);names(d)[1]<-'sample_id'
    m <- qc %>% filter(sample_id %in% d$sample_id) %>% arrange(match(sample_id,d$sample_id))
    if(tag=='zbyonly21_cat414') m <- m %>% left_join(mani %>% select(sample_id,basin,season),by='sample_id') %>% mutate(basin=factor(basin),season=factor(season))
    m$lake<-factor(m$lake)
    mat<-as.matrix(d[,-1]);rownames(mat)<-d$sample_id;mat<-mat[,colSums(mat)>0,drop=FALSE]
    stopifnot(identical(d$sample_id,m$sample_id))
    for(f in form_list) {
      r<-adonis2(update(f,mat ~ .),data=m,method='bray',by='margin',permutations=999)
      for(term in setdiff(rownames(r),c('Residual','Total'))) crossrows[[length(crossrows)+1L]] <<- data.frame(
        version=version,scope=tag,response=k,formula=deparse(f),term=term,n=nrow(m),df=r[term,'Df'],R2=r[term,'R2'],pseudo_F=r[term,'F'],P=r[term,'Pr(>F)'],
        status='historical_free_permutation_audit_only_mixed_support_and_repeated_sites')
    }
    g<-if(tag=='zbyonly21_cat414')m$basin else m$lake
    bd<-anova(betadisper(vegdist(mat,method='bray'),g))
    dispersion_rows[[length(dispersion_rows)+1L]] <<- data.frame(version=version,scope=tag,response=k,F=bd[1,'F value'],P=bd[1,'Pr(>F)'])
  }
}
run('crosslake26_cat414',list(~lake,~n_detected_mags+lake))
run('zbyonly21_cat414',list(~basin+season,~n_detected_mags+basin+season))
write_tsv(bind_rows(crossrows),file.path(target,'results/crosslake_PERMANOVA_historical_audit.tsv'))
write_tsv(bind_rows(dispersion_rows),file.path(target,'results/crosslake_dispersion_historical_audit.tsv'))

# Original nonparametric per-module cross-lake comparison, 49-test BH family.
d<-read_tsv(file.path(inputs,'crosslake26_cat414_response_graded_weighted.tsv'),show_col_types=FALSE);names(d)[1]<-'sample_id'
d<-left_join(d,qc %>% select(sample_id,lake),by='sample_id')
mods<-setdiff(names(d),c('sample_id','lake'))
res<-map_dfr(mods,function(m){
  x<-d[[m]][d$lake=='DXC'];y<-d[[m]][d$lake=='ZBY'];w<-suppressWarnings(wilcox.test(x,y,exact=FALSE))
  tibble(module=m,mean_DXC=mean(x),mean_ZBY=mean(y),diff=mean(x)-mean(y),delta=mean(outer(x,y,function(a,b)sign(a-b))),p=w$p.value)
}) %>% mutate(p_BH=p.adjust(p,'BH'),version=version,status='historical_unpaired_lake_contrast_audit_only')
write_tsv(res,file.path(target,'results/module_crosslake_univariate.tsv'))
capture.output(sessionInfo(),file=file.path(target,'logs/sessionInfo.txt'))
