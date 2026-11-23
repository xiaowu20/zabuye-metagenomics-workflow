suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
out <- args[2]
follow <- file.path(root, "deposition_qc/followup_20261112")
matrix_root <- if (length(args) >= 3) args[3] else file.path(follow, "fixed_reference_sensitivity")
tree_root <- if (length(args) >= 4) args[4] else file.path(follow, "phylogenetic_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ids<-c("baseline_original269","without_EXCLUDE_MAGs","without_all_action_MAGs","EXCLUDE_contigs_function_only")
retained <- vapply(ids, function(id) {
  as.integer(jsonlite::fromJSON(file.path(matrix_root, "catalogues", id, "definition.json"))$retained_MAGs)
}, integer(1))
labels <- paste(c("Original", "Omit EXCLUDE:", "Omit all action:", "Function-only:"), retained)
d<-bind_rows(lapply(ids,function(id)read.delim(file.path(matrix_root,"catalogues",id,"statistics/PERMANOVA.tsv"))))
d<-d%>%filter(model=="detected_MAG_count",response%in%c("GMP","SCA"))%>%
  mutate(catalogue=factor(catalogue,levels=rev(ids),labels=rev(labels)),
         term=factor(term,levels=c("season","basin"),labels=c("Sampling period","Basin")),
         nominal=factor(ifelse(P_exact<.05,"P < 0.05","P >= 0.05"),levels=c("P < 0.05","P >= 0.05")))
stopifnot(nrow(d)==16,all(d$bootstrap_valid==999))
theme_set(theme_classic(base_size=13)+theme(axis.text=element_text(colour="#25282B"),
  plot.title=element_text(size=14,face="bold"),strip.background=element_rect(fill="#F4F6F7",colour=NA),
  legend.position="bottom",plot.margin=margin(9,12,9,9)))
p1<-ggplot(d,aes(R2,catalogue))+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=catalogue),colour="#687482",linewidth=.6)+
  geom_point(aes(shape=nominal),size=2.7,colour="#3F6B99",stroke=.8)+
  geom_text(aes(x=.43,label=sprintf("P = %.4f",P_exact)),hjust=0,size=3.0)+
  facet_grid(response~term)+scale_x_continuous(limits=c(0,.56),breaks=c(0,.1,.2,.3,.4))+
  scale_shape_manual(values=c("P < 0.05"=16,"P >= 0.05"=1))+
  labs(x=expression(paste("Term-specific ",R^2," (95% station-bootstrap interval)")),y=NULL,shape=NULL,
    title="a  Fixed-reference ecological sensitivity")+
  theme(panel.spacing=unit(1.2,"lines"),legend.position="none")
f<-read.delim(file.path(tree_root,"PGLS_coefficients.tsv"))
trees<-c("archival_143_pruned","GTDB_masked_AA84","matched_AA84","matched_codon_NT84")
tree_labels<-c("Archived AA","GTDB masked AA","Matched AA","Codon-aligned NT")
f<-f%>%filter(term=="z_salt",startsWith(interval_method,"residual_df"))%>%
  mutate(tree=factor(tree,levels=rev(trees),labels=rev(tree_labels)),
    scenario=factor(scenario,levels=c("baseline84","omit_REVIEW_RARE_MAG83"),labels=c("Original 84","Omit REVIEW_RARE MAG: 83")))
stopifnot(nrow(f) == 8, all(is.finite(f$estimate)), all(f$conf_low <= f$conf_high))
f$y<-as.numeric(f$tree)+ifelse(f$scenario=="Original 84",.13,-.13)
colours<-c("Original 84"="#3F6B99","Omit REVIEW_RARE MAG: 83"="#D96B4A")
shapes<-c("Original 84"=1,"Omit REVIEW_RARE MAG: 83"=16)
stopifnot(setequal(levels(f$scenario),names(colours)),setequal(levels(f$scenario),names(shapes)))
p2<-ggplot(f,aes(estimate,y,colour=scenario))+
  geom_vline(xintercept=0,linetype=3,colour="#A9AFB5")+
  geom_segment(aes(x=conf_low,xend=conf_high,yend=y),linewidth=.6)+geom_point(aes(shape=scenario),size=2.3)+
  scale_colour_manual(values=colours)+scale_shape_manual(values=shapes)+
  scale_y_continuous(breaks=1:4,labels=rev(tree_labels),limits=c(.5,4.5))+
  scale_x_continuous(limits=c(-.04,.04),breaks=c(-.04,-.02,0,.02,.04))+
  labs(x="Standardised salt-composite coefficient (conditional 95% interval)",y=NULL,colour=NULL,shape=NULL,
    title="b  Halobacteriota sensitivity on pruned trees")+
  theme(legend.text=element_text(size=11))
figure<-p1/p2+plot_layout(heights=c(1.55,1))
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_FCS_inference_sensitivity.",ext)),figure,width=9,height=8.4,dpi=300)
write.table(d,file.path(out,"panel_a_source.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(f,file.path(out,"panel_b_source.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
