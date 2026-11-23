suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
root<-commandArgs(trailingOnly=TRUE)[1];follow<-file.path(root,"deposition_qc/followup_20261112")
out <- commandArgs(trailingOnly = TRUE)[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
q<-read.delim(file.path(follow,"tables/paired_CheckM2_comparison_6.tsv"))
m<-read.delim(file.path(follow,"tables/paired_candidate_manifest_6.tsv"))
pi<-read.delim(file.path(follow,"scientific_review/EXCLUDE_contig_protein_pI_sensitivity.tsv"))
ids <- sort(unique(as.character(m$mag_id)), method = "radix")
stopifnot(!anyDuplicated(q$mag_id), !anyDuplicated(m$mag_id), setequal(q$mag_id, ids))
pattern <- "^(.+)__(metabat2|maxbin2)_bin[.]([0-9]+.*)$"
parts <- regmatches(ids, regexec(pattern, ids))
stopifnot(all(lengths(parts) == 4L))
labels <- vapply(parts, function(x) {
  site <- sub("_water$", " water", sub("_sed$", " sed.", x[2]))
  paste(site, paste0(c(metabat2 = "MB", maxbin2 = "MX")[[x[3]]], x[4]))
}, character(1))
q$label<-factor(q$mag_id,levels=rev(ids),labels=rev(labels));m$label<-factor(m$mag_id,levels=rev(ids),labels=rev(labels))
colours<-c("Source control"="#3F6B99","EXCLUDE-only candidate"="#D96B4A")
theme_set(theme_classic(base_size=16)+theme(axis.text=element_text(colour="#25282B"),axis.title=element_text(size=13),
  plot.title=element_text(size=14,face="bold"),plot.subtitle=element_text(size=12),
  legend.position="bottom",plot.margin=margin(7,7,7,7)))
# Adapt the archived PlotCase bar example's pivot_longer/stacked geom_col component.
# Keep percent denominators explicit; this plot does not estimate biological contamination.
bars<-m%>%transmute(mag_id,label,Retained=candidate_bp/original_bp,EXCLUDE=removed_bp/original_bp)%>%
  pivot_longer(c(Retained,EXCLUDE),names_to="sequence_status",values_to="fraction")
p1<-ggplot(bars,aes(label,fraction,fill=sequence_status))+geom_col(width=.62,position="stack")+
  geom_text(data=m,inherit.aes=FALSE,mapping=aes(x=label,label=ifelse(100*removed_bp/original_bp<.1,sprintf("%.3f%%",100*removed_bp/original_bp),sprintf("%.2f%%",100*removed_bp/original_bp)),y=1.035),
    hjust=0,size=4.5)+
  scale_fill_manual(values=c(Retained="#DCE1E5",EXCLUDE="#D96B4A"))+
  scale_y_continuous(labels=scales::percent,breaks=c(0,.5,1),limits=c(0,1.45),expand=c(0,0))+
  coord_flip()+labs(x=NULL,y="Original genome fraction",fill=NULL,title="a  EXCLUDE contigs")+
  guides(fill="none")
paired<-function(field,title,xlab) {
  values<-q%>%select(mag_id,label,all_of(paste0(c("paired_source_","exclude_only_"),field)))%>%
    pivot_longer(-c(mag_id,label),names_to="stage",values_to="value")%>%
    mutate(stage=ifelse(startsWith(stage,"paired_source_"),"Source control","EXCLUDE-only candidate"))
  ggplot(values,aes(value,label))+geom_segment(data=q,aes(x=.data[[paste0("paired_source_",field)]],xend=.data[[paste0("exclude_only_",field)]],yend=label),colour="#A9AFB5",linewidth=.5)+
    geom_point(aes(colour=stage,shape=stage),size=2.5,stroke=.65)+
    scale_colour_manual(values=colours)+scale_shape_manual(values=c("Source control"=1,"EXCLUDE-only candidate"=16))+
    labs(x=xlab,y=NULL,title=title,colour=NULL,shape=NULL)
}
p2<-paired("Completeness","b  Completeness","Completeness (%)")+
  geom_vline(xintercept=90,linetype=3,colour="#A9AFB5")+scale_x_continuous(limits=c(89.5,96))
p3<-paired("Contamination","c  Contamination estimate","Contamination estimate (%)")+
  geom_vline(xintercept=5,linetype=3,colour="#A9AFB5")+scale_x_continuous(limits=c(0,5.2))
water <- ids[startsWith(ids, "ZBY_water__")]
pi <- pi[pi$mag_id %in% water, ]
stopifnot(setequal(pi$mag_id, water))
pi$label <- factor(pi$mag_id, levels = rev(water), labels = rev(labels[match(water, ids)]))
pi$delta<-pi$candidate_median_pI-pi$source_median_pI
p4<-ggplot(pi,aes(delta,label))+geom_vline(xintercept=0,linetype=3,colour="#A9AFB5")+
  geom_segment(aes(x=0,xend=delta,yend=label),colour="#A9AFB5",linewidth=.55)+geom_point(size=2.5,colour="#D96B4A")+
  labs(x="Change in predicted\nmedian protein pI",y=NULL,title="d  Predicted protein pI",subtitle="Candidate minus source\n(ZBY water only)")+
  scale_x_continuous(limits=c(-.005,.020),breaks=c(0,.01,.02))
figure<-(p1+p2)/(p3+p4)+plot_layout(guides="collect")&theme(legend.position="bottom")
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_FCS_quality_followup.",ext)),figure,width=9,height=6.6,dpi=300)
write.table(data.frame(mag_id=ids,display_label=labels),file.path(out,"figure_MAG_key.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(bars,file.path(out,"panel_a_sequence_fractions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(q,file.path(out,"panels_bc_paired_quality.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(pi,file.path(out,"panel_d_protein_pI.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
