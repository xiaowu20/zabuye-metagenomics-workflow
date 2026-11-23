suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
root<-commandArgs(trailingOnly=TRUE)[1]
source_dir<-file.path(root,"results/controlled_catalogue_comparisons")
stopifnot(file.exists(file.path(source_dir,"COMPLETE.txt")))
read_table<-function(name)read.delim(file.path(source_dir,paste0(name,".tsv")),check.names=FALSE)
ids<-c("binner","pooling","assembler")
labels<-c("2 vs 3 binners","Pooled vs single","MEGAHIT vs\nmetaSPAdes")
comp_label<-function(x)factor(x,levels=ids,labels=labels)
palette<-c(reference="#3F6B99",alternative="#D99521")
theme_set(theme_classic(base_size=16)+theme(axis.text=element_text(colour="#25282B"),
 legend.position="bottom",legend.text=element_text(size=14),
 axis.title=element_text(size=17),plot.margin=margin(14,20,14,16),
 strip.background=element_blank(),strip.text=element_text(face="bold")))

recovery<-read_table("recovery_summary")%>%filter(criterion=="match95_AF30")
stopifnot(nrow(recovery)==3,all(recovery$matched_reference_MAGs<=recovery$reference_MAGs),
 all(recovery$matched_alternative_MAGs<=recovery$alternative_MAGs))
a<-bind_rows(
 recovery%>%transmute(comparison,side="reference",matched=matched_reference_MAGs,total=reference_MAGs),
 recovery%>%transmute(comparison,side="alternative",matched=matched_alternative_MAGs,total=alternative_MAGs))%>%
 mutate(comparison=comp_label(comparison),percent=100*matched/total,label=paste0(matched,"/",total))
pd<-position_dodge(width=.45)
p1<-ggplot(a,aes(percent,comparison,colour=side,shape=side))+
 geom_point(position=pd,size=3)+geom_text(aes(label=label),position=pd,hjust=-.2,size=4,show.legend=FALSE)+
 scale_colour_manual(values=palette)+scale_shape_manual(values=c(reference=16,alternative=17))+
 scale_x_continuous(limits=c(0,130),breaks=c(0,25,50,75,100),expand=expansion(mult=c(.02,.02)))+
 labs(x="Representatives with\nreciprocal matches (%)",y=NULL,colour=NULL,shape=NULL)

samples<-read_table("sample_concordance")%>%mutate(comparison=comp_label(comparison))
stopifnot(nrow(samples)==42,all(table(samples$comparison)==14))
p2<-ggplot(samples,aes(comparison,abundance_Spearman))+
 geom_boxplot(width=.48,fill="#F4F6F7",colour="#687482",outlier.shape=NA)+
 geom_point(position=position_jitter(width=.12,height=0,seed=20261110),size=2,colour="#3B9189")+
 scale_y_continuous(limits=c(-1,1),breaks=c(-1,-.5,0,.5,1))+
 labs(x=NULL,y="Within-library abundance\nSpearman rho")+
 theme(axis.text.x=element_text(angle=20,hjust=1))
sample_n<-samples%>%group_by(comparison)%>%summarise(n=sum(is.finite(abundance_Spearman)),.groups="drop")
p2<-p2+geom_text(data=sample_n,aes(comparison,-.93,label=paste0("Defined n = ",n)),inherit.aes=FALSE,size=3.8)

fractions<-samples%>%select(comparison,sample_id,reference_abundance_fraction_matched,
 alternative_abundance_fraction_matched)%>%pivot_longer(ends_with("fraction_matched"),
 names_to="side",values_to="fraction")%>%mutate(side=if_else(grepl("^reference",side),"reference","alternative"))
stopifnot(all(is.finite(fractions$fraction)),all(fractions$fraction>=0),all(fractions$fraction<=1+1e-8))
p3<-ggplot(fractions,aes(comparison,100*fraction,colour=side))+
 geom_point(position=position_jitterdodge(jitter.width=.12,dodge.width=.5,seed=20261110),size=2,alpha=.85)+
 scale_colour_manual(values=palette)+scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
 labs(x=NULL,y="Relative abundance in\nmatched groups (%)",colour=NULL)+
 theme(axis.text.x=element_text(angle=20,hjust=1))

directions<-read_table("contrast_direction_summary")%>%mutate(comparison=comp_label(comparison),
 contrast=factor(contrast,levels=c("spring_minus_autumn","north_minus_south"),
 labels=c("Spring minus autumn","North minus south")),
 percent=100*fraction_same_direction,label=paste0(same_direction,"/",nonzero_contrasts_in_both))
stopifnot(nrow(directions)==6,all(directions$same_direction<=directions$nonzero_contrasts_in_both))
p4<-ggplot(directions,aes(percent,comparison,colour=contrast,shape=contrast))+
 geom_point(position=pd,size=3)+geom_text(aes(label=label),position=pd,hjust=-.2,size=4,show.legend=FALSE)+
 scale_colour_manual(values=c("#3F6B99","#3B9189"))+
 scale_x_continuous(limits=c(0,130),breaks=c(0,25,50,75,100))+
 labs(x="Matched groups with\nagreeing contrast signs (%)",y=NULL,colour=NULL,shape=NULL)

figure<-wrap_plots(p1,p2,p3,p4,ncol=2,widths=c(1,1),heights=c(1,1))+
 plot_annotation(tag_levels="a")
for(ext in c("pdf","svg","png"))ggsave(file.path(root,"figures",paste0("R1_Controlled_population_comparisons.",ext)),
 figure,width=14,height=10.5,dpi=300)
write.table(data.frame(comparison=gsub("\n"," ",as.character(samples$comparison),fixed=TRUE),sample_id=samples$sample_id,
 undefined_rho=!is.finite(samples$abundance_Spearman)),file.path(source_dir,"figure_undefined_correlations.tsv"),
 sep="\t",row.names=FALSE,quote=FALSE)
