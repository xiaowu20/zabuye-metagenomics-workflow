suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
root<-commandArgs(trailingOnly=TRUE)[1];src<-file.path(root,"results/assembly_sensitivity");out<-file.path(root,"figures")
ids<-c("archival_original269","remapped_original269","original_ZBY_water_das2_matched","original_ZBY_water_das3","pooled_megahit14","single_megahit14","single_metaspades14")
labs<-c("Original 269 / 414 refs","Original 269 / 269 refs","Archived 21 / 2 binners","Archived 21 / 3 binners","Pooled 14 / MEGAHIT","Single 14 / MEGAHIT","Single 14 / metaSPAdes")
label_catalogue<-function(x)factor(x,levels=rev(ids),labels=rev(labs))
read<-function(name)read.delim(file.path(src,paste0(name,".tsv")),check.names=FALSE)
theme_set(theme_classic(base_size=17)+theme(axis.text=element_text(colour="#25282B"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="bottom"))

# Existing project example labelled PlotCase 多层注释柱状图 supplies the
# pivot_longer and stacked geom_col component; all scientific summaries are new.
recovery<-read("catalogue_recovery") %>% filter(catalogue!="archival_original269")
p1<-recovery %>% select(catalogue,HQ,MQ) %>% pivot_longer(-catalogue,names_to="tier",values_to="n") %>%
  ggplot(aes(label_catalogue(catalogue),n,fill=tier))+geom_col(position="stack",width=.65)+
  scale_fill_manual(values=c(HQ="#3F6B99",MQ="#DCE1E5"))+
  scale_y_continuous(expand=expansion(mult=c(0,.05)))+coord_flip()+
  labs(x=NULL,y="95% ANI catalogue representatives",fill="Project quality tier")
ani<-read("ANI_recovery_summary") %>% mutate(criterion=factor(criterion,
  levels=c("match95_AF30","match99_AF30","match95_AF65"),labels=c("95% ANI, AF >= 0.30","99% ANI, AF >= 0.30","95% ANI, AF >= 0.65")))
p2<-ggplot(ani,aes(matched_original_MAGs/269*100,label_catalogue(catalogue),colour=criterion,shape=criterion))+
  geom_point(position=position_dodge(width=.45),size=2.7)+
  scale_colour_manual(values=c("#3F6B99","#3B9189","#D96B4A"))+
  scale_x_continuous(limits=c(0,100))+guides(colour=guide_legend(ncol=1),shape=guide_legend(ncol=1))+
  labs(x="Original MAGs with\nmatches (%)",y=NULL,colour=NULL,shape=NULL)
qc<-read("sample_coverage") %>% filter(breadth_cutoff==.5,catalogue!="archival_original269")
p3<-ggplot(qc,aes(label_catalogue(catalogue),prethreshold_catalogue_capture_pct))+
  geom_boxplot(width=.5,outlier.shape=NA,colour="#687482",fill="#F4F6F7")+
  geom_point(position=position_jitter(width=.13,height=0,seed=20261110),size=1.7,colour="#3B9189",alpha=.8)+
  coord_flip()+labs(x=NULL,y="Unthresholded catalogue capture (%)")
p4<-ggplot(qc,aes(label_catalogue(catalogue),n_detected_mags))+
  geom_boxplot(width=.5,outlier.shape=NA,colour="#687482",fill="#F4F6F7")+
  geom_point(position=position_jitter(width=.13,height=0,seed=20261110),size=1.7,colour="#D99521",alpha=.8)+
  coord_flip()+labs(x=NULL,y="Detected MAGs per library")
p<-(p1+p2)/(p3+p4)+plot_annotation(tag_levels="a")
p<-p & theme(text=element_text(size=20),axis.title=element_text(size=18),
  axis.text=element_text(size=17.6),legend.text=element_text(size=16),
  legend.title=element_text(size=17),plot.margin=margin(12,16,12,12))
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_Catalogue_recovery.",ext)),p,width=14,height=9.5,dpi=300)

stats<-read("PERMANOVA_all_catalogues") %>% filter(breadth==.5,technical_adjustment=="detected_MAG_count") %>%
  mutate(catalogue=label_catalogue(catalogue),response=factor(response,levels=c("MAG","GMP","SCA","DCP")),
    term=factor(term,levels=c("season","basin")),nominal=P_exact<.05)
p5<-ggplot(stats,aes(R2,catalogue))+
  geom_segment(aes(x=R2_CI_low,xend=R2_CI_high,yend=catalogue),linewidth=.55,colour="#A9AFB5")+
  geom_point(aes(fill=nominal),shape=21,size=2.4,colour="#3F6B99")+
  scale_fill_manual(values=c(`TRUE`="#3F6B99",`FALSE`="white"),labels=c(`TRUE`="P < 0.05",`FALSE`="P >= 0.05"))+
  facet_grid(response~term,scales="free_x")+labs(x="Term-specific R-squared (95% station bootstrap CI)",y=NULL,fill="Exact nominal test")
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_Assembly_ecological_sensitivity.",ext)),p5,width=11,height=11.5,dpi=300)

# Grid and filled-circle component adapted from the archived project example
# labelled PlotCase 多分组气泡热图. No case discovery or application-asset reading.
effects<-read("module_basin_contrasts") %>% filter(response=="GMP",breadth==.5)
top<-effects %>% filter(catalogue=="archival_original269") %>% arrange(desc(abs(north_minus_south))) %>% slice_head(n=12) %>% pull(module)
# Display names retained from the final Figure 3 source; identifiers stay in tables.
module_labels<-c(Salt_Trk_Ktr_core="Trk/Ktr marker (K03498)",
  TCA_step8_malate_dh="Malate dehydrogenase",C1_formate_ox="Formate dehydrogenase",
  N_urea_utilization="Urea hydrolysis (ureABC)",Salt_betaine_synthesis="Betaine synthesis (betAB)",
  N_nitrate_red_periplasmic="NapAB nitrate reduction",
  TCA_step5_succinylCoA_synthetase_panelSucCD="Succinyl-CoA synthetase",
  N_nitrate_red_membrane="NarGHI-like module",TCA_step2_aconitase="Aconitase",
  Salt_hydroxyectoine_ext="Hydroxyectoine (ectD)",TCA_step1_citrate_synthase="Citrate synthase",
  Salt_betaine_transport="Choline/betaine transport")
stopifnot(all(top%in%names(module_labels)))
e<-effects %>% filter(module%in%top) %>% mutate(catalogue=factor(catalogue,levels=ids,labels=labs),
  module=factor(module,levels=rev(top)),difference_pp=100*north_minus_south,magnitude=abs(difference_pp))
p6<-ggplot(e,aes(catalogue,module,size=magnitude,fill=difference_pp))+
  geom_tile(colour="#DCE1E5",fill="white",linewidth=.3,show.legend=FALSE)+
  geom_point(shape=21,colour="#25282B",stroke=.25)+
  scale_size_area(max_size=7,breaks=c(10,30,50),name="Absolute difference (pp)")+
  scale_fill_gradient2(low="#E5A442",mid="white",high="#3F6B99",midpoint=0,
    name="North minus south (pp)",n.breaks=5)+
  scale_y_discrete(labels=module_labels)+
  scale_x_discrete(labels=function(x)gsub(" / ","\n",x,fixed=TRUE))+
  guides(size=guide_legend(order=1,title.position="top",nrow=1),
    fill=guide_colourbar(order=2,title.position="top",barwidth=grid::unit(6,"cm"),
      barheight=grid::unit(.45,"cm")))+
  labs(x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=35,hjust=1,size=14),
    axis.text.y=element_text(size=15),axis.ticks=element_blank(),axis.line=element_blank(),
    legend.box="vertical",legend.justification="center",legend.title=element_text(size=15),
    legend.text=element_text(size=14),legend.spacing.y=grid::unit(.15,"cm"),
    plot.margin=margin(16,20,16,16))
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_Module_direction_sensitivity.",ext)),p6,width=12,height=8.5,dpi=300)

# The figure runner invokes the separate controlled-comparisons renderer for Figure S26.
