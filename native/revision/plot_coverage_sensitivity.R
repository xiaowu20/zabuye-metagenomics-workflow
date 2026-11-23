suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
root<-commandArgs(trailingOnly=TRUE)[1]
src<-file.path(root,"results/assembly_sensitivity");dir.create(src,recursive=TRUE,showWarnings=FALSE)
out<-file.path(root,"figures")
theme_set(theme_classic(base_size=17)+theme(axis.text=element_text(colour="#25282B"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="bottom"))
np<-bind_rows(lapply(c("read_species","archival_original269"),function(g) {
  read.delim(file.path(root,"results/coverage_sensitivity",g,"PERMANOVA_Nonpareil_sensitivity.tsv"))
}))
np<-np %>% filter(catalogue%in%c("read_species","archival_original269")) %>%
  mutate(response=factor(response,levels=c("read_species","MAG","GMP","SCA","DCP")),
    model=factor(model,levels=c("unadjusted","detected_MAG_count","Nonpareil_coverage","both"),
      labels=c("Season + basin","+ detected MAG count","+ Nonpareil coverage","+ both diagnostics")),
    scope=factor(scope,levels=c("paired14","omit_station_N3"),labels=c("Seven station pairs","Six pairs, N3 omitted")),
    term=factor(term,levels=c("season","basin")),
    plot_y=as.numeric(model)+ifelse(scope=="Seven station pairs",-.10,.10),nominal=P_exact<.05)
# Apply the same explicit vertical offset to both ends of every CI.
# Position-dodging geom_segment shifts y without consistently shifting yend.
p7<-ggplot(np,aes(R2,plot_y,colour=scope))+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=plot_y),linewidth=.5)+
  geom_point(aes(shape=nominal),size=2.6,stroke=.7)+
  scale_y_continuous(breaks=seq_len(nlevels(np$model)),labels=levels(np$model),limits=c(.65,4.35),expand=expansion(mult=0))+
  facet_grid(response~term,scales="free_x",labeller=labeller(
    response=c(read_species="Species",MAG="MAG",GMP="GMP",SCA="SCA",DCP="DCP"),
    term=c(season="Season",basin="Basin")))+
  scale_colour_manual(values=c("#3F6B99","#D96B4A"))+
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),labels=c(`TRUE`="P < 0.05",`FALSE`="P >= 0.05"))+
  guides(colour=guide_legend(order=1,nrow=1),shape=guide_legend(order=2,nrow=1))+
  labs(x="Term-specific R-squared (95% station bootstrap CI)",y=NULL,colour=NULL,shape="Exact nominal test")+
  theme(legend.box="vertical")
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_Coverage_model_sensitivity.",ext)),p7,width=11,height=12.5,dpi=300)
write.table(np,file.path(src,"coverage_sensitivity_figure_data.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
