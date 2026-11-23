suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
root <- commandArgs(trailingOnly=TRUE)[1]
out <- file.path(root,"results/nonpareil")
stats <- read.delim(file.path(out,"sample_coverage_summary.tsv"))
pred <- read.delim(file.path(out,"fitted_coverage_curves.tsv"))
raw <- read.delim(file.path(out,"rarefied_coverage_points.tsv"))
pal<-c(autumn="#7A63A8",spring="#55A47B")
theme_set(theme_classic(base_size=13)+theme(axis.text=element_text(colour="#25282B"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="bottom"))
pred$season<-factor(pred$season,levels=c("autumn","spring"))
stats$season<-factor(stats$season,levels=c("autumn","spring"))
pred$station<-factor(pred$station,levels=c("N1","N2","N3","N4","S1","S2","S3"))
stats$station<-factor(stats$station,levels=levels(pred$station))
p1<-ggplot(pred,aes(effort_bp/1e9,coverage*100,colour=season,group=sample))+
  geom_hline(yintercept=95,linetype=3,colour="#A9AFB5",linewidth=.35)+
  geom_line(aes(linetype=range,group=interaction(sample,range)),linewidth=.55)+
  geom_point(data=stats,aes(x=sequencing_effort_bp/1e9,y=observed_average_coverage*100),size=1.9)+
  facet_wrap(~station,nrow=2)+scale_x_log10(breaks=c(.001,.1,10),labels=c("0.001","0.1","10"))+
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+
  scale_colour_manual(values=pal)+scale_linetype_manual(values=c("Within input effort"="solid","Extrapolation"="dashed"))+
  labs(x="Sequencing effort from R1 (Gbp)",y="Nonpareil average coverage (%)",colour=NULL,linetype=NULL)
p2<-ggplot(stats,aes(season,observed_average_coverage*100,group=station))+
  geom_line(colour="#A9AFB5",linewidth=.45)+geom_point(aes(colour=season,shape=basin),size=2.7)+
  scale_colour_manual(values=pal)+labs(x=NULL,y="Coverage at input (%)",colour=NULL,shape="Basin")+
  guides(colour="none",shape=guide_legend(override.aes=list(colour="#25282B")))
p3<-ggplot(stats,aes(season,Nonpareil_diversity_Nd,group=station))+
  geom_line(colour="#A9AFB5",linewidth=.45)+geom_point(aes(colour=season,shape=basin),size=2.7)+
  scale_colour_manual(values=pal)+labs(x=NULL,y="Sequence diversity (Nd)",colour=NULL,shape="Basin")+
  guides(colour="none",shape=guide_legend(override.aes=list(colour="#25282B")))
p<-p1/(p2+p3)+plot_layout(heights=c(2,1.15),guides="collect")+plot_annotation(tag_levels="a")
fig<-file.path(root,"figures");dir.create(fig,showWarnings=FALSE)
ggsave(file.path(fig,"R1_Nonpareil14.pdf"),p,width=10,height=8.2)
ggsave(file.path(fig,"R1_Nonpareil14.svg"),p,width=10,height=8.2)
ggsave(file.path(fig,"R1_Nonpareil14.png"),p,width=10,height=8.2,dpi=300)
