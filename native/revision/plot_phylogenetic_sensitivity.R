suppressPackageStartupMessages({library(ape);library(phangorn);library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
root<-commandArgs(trailingOnly=TRUE)[1]
src<-file.path(root,"results/phylogenetic_sensitivity");out<-file.path(root,"figures")
fit<-read.delim(file.path(src,"PGLS_coefficients_by_tree.tsv"))
sig<-read.delim(file.path(src,"Pagel_lambda_by_tree.tsv"))
support<-read.delim(file.path(src,"inferred_branch_supports.tsv"))
cmp<-read.delim(file.path(src,"topology_comparisons.tsv"))
ids<-c("archival_143_pruned","GTDB_masked_AA84","matched_AA84","matched_codon_NT84")
labs<-c("Archived AA","GTDB masked AA","Matched AA","Codon-aligned NT")
stopifnot(setequal(unique(fit$tree),ids))
f<-fit %>% filter(term=="z_salt",interval_method=="archival_ML_scale_conditional_t79") %>%
  mutate(tree=factor(tree,levels=rev(ids),labels=rev(labs)))
theme_set(theme_classic(base_size=16)+theme(axis.text=element_text(colour="#25282B"),
  plot.title=element_text(face="bold"),legend.position="bottom",
  plot.margin=margin(10,14,10,10)))
p1<-ggplot(f,aes(estimate,tree))+
  geom_vline(xintercept=0,linetype=3,colour="#A9AFB5")+
  geom_segment(aes(x=conf_low,xend=conf_high,yend=tree),colour="#687482",linewidth=.7)+
  geom_point(size=2.5,colour="#3F6B99")+
  labs(x="Standardised salt-composite effect\n(95% CI)",y=NULL)

# Stacked-bar component adapted from the project's archived example labelled
# PlotCase 多层注释柱状图: pivot_longer + stacked geom_col + percent scale.
# The example was an existing project input; no new PlotCase retrieval occurred.
df<-support %>% filter(tree!="archival_143_pruned") %>%
  mutate(category=ifelse(UFBoot>=95 & SH_aLRT>=80,"Both thresholds met","Below either threshold")) %>%
  count(tree,category) %>% group_by(tree) %>% mutate(value=n/sum(n)) %>% ungroup() %>%
  select(tree,category,value) %>% pivot_wider(names_from=category,values_from=value,values_fill=0)
p2<-df %>% pivot_longer(-tree) %>%
  mutate(tree=factor(tree,levels=ids[-1],labels=labs[-1])) %>%
  ggplot(aes(tree,value,fill=name))+
  geom_col(position="stack",width=.65)+scale_fill_manual(values=c("Both thresholds met"="#3F6B99","Below either threshold"="#DCE1E5"),
    labels=c("Both thresholds met"="Both thresholds met","Below either threshold"="At least one below"))+
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1))+
  coord_flip()+labs(x=NULL,y="Internal branches (81 per tree)",fill=NULL)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

rf<-bind_rows(cmp %>% transmute(tree1,tree2,value=normalized_RF),
  data.frame(tree1=cmp$tree2,tree2=cmp$tree1,value=cmp$normalized_RF),data.frame(tree1=ids,tree2=ids,value=0)) %>%
  mutate(tree1=factor(tree1,levels=ids,labels=c("Archived","GTDB AA","Matched AA","Nucleotide")),
    tree2=factor(tree2,levels=rev(ids),labels=rev(c("Archived","GTDB AA","Matched AA","Nucleotide"))))
stopifnot(nrow(rf)==16,!anyDuplicated(paste(rf$tree1,rf$tree2)))
p3<-ggplot(rf,aes(tree1,tree2,fill=value))+geom_tile(colour="white",linewidth=.7)+
  geom_text(aes(label=sprintf("%.3f",value)),size=4.5)+
  scale_fill_gradient(low="#F4F6F7",high="#E5A442",limits=c(0,max(rf$value)),breaks=c(0,.08,.16))+
  coord_equal()+labs(x=NULL,y=NULL,fill="Normalised RF")+
  guides(fill=guide_colourbar(barwidth=unit(3,"cm"),barheight=unit(.35,"cm"),title.position="top"))+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=35,hjust=1))
s<-sig %>% mutate(tree=factor(tree,levels=rev(ids),labels=rev(labs)))
p4<-ggplot(s,aes(lambda,tree))+geom_point(size=2.5,colour="#3B9189")+
  scale_x_continuous(limits=c(0,1.01),breaks=c(0,.5,1))+
  labs(x="Pagel's lambda for median pI",y=NULL)
p<-wrap_plots(p1,p2,p3,p4,ncol=2,widths=c(1,1),heights=c(1,1))+plot_annotation(tag_levels="a")
for(ext in c("pdf","svg","png"))ggsave(file.path(out,paste0("R1_Phylogenetic_sensitivity.",ext)),p,width=12,height=10,dpi=300)

obj<-readRDS(file.path(src,"comparison.rds"))
key<-obj$data[,c("mag_id","gtdb_genus")]
key$display_label<-paste(gsub("ZBY_water__","",key$mag_id),key$gtdb_genus,sep=" | ")
write.table(key,file.path(src,"tree_tip_key.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
for(id in ids[-1]) {
  tr<-ladderize(obj$trees[[id]])
  tr$tip.label<-key$display_label[match(tr$tip.label,key$mag_id)]
  for(ext in c("pdf","svg","png")) {
    path<-file.path(out,paste0("R1_tree_",id,".",ext))
    if(ext=="pdf")pdf(path,width=15,height=20,useDingbats=FALSE)
    if(ext=="svg")svglite::svglite(path,width=15,height=20)
    if(ext=="png")png(path,width=3000,height=4000,res=200)
    par(mar=c(5,1,3,1))
    plot.phylo(tr,cex=.68,no.margin=FALSE,edge.width=.65)
    nodelabels(text=tr$node.label,frame="none",cex=.40,adj=c(1,-.3),col="#687482")
    add.scale.bar(cex=.8)
    title(main=labs[match(id,ids)])
    mtext("Node labels: SH-aLRT / ultrafast bootstrap (%); midpoint-rooted display",side=1,line=2,cex=.85)
    mtext("Branch-length scale: expected substitutions per site",side=1,line=3.2,cex=.85)
    dev.off()
  }
}
