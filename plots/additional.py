"""Data-driven redraws for the genome traits, environmental models and carriers."""
import argparse
from functools import lru_cache
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import linregress, t

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("books",type=Path)
parser.add_argument("local_inputs",type=Path)
parser.add_argument("output",type=Path)
parser.add_argument("--paper",action="store_true",help="Render the seven panels whose numerical inputs are distributed in the supplements")
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})
red,blue="#ae513b","#3c6c91"

@lru_cache(maxsize=None)
def sheet(name):
    return pd.read_excel(args.books/"Supplementary_Data_S4.xlsx",sheet_name=name)

def extra(name):
    return pd.read_csv(args.local_inputs/"extra"/(name+".tsv"),sep="\t")

def save(fig, name):
    for extension in ("png","pdf","svg"):
        fig.savefig(args.output/f"Figure_{name}.{extension}",dpi=220,bbox_inches="tight")
    plt.close(fig)

def density(ax,tab):
    for (group,g),colour in zip(tab.groupby("group",sort=True),(red,blue)):
        g=g.sort_values("x")
        ax.plot(g.x,g.y,label=group,color=colour)
        if {"lo","hi"}<=set(g):ax.fill_between(g.x,g.lo,g.hi,color=colour,alpha=.15)
    ax.set(xlabel="Predicted protein pI",ylabel="Equal-MAG-weight density");ax.legend(fontsize=7)

def forest(ax,rows,labels,title):
    y=np.arange(len(rows))
    ax.axvline(0,ls=":",lw=.8,color="gray")
    ax.hlines(y,rows.conf_low,rows.conf_high,color=blue)
    ax.scatter(rows.estimate,y,s=15,color=blue)
    ax.set_yticks(y,labels,fontsize=7);ax.invert_yaxis()
    ax.set_title(title,loc="left");ax.set_xlabel("Coefficient (95% CI)")

def fitline(ax,x,y,ci=True):
    x=np.asarray(x,float);y=np.asarray(y,float);fit=linregress(x,y)
    xp=np.linspace(x.min(),x.max(),100);yp=fit.intercept+fit.slope*xp
    ax.plot(xp,yp,color=blue,lw=1)
    if ci:
        se=np.sqrt(np.sum((y-fit.intercept-fit.slope*x)**2)/(len(y)-2))
        width=t.ppf(.975,len(y)-2)*se*np.sqrt(1/len(y)+(xp-x.mean())**2/np.sum((x-x.mean())**2))
        ax.fill_between(xp,yp-width,yp+width,color=blue,alpha=.16)

def tree_panel(ax):
    tree=Phylo.read(args.local_inputs/"main_figure_inputs/archival_midpoint_halo84.nwk","newick")
    tips=tree.get_terminals();depth=tree.depths()
    angles={node:i*2*np.pi/len(tips) for i,node in enumerate(tips)}
    radius=max(depth.values())
    for node in tree.get_nonterminals(order="postorder"):
        child_angles=[angles[c] for c in node.clades];angles[node]=np.mean(child_angles)
        ts=np.linspace(min(child_angles),max(child_angles),60)
        ax.plot(ts,np.full_like(ts,depth[node]/radius),color="#64686a",lw=.5)
        for child in node.clades:
            ax.plot([angles[child]]*2,[depth[node]/radius,depth[child]/radius],color="#64686a",lw=.5)
    trait=sheet("pI_model_input_269").set_index("mag_id").loc[[t.name for t in tips]]
    markers=sheet("Marker_pI_status").set_index("mag_id").loc[trait.index]
    genus_codes=pd.Categorical(trait.gtdb_genus).codes
    theta=[angles[t] for t in tips]; width=2*np.pi/len(tips)
    rings=[np.where(trait.quality_tier.astype(str).str.startswith("HQ"),blue,"#ddd"),
           plt.cm.tab20((genus_codes%20)/19),plt.cm.viridis(Normalize()(trait.median_pI)),
           np.where(markers.has498.eq(True),red,"#eee"),np.where(markers.has499.eq(True),blue,"#eee")]
    for i,colors in enumerate(rings):ax.bar(theta,[.065]*len(tips),width=width,bottom=1.02+.07*i,color=colors,linewidth=0)
    comp=trait.salt_composite.to_numpy()
    ax.bar(theta,.20*comp/max(1,comp.max()),width=width,bottom=1.39,color="#71747b",linewidth=0)
    ax.set_axis_off();ax.set_title("a  Organismal tree and genome traits",fontsize=9)
    ax.text(.5,-.04,"Rings: quality, genus, pI, K03498, K03499. Bars: salt score.",
            transform=ax.transAxes,ha="center",fontsize=6)

def protein_models():
    models = sheet("Salt_composite_models").query("term=='z_salt'").copy()
    modules = sheet("Module_pI_models_10").query("term=='z_module'")
    return models, modules

def density_figures():
    tab=sheet("ZBY_proteome_pI")
    halo=tab[tab.gtdb_phylum=="Halobacteriota"]
    other=tab[tab.gtdb_phylum!="Halobacteriota"]
    tests=sheet("pI_tests").set_index("subset")
    density_all=extra("all_MAG_density")
    fig=plt.figure(figsize=(12,9),layout="constrained")
    gs=fig.add_gridspec(2,6)
    tree_panel(fig.add_subplot(gs[0,:2],projection="polar"))
    ax=fig.add_subplot(gs[0,2:4]);density(ax,density_all);ax.set_title("b  Predicted proteome acidity")
    inset=ax.inset_axes([.63,.45,.33,.38])
    vals=[halo.median_pI,other.median_pI];inset.violinplot(vals,showextrema=False)
    rng=np.random.default_rng(13)
    for i,(v,c) in enumerate(zip(vals,(red,blue)),1):
        inset.scatter(i+rng.uniform(-.12,.12,len(v)),v,s=2,alpha=.45,color=c)
        inset.scatter(i,v.median(),marker="D",s=13,color="black")
    inset.set_xticks([1,2],[f"Halo.\n{len(halo)} MAGs",f"Other\n{len(other)} MAGs"],fontsize=6);inset.tick_params(labelsize=6)
    inset.set_title(f"Medians {halo.median_pI.median():.2f} / {other.median_pI.median():.2f}\nWilcoxon P = {tests.at['Halobacteriota vs other phyla','p_value']:.3g}",fontsize=6)
    ax=fig.add_subplot(gs[0,4:])
    for tier,col in (("MQ","#b8a4d8"),("HQ","#554475")):
        h=halo[halo.quality_tier==tier]
        ax.scatter(h.acidic_minus_basic,h.median_pI,s=12,color=col,label=tier)
    ax.legend(loc="lower left",fontsize=7)
    fitline(ax,halo.acidic_minus_basic,halo.median_pI)
    balance=tests.loc["Halobacteriota pI vs acidic-minus-basic balance"]
    ax.text(.02,.97,f"Spearman ρ = {balance.statistic:.3f}, P = {balance.p_value:.3g}",va="top",transform=ax.transAxes,fontsize=7)
    ax.set(xlabel="Acidic minus basic residues (percentage points)",ylabel="MAG-median pI",title="c  Residue balance")
    m, modules = protein_models()
    pgls=sheet("Halobacteriota_PGLS").query("term=='z_salt'").assign(model="Halobacteriota PGLS")
    allmodels=pd.concat([m,pgls],ignore_index=True)
    signal=sheet("Phylogenetic_signal").iloc[0]
    ax=fig.add_subplot(gs[1,:3])
    forest(ax,allmodels,["Unadjusted","Phylum contrast + quality","Phylum intercept + quality","Additional N50","Within-lineage PGLS"],f"d  Salt-module composite\nWithin-lineage λ = {signal['lambda']:.4f}, PGLS P = {pgls.p_value.iloc[0]:.3f}")
    forest(fig.add_subplot(gs[1,3:]),modules,modules.module.str.replace("Salt_","",regex=False),f"e  Individual module estimates\nMinimum BH q = {modules.p_BH.min():.3f}")
    save(fig,"4")
    fig,ax=plt.subplots();density(ax,extra("HQ_density"))
    hq=tests.loc["HQ-only Halobacteriota vs other phyla"]
    ax.set_title(f"High-quality subset ({sum(halo.quality_tier=='HQ')} Halobacteriota and {sum(other.quality_tier=='HQ')} other MAGs)\nWilcoxon W = {hq.statistic:g}, P = {hq.p_value:.3g}")
    save(fig,"S9")

def traits():
    tab=sheet("ZBY_proteome_pI")
    rng=np.random.default_rng(13)
    m, modules = protein_models()
    fig,ax=plt.subplots(figsize=(7,5));groups=[(name,g) for name,g in tab.groupby("gtdb_phylum") if len(g)>=3]
    groups.sort(key=lambda x:x[1].median_pI.median())
    for i,(name,g) in enumerate(groups):
        ax.boxplot(g.median_pI,positions=[i],vert=False,showfliers=False)
        ax.scatter(g.median_pI,i+rng.uniform(-.15,.15,len(g)),s=5,color=red if name=="Halobacteriota" else blue)
    ax.set_yticks(range(len(groups)),[f"{name} ({len(g)})" for name,g in groups]);ax.set_xlabel("MAG-median predicted pI")
    ax.axvline(4.9,color=red,ls="--",lw=.7,label="Halobacterium NRC-1 reference")
    ax.legend(fontsize=7)
    save(fig,"S7")
    fig,axs=plt.subplots(2,2,figsize=(10,7),layout="constrained")
    hits=sheet("Confirmed_protein_hits").groupby(["KO","mag_id"]).pI_delta.median().reset_index()
    kos=list(hits.KO.unique());axs[0,0].boxplot([hits.loc[hits.KO==k,"pI_delta"] for k in kos],tick_labels=kos)
    axs[0,0].axhline(0,lw=.8,ls=":",color="gray");axs[0,0].set_ylabel("Within-MAG marker pI offset")
    forest(axs[0,1],m,["Unadjusted","Phylum contrast","Phylum intercept","Additional N50"],"Salt composite")
    forest(axs[1,0],modules,modules.module.str.replace("Salt_","",regex=False),"Individual modules")
    h=sheet("pI_model_input_269").query("gtdb_phylum=='Halobacteriota'")
    for name,g in h.groupby("gtdb_genus"):axs[1,1].scatter(g.salt_composite,g.median_pI,s=10,label=name)
    fitline(axs[1,1],h.salt_composite,h.median_pI);axs[1,1].set(xlabel="Salt-module composite",ylabel="MAG-median pI")
    save(fig,"S10")

def lakeplot(ax):
    lakes=sheet("QXLSG_lake_means")
    ax.scatter(lakes.log10_salinity,lakes.mean_MAG_median_pI,s=np.sqrt(lakes.n_attributable_MAGs)*3,color=blue)
    fitline(ax,lakes.log10_salinity,lakes.mean_MAG_median_pI)
    ax.set(xlabel="log₁₀ lake salinity (%)",ylabel="Lake mean of MAG-median pI")

def other():
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout="constrained");lakeplot(axs[0])
    groups=sheet("QXLSG_proteome_catalogue").groupby("Phylum").median_pI.agg(["median","size"]).query("size>=5").sort_values("median")
    axs[1].scatter(groups["median"],np.arange(len(groups)),s=10,color=[red if ph=="Halobacteriota" else blue for ph in groups.index])
    axs[1].set_yticks(range(len(groups)),groups.index,fontsize=4);axs[1].set_xlabel("Phylum median pI")
    records=sheet("Fig5d_MAG_records");arr=list(records.groupby(["salinity_class","lineage"],sort=True))
    for i,(_,g) in enumerate(arr):
        if len(g)>2:axs[2].boxplot(g.median_pI,positions=[i],showfliers=False)
        else:axs[2].scatter(np.full(len(g),i),g.median_pI)
    axs[2].set_xticks(range(len(arr)),["\n".join(map(str,k)) for k,_ in arr],rotation=25,fontsize=7)
    z=sheet("ZBY_proteome_pI");h=z.gtdb_phylum=="Halobacteriota"
    for v,c in ((z.loc[h,"median_pI"],red),(z.loc[~h,"median_pI"],blue)):axs[2].axhline(v.median(),ls="--",lw=.7,color=c)
    save(fig,"S12")
    fig,axs=plt.subplots(1,2,figsize=(9,4),layout="constrained");lakeplot(axs[0])
    cf=sheet("QXLSG_model_coefficients");cf=cf[cf.term.astype(str).str.contains("sal",case=False)]
    forest(axs[1],cf,cf.model,"Salinity slopes");save(fig,"S13")
    fig,axs=plt.subplots(1,2,figsize=(9,4),layout="constrained")
    models=sheet("DOM_rebuilt_models");labels=models.version.astype(str)+" / "+models.method.astype(str)
    axs[0].barh(range(len(models)),models.adjusted_R2,color=blue);axs[0].set_yticks(range(len(models)),labels,fontsize=6)
    axs[0].set_xlabel("Adjusted R²")
    for col in models.filter(regex="^P_").columns:axs[1].scatter(models[col],range(len(models)),label=col[2:],s=20)
    axs[1].axvline(.05,ls="--",lw=.7,color="gray")
    axs[1].set_yticks([]);axs[1].set_xlabel("Permutation P");axs[1].legend(fontsize=6);save(fig,"S15")
    from environment_carriers import render
    render(args.local_inputs,args.output)

if __name__=="__main__":
    if not args.paper:
        density_figures()
    traits();other()
