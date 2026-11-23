"""Draw the submitted ion and carrier panels from the local numerical sources."""
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


def render(inputs, output):
    def read(name, index=None):
        return pd.read_csv(inputs/"extra"/(name+".tsv"),sep="\t",index_col=index)

    def export(fig, stem):
        for suffix in ("pdf","svg","png"):
            fig.savefig(output/(stem+"."+suffix),dpi=300,bbox_inches="tight")
        plt.close(fig)

    P=json.loads((inputs/"si_display/palette.json").read_text())
    N,S=P["col_north"],P["col_south"]
    TEAL,PALE,MID,INK=P["col_teal"],P["col_pale"],P["col_mid"],P["col_ink"]
    def clean(ax, grid=True):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2.6, width=.5)
        if grid:
            ax.set_axisbelow(True)
            ax.grid(axis="x", color=P["col_light"], lw=.4)

    def heading(fig, letter, title, x, y, subtitle=None):
        fig.text(x, y, letter, fontsize=12, weight="bold", va="top")
        fig.text(x + .035, y-.002, title, fontsize=8.5, weight="bold", va="top")
        if subtitle:
            fig.text(x + .035, y-.024, subtitle, fontsize=6.7, color=MID, va="top")

    RESP = [
        ("MAG_relative_abundance_269", "MAG composition · 269"),
        ("graded_modules_269", "GMP · 269"),
        ("strict_modules_269", "SCA · 269"),
        ("carrier_fraction_269", "Prevalence · 269"),
        ("graded_modules_414", "GMP · 414"),
        ("strict_modules_414", "SCA · 414"),
        ("carrier_fraction_414", "Prevalence · 414"),
    ]
    ec=read("environment_n14")
    ld=read("ion_PCA_loadings")
    ions={"K":r"K$^+$","Br":r"Br$^-$","Cl":r"Cl$^-$","Na":r"Na$^+$","F":r"F$^-$",
          "SO4":r"SO$_4^{2-}$","Li":r"Li$^+$","Ca":r"Ca$^{2+}$","Mg":r"Mg$^{2+}$"}
    fig=plt.figure(figsize=(7.48,7.60))
    heading(fig,"a","Salinity and ionic composition",.02,.977,
            "14 paired libraries; lines join stations")
    ax=fig.add_axes([.105,.586,.335,.315])
    for site,g in ec.groupby("site"):
        ax.plot(g.salinity_z,g.ion_PC1_z,color="#B4BEC5",lw=.8,zorder=1)
    for _,r in ec.iterrows():
        ax.scatter(r.salinity_z,r.ion_PC1_z,c=N if r.basin=="north" else S,
                   marker="o" if r.season=="spring" else "s",s=24,ec="white",lw=.4,zorder=3)
    ax.set_xlabel("Salinity (standardized)");ax.set_ylabel("Ion PC1 (standardized)")
    clean(ax)
    ax.legend(handles=[Line2D([],[],color=N,marker="o",ls="",label="NZ"),
                       Line2D([],[],color=S,marker="o",ls="",label="SZ"),
                       Line2D([],[],color=MID,marker="o",ls="",label="Spring"),
                       Line2D([],[],color=MID,marker="s",ls="",label="Autumn")],
              loc="lower left",bbox_to_anchor=(-.18,-.29),ncol=4,fontsize=6.3,
              handletextpad=.25,columnspacing=.6)
    pc1 = read("ion_PCA_variance").set_index("axis").at["PC1","variance_fraction"]
    heading(fig,"b","Ionic-composition loadings",.51,.977,
            f"PC1 explains {pc1:.1%} of CLR ion variance")
    ax=fig.add_axes([.60,.586,.37,.315])
    ld=ld.set_index("ion").loc[["K","Br","Cl","Na","F","SO4","Li","Ca","Mg"]]
    ax.barh(range(9),ld.PC1,color=[S if x>0 else N for x in ld.PC1],height=.65)
    ax.axvline(0,color=MID,lw=.6);ax.set_yticks(range(9),[ions[x] for x in ld.index])
    ax.set_xticks([-.6,-.4,-.2,0,.2,.4]);ax.set_xlim(-.71,.42)
    ax.invert_yaxis();ax.set_xlabel("Loading on CLR-derived Ion PC1");clean(ax)
    heading(fig,"c","Conditional ion-axis associations",.02,.462,
            "Exploratory; BH correction across seven tests")
    ax=fig.add_axes([.245,.095,.255,.282])
    ionrows=read("FigureS6c_ion_effects")
    assert len(ionrows)==7
    for i,(key,label) in enumerate(RESP):
        r=ionrows[ionrows.response==key].iloc[0]
        color=TEAL if key.endswith("269") else P["col_autumn"]
        ax.plot([0,r.R2],[i,i],color="#DFE3E7",lw=1.2)
        ax.scatter(r.R2,i,marker="o" if key.endswith("269") else "s",s=25,c=color,zorder=3)
        ax.text(.278,i,f"{r.P_BH_within_model_term:.3f}",ha="right",va="center",fontsize=6.7)
    ax.set_yticks(range(7),[label for _,label in RESP],fontsize=6.7)
    ax.set_ylim(6.6,-.6);ax.set_xlim(0,.285);ax.set_xlabel("Ion-PC1 term-specific $R^2$")
    ax.text(.278,-.96,"BH q",ha="right",fontsize=6.7,color=MID)
    ax.tick_params(axis="y",length=0);clean(ax)
    heading(fig,"d","Basin term across conditioning models",.54,.462,
            "Cells show R² / nominal restricted P")
    ax=fig.add_axes([.655,.095,.32,.282])
    mods=["baseline","salinity_conditioned","salinity_ion_PC1_conditioned"]
    ordered_keys = [key for key, label in RESP]
    vals = read("FigureS6d_basin_R2", 0).loc[ordered_keys,mods].to_numpy()
    pvals = read("FigureS6d_basin_P", 0).loc[ordered_keys,mods].to_numpy()
    im=ax.pcolormesh(np.arange(4)-.5,np.arange(8)-.5,vals,cmap="Greys",vmin=0,vmax=.5,rasterized=False)
    for i in range(7):
        for j in range(3):
            ax.text(j,i,f"{vals[i,j]:.3f}\n{pvals[i,j]:.3f}",ha="center",va="center",fontsize=6.4)
    ax.set_yticks(range(7),[label for _,label in RESP],fontsize=6.7)
    ax.set_xticks(range(3),["Baseline","+ salinity","+ salinity\n+ Ion PC1"],fontsize=6.6)
    ax.set_ylim(6.5,-.5);ax.tick_params(length=0,pad=5)
    ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,7,1),minor=True)
    ax.grid(which="minor",color="white",lw=1);ax.tick_params(which="minor",length=0)
    vif = read("predictor_VIF_n14").query("term == 'ion_PC1_z'").set_index("model").VIF
    fig.text(.07,.020,f"Ion-PC1 VIF = {vif['ion_PC1_conditioned']:.2f} in the ion model and {vif['salinity_ion_PC1_conditioned']:.2f} with salinity.",
             fontsize=6.7,color=MID)
    export(fig,"Figure_S6")

    ct=read("display_top16_carriers").set_index("feature")
    magids=ct.index.tolist()
    meta=read("paired_manifest_14").set_index("sample_id")
    ab=read("MAG_abundance_sample14",0).loc[meta.index,magids]*100
    fig=plt.figure(figsize=(7.48,6.20))
    heading(fig,"","MAG-level carriers of the displayed contrasts",.015,.983,
            "Area = relative abundance; circles = spring, squares = autumn")
    ax=fig.add_axes([.348,.18,.444,.67])
    aq=fig.add_axes([.842,.18,.063,.67],sharey=ax)
    ao=fig.add_axes([.936,.18,.045,.67],sharey=ax)
    scale=115/float(ab.max().max())
    maglabels=[]
    for i,mag in enumerate(magids):
        for z in [ax,aq,ao]:
            if i%2==0:z.axhspan(i-.5,i+.5,color=PALE)
        r=ct.loc[mag]
        genus=r.gtdb_genus
        if genus=="Unclassified":
            genus="Uncl. "+str(r.gtdb_family if r.gtdb_family!="Unclassified" else r.gtdb_order)
        maglabels.append(f"{genus} | {mag.split('__')[1].replace('metabat2_','')}")
        for j,sample in enumerate(meta.index):
            value=ab.loc[sample,mag]
            if value<=0:continue
            spring=meta.loc[sample,"season"]=="spring"
            ax.scatter(j,i,s=value*scale*(4/np.pi if spring else 1),
                       marker="o" if spring else "s",
                       c=N if meta.loc[sample,"basin"]=="north" else S,
                       ec="white",lw=.3,alpha=.92,zorder=3)
        ax.add_patch(Rectangle((-.97,i-.5),.17,1,
                              fc=P["col_autumn"] if r.domain_gtdb=="Archaea" else TEAL,
                              ec="white",lw=.2,clip_on=False))
        aq.barh(i,r.completeness,height=.55,color=TEAL)
        aq.text(50,i,f"{r.completeness:.0f}",ha="center",va="center",fontsize=6.1,color="white")
        ao.text(.5,i,f"{int(r.detected_stations)}/7",ha="center",va="center",fontsize=6.6)
    ax.set_yticks(range(16),maglabels,fontsize=6.5)
    ax.tick_params(axis="y",length=0,pad=10)
    ax.set_xticks(range(14),meta.index,rotation=55,ha="right",fontsize=6.6)
    ax.set_xlim(-.55,13.55)
    ax.set_xticks(np.arange(-.5,14,1),minor=True);ax.set_yticks(np.arange(-.5,16,1),minor=True)
    ax.grid(which="minor",color="#E0E5E8",lw=.4);ax.tick_params(which="minor",length=0)
    ax.axvline(7.5,color="white",lw=2)
    for j in range(7):
        ax.add_patch(Rectangle((2*j-.5,-1.18),2,.5,fc=N if j<4 else S,ec="white",lw=.5,clip_on=False))
        ax.text(2*j+.5,-.93,meta.station.iloc[2*j],color="white",fontsize=6.6,ha="center",va="center")
    for z in [ax,aq,ao]:z.set_ylim(15.6,-.6)
    aq.set_xlim(0,100);aq.set_xticks([0,100],["0","100"],fontsize=6)
    aq.tick_params(axis="y",left=False,labelleft=False)
    aq.spines[["top","right","left"]].set_visible(False)
    aq.text(.5,1.035,"Completeness\n(%)",transform=aq.transAxes,ha="center",fontsize=6.6)
    ao.axis("off");ao.set_xlim(0,1)
    ao.text(.5,1.035,"Detected\nstations",transform=ao.transAxes,ha="center",fontsize=6.6)
    fig.legend(handles=[Line2D([],[],color=P["col_autumn"],lw=5,label="Archaea"),
                        Line2D([],[],color=TEAL,lw=5,label="Bacteria")],
               loc="lower left",bbox_to_anchor=(.02,.045),fontsize=6.7)
    sizehandles=[plt.scatter([],[],s=x*scale*4/np.pi,fc="#89949C",ec="white",lw=.3,label=f"{x}%") for x in [1,10,30]]
    fig.legend(handles=sizehandles,title="Relative abundance",loc="lower left",bbox_to_anchor=(.34,.025),
               ncol=3,fontsize=6.7,title_fontsize=6.7,handletextpad=.7,columnspacing=1.2)
    tests=read("MAG_differential")
    fig.text(.70,.059,f"All {len(tests)} MAG tests:\nBH q ≥ {tests.q_BH.min():.3f}",fontsize=6.7,color=MID)
    export(fig,"Figure_S19")
