from pathlib import Path
import argparse, hashlib
from functools import lru_cache
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from scipy.spatial.distance import pdist, squareform
from scipy.stats import t
p=argparse.ArgumentParser()
p.add_argument('books', type=Path)
p.add_argument('output', type=Path)
args=p.parse_args()
F=args.output; F.mkdir(parents=True,exist_ok=True)
FD=F/'source_data'; FD.mkdir(exist_ok=True)
@lru_cache(maxsize=None)
def df(k,n):
    return pd.read_excel(args.books/f'Supplementary_Data_S{k}.xlsx', sheet_name=n)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
C={'north':'#326b91','south':'#ba6248','spring':'#398866','autumn':'#997245','Halobacteriota':'#c25142','Other phyla':'#537b91'}
manifest=[]
def save(fig,name,sources,note=''):
    fig.canvas.draw()
    for ax in fig.axes:
        if not hasattr(ax,'_station_labels'):continue
        renderer=fig.canvas.get_renderer();used=[]
        legend=ax.get_legend()
        if legend is not None:used.append(legend.get_window_extent(renderer))
        positions=ax.transData.transform([p for p,_ in ax._station_labels])
        points=[Bbox.from_extents(x-4,y-4,x+4,y+4) for x,y in positions]
        for xy,label in ax._station_labels:
            placed=False
            for distance in [5,10,16,24,34,44]:
                for dx,dy in [(1,1),(-1,1),(1,-1),(-1,-1),(0,1),(0,-1),(1,0),(-1,0)]:
                    txt=ax.annotate(label,xy,xytext=(dx*distance,dy*distance),textcoords='offset points',fontsize=6.5,ha='left' if dx>0 else 'right' if dx<0 else 'center',va='bottom' if dy>=0 else 'top')
                    txt.set_in_layout(False);bb=txt.get_window_extent(renderer).expanded(1.1,1.15)
                    inside=ax.bbox.contains(bb.x0,bb.y0) and ax.bbox.contains(bb.x1,bb.y1)
                    if inside and not any(bb.overlaps(other) for other in used+points):
                        used.append(bb);placed=True
                        if distance>5:
                            lead=ax.annotate('',xy,xytext=(dx*distance,dy*distance),textcoords='offset points',arrowprops={'arrowstyle':'-','lw':.35,'color':'#999999'},zorder=1);lead.set_in_layout(False)
                        break
                    txt.remove()
                if placed:break
            assert placed,(name,label,'unplaced station label')
    fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight',pad_inches=.08)
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',pad_inches=.08)
    fig.savefig(F/(name+'.svg'),bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
    manifest.append({'figure':name,'action':'regenerated from source data','sources':sources,'note':note,'sha256':hashlib.sha256((F/(name+'.png')).read_bytes()).hexdigest()})
def panel(ax,label,title):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.14,1.03,label,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')
def export(data,name):data.to_csv(FD/(name+'.tsv'),sep='\t',index=False)
roles=df(1,'Sample_roles').set_index('sample_id');paired=roles.index[roles.paired_primary_n14].tolist()
def matrix(k,name,ids=paired,transpose=False):
    dat=df(k,name).set_index(df(k,name).columns[0])
    if transpose:dat=dat.T
    return dat.loc[ids].astype(float)
def pcoa(mat):
    a=np.asarray(mat,dtype=float)
    a=a[:,np.ptp(a,axis=0)>0]
    dis=squareform(pdist(a,'braycurtis'));n=len(a);J=np.eye(n)-np.ones((n,n))/n
    G=-.5*J@dis**2@J
    vals,vec=np.linalg.eigh(G);idx=np.argsort(vals)[::-1];vals=vals[idx];vec=vec[:,idx]
    coord=vec[:,:2]*np.sqrt(np.maximum(vals[:2],0))
    return coord,100*vals[:2]/vals[vals>0].sum(),G
def r2(G,x,reduced):
    H=x@np.linalg.pinv(x);HR=reduced@np.linalg.pinv(reduced)
    return np.trace((H-HR)@G)/np.trace(G)
def draw_pcoa(ax,mat,ids,title,paired_lines=True):
    xy,per,G=pcoa(mat)
    info=roles.loc[ids]
    if paired_lines:
        for site in info.site_corrected.unique():
            take=np.flatnonzero(info.site_corrected.values==site)
            if len(take)==2:
                # Direction is autumn to spring, independent of row order.
                a=next(i for i in take if info.season.iloc[i]=='autumn');b=next(i for i in take if info.season.iloc[i]=='spring')
                ax.annotate('',xy[b],xy[a],arrowprops={'arrowstyle':'->','lw':.6,'color':'#aeb5b9'},zorder=1)
    marker={'spring':'o','autumn':'s','summer':'^','winter':'D'}
    for i,sid in enumerate(ids):
        row=info.loc[sid]
        ax.scatter(*xy[i],s=29,marker=marker[row.season],c=C[row.basin],edgecolors='white',linewidths=.4,zorder=3)
    if paired_lines:ax._station_labels=[(xy[i],info.loc[sid].site_corrected) for i,sid in enumerate(ids)]
    ax.set(xlabel=f'PCoA1 ({per[0]:.1f}%)',ylabel=f'PCoA2 ({per[1]:.1f}%)')
    ax.set_title(title,loc='left',pad=8)
    ax.margins(.2)
    handles=[Line2D([],[],marker='o',ls='',color=C[b],label=b.title()) for b in ['north','south']]
    handles += [Line2D([],[],marker=marker[s],ls='',color='#555555',label=s.title()) for s in info.season.unique()]
    ax.legend(handles=handles,frameon=False,ncol=2,loc='best',handletextpad=.4,columnspacing=.8)
    out=pd.DataFrame(xy,columns=['PCoA1','PCoA2']);out.insert(0,'sample_id',ids)
    return out,G

# Source-defined station coordinates and analytical sample units replace an unsupported workflow.
fig=plt.figure(figsize=(7.1,6.7),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,.9])
sm=df(1,'Sample_manifest_35')
for j,lake in enumerate(['ZBY','DXC']):
    ax=fig.add_subplot(gs[0,j]);s=sm[(sm.lake==lake)&(sm.matrix=='water')&(sm.pooling_status=='Not pooled')].drop_duplicates('resolved_station')
    for _,r in s.iterrows():
        station=r.resolved_station;col=C['north' if station.startswith('N') else 'south'] if lake=='ZBY' else '#398866'
        ax.scatter(r.longitude_degE,r.latitude_degN,c=col,s=33,zorder=3)
        off={'N1':(-14,6),'N2':(4,-10),'N3':(4,5),'N4':(3,6),'S1':(-15,-9),'S2':(-14,6),'S3':(3,-10),'S4':(3,6),'S5':(4,-7)}.get(station,(4,5))
        label=station+('*' if station in ['N5','S4','S5'] else '')
        ax.annotate(label,(r.longitude_degE,r.latitude_degN),xytext=off,textcoords='offset points',fontsize=8)
    ax.ticklabel_format(useOffset=False);ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)');ax.grid(alpha=.2)
    ax.margins(.25);panel(ax,chr(97+j),f'{lake} station coordinates')
    if lake=='ZBY':ax.text(.02,.03,'* Autumn metagenome only',transform=ax.transAxes,fontsize=7)
ax=fig.add_subplot(gs[1,:]);ax.axis('off');panel(ax,'c','Library design and analysed response layers')
catalogue = df(4,'MAG_catalogue_414')
group_order = [('ZBY','water'),('ZBY','sediment'),('DXC','water'),('DXC','sediment')]
library_counts = sm.groupby(['lake','matrix']).size()
genome_counts = catalogue.source_assembly.value_counts()
quality = catalogue.quality_tier.value_counts()
water_units = roles.loc[roles.paired_primary_n14,'site_corrected'].nunique()
dom_matched = sm.recommended_DOM_EEM_set_15.eq('Yes').sum()
grid=[
[f'{len(sm)} physical libraries',f'{len(genome_counts)} coassemblies',f'{len(catalogue)}-MAG catalogue'],
['\n'.join(' / '.join(f'{library_counts[lk,mt]} {lk} {mt}' for lk,mt in group_order[i:i+2]) for i in (0,2)),
 'ZBY water / ZBY sediment\nDXC water / DXC sediment',
 ' / '.join(str(genome_counts[group]) for group in ('ZBY_water','ZBY_sed','DXC_water','DXC_sed'))+
 f' MAGs\n{quality["HQ"]} project HQ / {quality["MQ"]} MQ'],
['Within-ZBY observations','Distinct response definitions','Trait and external analyses'],
[f'{library_counts["ZBY","water"]} descriptive water libraries\n{len(paired)} libraries in {water_units} station pairs\n{dom_matched} matched optical-DOM libraries','Read species and community KOs\nMAG composition and weighted modules\nStrict carriers and carrier prevalence',f'MAG-median predicted proteome pI\nMarker and phylogenetic models\n{len(df(4,"QXLSG_lake_means"))} external lake means'],
]
tab=ax.table(cellText=grid,cellLoc='left',colWidths=[.33,.34,.33],bbox=[0,.08,1,.86])
tab.auto_set_font_size(False);tab.set_fontsize(8)
for (r,c),cell in tab.get_celld().items():
    cell.set_edgecolor('#b8c1c5');cell.set_linewidth(.5);cell.PAD=.04
    if r in [0,2]:cell.set_facecolor('#e5eeef');cell.set_text_props(weight='bold')
    else:cell.set_height(.22)
ax.text(0,0,'Counts denote libraries, genomes or lakes as specified; they are not interchangeable replicates.',fontsize=7,transform=ax.transAxes)
export(sm,'Figure_1_sample_manifest');save(fig,'Figure_1',['S1:Sample_manifest_35','S4:MAG_catalogue_414'],'Coordinate panels use the confirmed sample manifest; no inferred lake boundaries or flow paths.')

# Species paired ordination and station omissions.
sp=matrix(1,'Species_counts_35',transpose=True);sp=sp.loc[:,(sp>0).sum()>=2];sp=sp.div(sp.sum(axis=1),axis=0)
assert sp.shape[0]==len(paired) and sp.shape[1]>0
fig,axs=plt.subplots(1,2,figsize=(7.1,3.2),layout='constrained',gridspec_kw={'width_ratios':[1,1.05]})
coords,G=draw_pcoa(axs[0],sp,paired,'Paired species profiles')
season=(roles.loc[paired].season=='spring').astype(float).values;basin=(roles.loc[paired].basin=='south').astype(float).values
calc=r2(G,np.c_[np.ones(14),basin,season],np.c_[np.ones(14),basin]);assert np.isfinite(calc)
los=df(1,'Leave_one_station');ys=np.arange(len(los));a=axs[1]
a.scatter(los.R2,ys,c=['#23282b']+[C['north']]*4+[C['south']]*3,s=30)
for i,r in los.iterrows():a.text(.385,i,f"P = {r.P_exact:.6f}".rstrip('0'),va='center',fontsize=7)
a.axvline(los.R2.iloc[0],ls=':',lw=.7,c='#555555');a.set_yticks(ys,['Full set']+[f'Omit {x}' for x in los.omitted_site.iloc[1:]])
a.set(xlim=(.18,.55),xlabel='Seasonal term-specific R²',ylim=(7.7,-.7));panel(axs[0],'a','Paired species profiles');panel(a,'b','Station-omission sensitivity')
export(coords,'Figure_S1_PCoA');export(los,'Figure_S1_omissions');save(fig,'Figure_S1',['S1:Species_counts_35','S1:Sample_roles','S1:Leave_one_station'],'PCoA percentages use positive eigenvalues; R² independently reproduced.')

# Corrected paired KO ordination and complete technical models.
ko=matrix(3,'KO_abundance_35');ko=ko.loc[:,np.ptp(ko,axis=0)>0]
assert ko.shape[1]==41
fig,axs=plt.subplots(2,2,figsize=(7.1,6.0),layout='constrained')
at=df(3,'KO_taxonomic_attribution');a=axs[0,0]
bottom=np.zeros(2)
for cl,col in [('Named taxa','#ba6248'),('Unclassified','#b8c2c6')]:
    vals=[at[(at['class']==cl)&(at.basin==b)].mean_CPM.iloc[0] for b in ['north','south']]
    a.bar(['North','South'],vals,bottom=bottom,color=col,label=cl,width=.55);bottom+=vals
a.set_ylabel('Mean module abundance (CPM)');a.legend(frameon=False);panel(a,'a','Taxonomic attribution')
coords,G=draw_pcoa(axs[0,1],ko,paired,'Paired KO profiles')
calc=r2(G,np.c_[np.ones(14),basin,season],np.c_[np.ones(14),basin]);assert np.isfinite(calc)
panel(axs[0,1],'b','Paired KO profiles')
models=df(3,'PERMANOVA_models');a=axs[1,0]
for term,off,col in [('season',-.12,'#398866'),('basin',.12,'#ba6248')]:
    tt=models[models.term==term]
    for j,(_,r) in enumerate(tt.iterrows()):
        y=j+off;a.plot([r.CI95_low,r.CI95_high],[y,y],c=col,lw=1)
        a.scatter(r.R2,y,facecolor=col if r.P_exact<.05 else 'white',edgecolor=col,s=22,zorder=3)
        a.text(.70,y,f'{r.P_exact:.4f}',va='center',fontsize=6,color=col)
a.set_yticks(range(6),['M0: unadjusted','M1: catalogue detection','M2: mapping rate','M3: detected KOs','M4: mapped fraction','M5: M3 + M4']);a.set(xlim=(0,.91),ylim=(5.7,-.7),xlabel='Term-specific R² with 95% CI')
a.set_ylim(6.15,-.7)
a.legend(handles=[Line2D([],[],color='#398866',marker='o',label='Season'),Line2D([],[],color='#ba6248',marker='o',label='Basin')],frameon=False,loc='lower right',ncol=2);a.text(.70,-.53,'Exact P',fontsize=6)
panel(a,'c','Six model specifications')
qc=df(3,'KO_QC_35');a=axs[1,1]
for (lake,hab),rows in qc.groupby(['lake','matrix']):
    a.scatter(rows.mapped_fraction_pct,rows.detected_panel_KOs,label=f'{lake} {hab}',s=22,marker='o' if hab=='water' else '^',alpha=.8)
r=qc[qc.sample_id=='ZC1N2'].iloc[0];a.annotate('ZC1N2',(r.mapped_fraction_pct,r.detected_panel_KOs),xytext=(6,6),textcoords='offset points',fontsize=6)
a.set(xlabel='HUMAnN mapped fraction (%)',ylabel='Detected panel KOs');a.legend(frameon=False,fontsize=6);panel(a,'d','Coverage across 35 libraries')
export(coords,'Figure_S2_PCoA');export(models,'Figure_S2_models');save(fig,'Figure_S2',['S3:KO_abundance_35','S3:PERMANOVA_models','S3:KO_taxonomic_attribution','S3:KO_QC_35'],'Ordination uses unstratified target CPM, without an extra TSS; paired seasonal R² independently reproduced.')

mag=matrix(3,'MAG_abundance_n21',ids=roles.index.tolist(),transpose=True)
for num,ids,title in [(3,roles.index.tolist(),'All 21 ZBY water libraries'),(4,paired,'Seven paired ZBY stations')]:
    fig,a=plt.subplots(figsize=(5.8,4.0),layout='constrained');coords,_=draw_pcoa(a,mag.loc[ids],ids,title,num==4)
    export(coords,f'Figure_S{num}_PCoA');save(fig,f'Figure_S{num}',['S3:MAG_abundance_n21','S1:Sample_roles'],'Fixed 269-MAG primary catalogue; new ordination, not relabelled legacy 414-MAG coordinates.')

disp=df(3,'MAG_dispersion');fig,axs=plt.subplots(1,2,figsize=(7.1,3.6),layout='constrained')
names={'mag_relative_abundance_269':'MAG composition (269)','graded_weighted_modules_269':'Graded modules (269)','strict_carrier_abundance_269':'Strict carriers (269)','detected_carrier_fraction_269':'Carrier prevalence (269)','graded_weighted_modules_414':'Graded modules (414)','strict_carrier_abundance_414':'Strict carriers (414)','detected_carrier_fraction_414':'Carrier prevalence (414)'}
for a,term,letter in zip(axs,['season','basin'],'ab'):
    tt=disp[disp.term==term]
    for j,(_,r) in enumerate(tt.iterrows()):
        a.barh(j,r.dispersion_F,color='#398866' if term=='season' else '#ba6248',height=.5)
        a.text(max(r.dispersion_F+.08,.18),j,f'P = {r.dispersion_P_exact:.4f}',va='center',fontsize=7)
    a.set_yticks(range(7),[names[x] for x in tt.response]);a.invert_yaxis();a.set_xlabel('Dispersion F statistic')
    a.set_xlim(0,11 if term=='season' else .19);panel(a,letter,term.title()+' dispersion')
    if term=='basin':
        a.set_xlim(0,.24)
        for tx in a.texts:
            if tx.get_text().startswith('P ='):tx.set_x(.115)
save(fig,'Figure_S5',['S3:MAG_dispersion'],'Displays the exact paired diagnostic statistics, not obsolete four-season distances.')

# External figure with fixed group order and actual points for the two-MAG group.
lm=df(4,'QXLSG_lake_means');mods=df(4,'QXLSG_models');records=df(4,'Fig5d_MAG_records')
fig,axs=plt.subplots(2,2,figsize=(7.1,6.2),layout='constrained')
a=axs[0,0];x=lm.log10_salinity.to_numpy();y=lm.mean_MAG_median_pI.to_numpy();X=np.c_[np.ones(10),x];beta=np.linalg.lstsq(X,y,rcond=None)[0];s2=np.sum((y-X@beta)**2)/8
xp=np.linspace(x.min()-.1,x.max()+.1,150);Xp=np.c_[np.ones(150),xp];pr=Xp@beta;ci=t.ppf(.975,8)*np.sqrt(s2*np.einsum('ij,jk,ik->i',Xp,np.linalg.inv(X.T@X),Xp))
a.fill_between(xp,pr-ci,pr+ci,color='#326b91',alpha=.12);a.plot(xp,pr,c='#326b91',lw=1)
a.scatter(x,y,s=lm.n_attributable_MAGs*.16,c='#326b91',edgecolor='white',linewidth=.5,zorder=3)
for _,r in lm.iterrows():a.annotate(r.Abbreviation,(r.log10_salinity,r.mean_MAG_median_pI),xytext=(4,5),textcoords='offset points',fontsize=7)
a.set(xlabel=r'$\log_{10}$ salinity (%)',ylabel='Lake mean of MAG-median pI',ylim=(5.38,6.72));panel(a,'a','Ten equally weighted lake means')
a.text(.03,.04,f"OLS P = {mods.p_value.iloc[0]:.3f}",transform=a.transAxes,color='#326b91',fontsize=7)
for n in [100,300,600]:a.scatter([],[],s=n*.16,c='#326b91',label=str(n))
a.legend(title='MAGs',frameon=False,loc='lower left',bbox_to_anchor=(.01,.1),labelspacing=.7)
a=axs[0,1]
for i,r in mods.iterrows():a.errorbar(r.estimate,i,xerr=[[r.estimate-r.conf_low],[r.conf_high-r.estimate]],fmt='o',color='#326b91' if i==0 else '#ba6248',capsize=2)
a.axvline(0,c='#555555',ls='--',lw=.7);a.set_yticks([0,1],['Lake-mean OLS','MAG-level LMM']);a.set(xlim=(-.48,.04),ylim=(1.6,-.6),xlabel=r'pI slope per $\log_{10}$ salinity');panel(a,'b','Primary and clustering models')
for i,r in mods.iterrows():a.text(-.46,i+.25,f'P = {r.p_value:.3f}',fontsize=7)
a=axs[1,0];lo=df(4,'QXLSG_LOLO_OLS');ll=df(4,'QXLSG_LOLO_LMM');lakes=lm.Abbreviation.tolist()
omitted_col=next(c for c in lo.columns if 'omit' in c.lower() or 'exclude' in c.lower())
est_col=next(c for c in lo.columns if c in ['estimate','slope','OLS_slope','ols_estimate'])
for dat,key,est,off,col,label,base in [(lo,omitted_col,est_col,-.12,'#326b91','OLS',mods.estimate.iloc[0]),(ll,'excluded_lake','lmm_estimate',.12,'#ba6248','LMM',mods.estimate.iloc[1])]:
    yy=[lakes.index(v)+off for v in dat[key]];a.scatter(dat[est],yy,c=col,s=24,label=label);a.axvline(base,c=col,ls=':',lw=.8)
a.set_yticks(range(10),lakes);a.set(xlabel='Leave-one-lake-out slope',ylim=(10.4,-.7));a.legend(frameon=False,ncol=2,loc='lower right');panel(a,'c','Lake-omission estimates')
a=axs[1,1];rng=np.random.default_rng(20260910)
groups=[('<=4% salinity','Other phyla'),('<=4% salinity','Halobacteriota'),('>4% salinity','Other phyla'),('>4% salinity','Halobacteriota')]
counts=[]
for i,(sc,lin) in enumerate(groups):
    vals=records[(records.salinity_class==sc)&(records.lineage==lin)].median_pI.to_numpy();counts.append(len(vals));col=C[lin]
    if len(vals)>2:
        v=a.violinplot([vals],positions=[i],widths=.8,showextrema=False)
        for b in v['bodies']:b.set(facecolor=col,edgecolor=col,alpha=.2)
        a.boxplot([vals],positions=[i],widths=.17,showfliers=False,medianprops={'color':'#1d272c'},boxprops={'color':col},whiskerprops={'color':col},capprops={'color':col})
    jitter=rng.uniform(-.13,.13,len(vals)) if len(vals)>2 else np.array([-.04,.04])
    a.scatter(i+jitter,vals,c=col,s=3 if len(vals)>20 else 18,alpha=.16 if len(vals)>20 else .85,zorder=2,rasterized=True)
    a.text(i,9.65,f'n = {len(vals):,}',ha='center',fontsize=6.5)
assert sum(counts)==len(records)
a.axhline(df(4,'ZBY_proteome_pI').query("gtdb_phylum=='Halobacteriota'").median_pI.median(),c='#c25142',ls='--',lw=.7);a.set_xticks(range(4),['≤4%\nOther','≤4%\nHalo.','>4%\nOther','>4%\nHalo.']);a.set(xlim=(-.6,3.6),ylim=(4.1,9.9),ylabel='MAG-median pI');panel(a,'d','Salinity class and lineage')
save(fig,'Figure_5',['S4:QXLSG_lake_means','S4:QXLSG_models','S4:QXLSG_LOLO_OLS','S4:QXLSG_LOLO_LMM','S4:Fig5d_MAG_records'],'Group counts asserted; OLS curve refitted from ten means; point sizes retain MAG count.')

# DOM sensitivity keeps both denominators and their implemented permutations visible.
dom=df(4,'DOM_conditional_models');fig,axs=plt.subplots(1,2,figsize=(6.5,3.0),layout='constrained')
for a,metric,letter in zip(axs,['adjusted_R2','P'],'ab'):
    for i,r in dom.iterrows():
        a.barh(i,r[metric],height=.5,color='#398866' if r.n==15 else '#bc8062')
        a.text(r[metric]+(.003 if metric=='adjusted_R2' else .015),i,f'{r[metric]:.4f}',va='center',fontsize=7)
    a.set_yticks(range(4),['n = 15 / Bray–Curtis','n = 15 / Hellinger','n = 19 / Bray–Curtis','n = 19 / Hellinger']);a.invert_yaxis()
    a.set_xlabel('Adjusted R²' if metric=='adjusted_R2' else 'Permutation P')
    a.set_xlim(0,.105 if metric=='adjusted_R2' else .62);panel(a,letter,'Joint DOM terms' if metric=='adjusted_R2' else 'Design-specific tests')
    if metric=='P':a.axvline(.05,c='#555555',ls=':',lw=.8)
save(fig,'Figure_S14',['S4:DOM_conditional_models','S4:DOM_archived_ANOVA'],'n15 uses 127 supplied site-block configurations; n19 uses 9999 unrestricted permutations.')

# Cross-lake source attribution and the actual adjusted archived results; no old within-ZBY test.
cross=df(4,'Crosslake_model_outputs');fig,axs=plt.subplots(1,3,figsize=(7.1,3.1),layout='constrained')
att=df(4,'Crosslake_source_attribution');att=att[att.samp_matrix=='water'];a=axs[0]
for i,lk in enumerate(['ZBY','DXC']):
    vals=att[att.samp_lake==lk].pct_own_lake.to_numpy();a.scatter(i+np.random.default_rng(i).uniform(-.12,.12,len(vals)),vals,s=15,c='#326b91' if i==0 else '#ba6248');a.plot([i-.18,i+.18],[np.median(vals)]*2,c='#333333')
a.set_xticks([0,1],['ZBY','DXC']);a.set(xlim=(-.5,1.5),ylabel='Same-lake source abundance (%)');panel(a,'a','Catalogue attribution')
rr=['MAG composition','graded_weighted','strict_carrier_abundance','detected_carrier_fraction'];labels=['MAG composition','Graded modules','Strict carriers','Carrier prevalence']
for i,rsp in enumerate(rr):
    ss=cross[(cross.response==rsp)&(cross.term=='lake')&cross.model.str.contains('n_detected')].iloc[0]
    axs[1].scatter(ss.R2,i,s=27,c='#326b91');axs[1].text(.255,i,f'P = {ss.P:.3f}',va='center',fontsize=6)
    dd=cross[(cross.response==rsp)&(cross.term=='Groups')].iloc[0]
    axs[2].scatter(dd.F,i,s=27,c='#ba6248');axs[2].text(.14,.90-i*.235,f'{labels[i]}\nP = {dd.P:.3g}',transform=axs[2].transAxes,fontsize=6,va='top')
axs[1].set_yticks(range(4),labels);axs[1].set(xlim=(0,.41),ylim=(3.6,-.6),xlabel='Adjusted-model lake R²');panel(axs[1],'b','Sampled-group separation')
axs[2].set_yticks(range(4),['']*4);axs[2].set(xlim=(0,450),ylim=(3.6,-.6),xlabel='Dispersion F statistic');panel(axs[2],'c','Unequal cross-lake dispersion')
save(fig,'Figure_S16',['S4:Crosslake_source_attribution','S4:Crosslake_model_outputs'],'Cross-lake P values are archived descriptive free-permutation fits, not replicated lake inference.')

# Preserve count annotations with readable contrast across the colour scale.
neighbour=df(4,'Neighbourhood_rules')
rules=['diff1','diff2','same_contig','mag_cooccur']
sets=['all_hits','kofam187','pfam189','intersect182','pfam_conf164']
mag=neighbour.pivot(index='rule',columns='cset',values='n_MAG').loc[rules,sets]
contig=neighbour.pivot(index='rule',columns='cset',values='n_contig').loc[rules,sets]
assert mag.shape==(len(rules),len(sets))
fig,a=plt.subplots(figsize=(7.1,4.45),layout='constrained')
im=a.imshow(mag.to_numpy(dtype=float),cmap='plasma_r',vmin=100,vmax=175,aspect='auto')
for i in range(4):
    for j in range(5):
        rgb=np.asarray(im.cmap(im.norm(float(mag.iloc[i,j])))[:3])
        linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
        lum=float(linear@np.array([.2126,.7152,.0722]))
        foreground='#151515' if lum>.18 else 'white'
        c='NA' if pd.isna(contig.iloc[i,j]) else str(int(contig.iloc[i,j]))
        a.text(j,i,f'{int(mag.iloc[i,j])} MAGs\n{c} contigs',ha='center',va='center',fontsize=8,color=foreground,linespacing=1.5)
a.set_xticks(range(5),sets,rotation=25,ha='right')
a.set_yticks(range(4),rules)
a.set_xticks(np.arange(-.5,5,1),minor=True);a.set_yticks(np.arange(-.5,4,1),minor=True)
a.grid(which='minor',color='white',linewidth=2);a.tick_params(which='minor',bottom=False,left=False)
a.set(xlabel='Carrier evidence set',ylabel='Neighbour rule')
a.set_title('Sensitivity of potassium-marker neighbourhood counts',loc='left',pad=10)
fig.colorbar(im,ax=a,shrink=.72,pad=.025,label='MAG count')
export(neighbour,'Figure_S11_neighbourhood_counts')
save(fig,'Figure_S11',['S4:Neighbourhood_rules'],'All 20 MAG and contig count cells retained; text colour uses cell luminance for contrast.')
