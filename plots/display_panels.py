"""Redraw S8/S17/S18 from frozen numerical sources; no model is refitted.

Historical module IDs are retained. Only the ureABC display label changes.
The source grid, ranked contrasts, intervals, P values and q values are checked.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, textwrap, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

parser=argparse.ArgumentParser()
parser.add_argument('inputs',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
SRC=args.inputs
FIG=args.output
FIG.mkdir(parents=True,exist_ok=True)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
denpath=SRC/'Figure_S8_density.tsv';rankpath=SRC/'Figure_S17_ranked_modules.tsv';palette=SRC/'palette.json'

P = json.loads(palette.read_text())
N, S, INK, MID, PALE = [P[k] for k in ['col_north', 'col_south', 'col_ink', 'col_mid', 'col_pale']]
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'axes.labelsize':8,
    'axes.titlesize':9, 'axes.titleweight':'bold', 'axes.spines.top':False,
    'axes.spines.right':False, 'axes.linewidth':.6, 'xtick.labelsize':7,
    'ytick.labelsize':7, 'legend.fontsize':7, 'legend.frameon':False,
    'pdf.fonttype':42, 'svg.fonttype':'none', 'savefig.facecolor':'white'})
exports = []

def export(fig, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Matplotlib creates tick Text objects beyond the view limits but does not draw them.
    undrawn = set()
    for ax in fig.axes:
        for axis in [ax.xaxis, ax.yaxis]:
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not ax.axison or not axis.get_visible() or not lo-1e-9 <= tick.get_loc() <= hi+1e-9:
                    undrawn.update([tick.label1,tick.label2])
            if not ax.axison or not axis.get_visible():
                undrawn.update([axis.label,axis.offsetText])
    visible = [t for t in fig.findobj(matplotlib.text.Text) if t not in undrawn and t.get_visible() and t.get_text().strip()]
    outside = []
    for t in visible:
        box = t.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1:
            outside.append(t.get_text())
    assert not outside, (name, outside)
    for suffix in ['png','pdf','svg']:
        fig.savefig(FIG / f'{name}_revised.{suffix}', dpi=600)
    fig.savefig(FIG / f'{name}_revised_preview.png', dpi=160)
    exports.append(dict(figure=name, sha256=digest(FIG/f'{name}_revised.png'),
                        outside_canvas=outside, size_inches=list(fig.get_size_inches())))
    plt.close(fig)

den = pd.read_csv(denpath, sep='\t')
fig, ax = plt.subplots(figsize=(5.6,3.55), layout='constrained')
for group, dd in den.groupby('group'):
    ax.plot(dd.x, dd.y, label=group, color={'Halobacteriota':'#c25142','Other phyla':'#537b91'}.get(group,'#555555'), lw=1.6)
ax.set(xlabel='Predicted protein pI', ylabel='Pooled-protein density', xlim=(2,11))
ax.set_title('Pooled proteins\nwithout equal-MAG weighting', loc='left', pad=9)
ax.legend()
export(fig, 'Figure_S8')

ct = pd.read_csv(SRC/'Crosslake_module_tests.tsv', sep='\t')
cg = pd.read_csv(SRC/'Crosslake_graded.tsv', sep='\t')
top = ct.loc[ct.delta.abs().nlargest(18).index].sort_values('delta')
pd.testing.assert_frame_equal(top.reset_index(drop=True), pd.read_csv(rankpath,sep='\t'), check_dtype=False)
fig = plt.figure(figsize=(7.48,8.1), layout='constrained')
gs = fig.add_gridspec(3,4,height_ratios=[3.3,1,1])
ax = fig.add_subplot(gs[0,:])
colors = dict(zip(ct.category.unique(),['#498197','#8370a4','#ba7650','#459072','#b45a63','#819b50']))
labels = top.label.str.replace('urea utilization','urea hydrolysis',case=False,regex=False)
assert (labels != top.label).sum() == 1
for j, (_, row) in enumerate(top.iterrows()):
    ax.plot([0,row.delta],[j,j], color=colors[row.category], lw=1)
    ax.scatter(row.delta,j, color=colors[row.category], s=22)
ax.set_yticks(range(len(top)),[textwrap.fill(s,37) for s in labels])
ax.invert_yaxis(); ax.axvline(0,color='#888888',lw=.7)
ax.set(xlim=(-1.08,1.08), xlabel="Cliff's delta (positive: higher in DXC)")
ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=k) for k,c in colors.items()],loc='upper right',fontsize=7,handletextpad=.4,labelspacing=.35)
ax.set_title('a  Exploratory module ranking', loc='left', pad=9)
salt = ['Salt_Kdp_complex','Salt_Mrp_complex','Salt_Trk_Ktr_regulatory','Salt_Trk_Ktr_core','Salt_ProU_transporter','Salt_ectoine_synthesis','Salt_betaine_transport','Salt_trehalose_synthesis']
titles = ['Kdp complex','Mrp complex','Trk/Ktr regulatory','Trk/Ktr core','ProU transporter','Ectoine synthesis','Betaine transport','Trehalose synthesis']
for j,(module,title) in enumerate(zip(salt,titles)):
    ax=fig.add_subplot(gs[1+j//4,j%4]);ax.set_title(('b  ' if j==0 else '')+title,fontsize=8)
    for i,lake in enumerate(['DXC','ZBY']):
        vals=cg.loc[cg.sample_id.str.startswith('D' if lake=='DXC' else 'Z'),module].to_numpy()
        assert len(vals)==(5 if lake=='DXC' else 21)
        ax.boxplot([vals],positions=[i],widths=.45,showfliers=False,medianprops={'color':'#333333'},boxprops={'color':'#777777'},whiskerprops={'color':'#777777'},capprops={'color':'#777777'})
        ax.scatter(i+np.random.default_rng(j*2+i).uniform(-.11,.11,len(vals)),vals,s=9,c='#ba6248' if i==0 else '#326b91',alpha=.7)
    ax.set_xticks([0,1],['DXC','ZBY']); ax.set_ylim(-.05,1.05)
    if j%4==0:ax.set_ylabel('Weighted score')
export(fig,'Figure_S17')

def read(name,index=None):
    return pd.read_csv(SRC/(name+'.tsv'),sep='\t',index_col=index)

LABELS={'Salt_Trk_Ktr_core':'Trk/Ktr marker (K03498)',
 'TCA_step8_malate_dh':'Malate dehydrogenase', 'C1_formate_ox':'Formate dehydrogenase',
 'N_urea_utilization':'Urea hydrolysis (ureABC)', 'Salt_betaine_synthesis':'Betaine synthesis (betAB)',
 'N_nitrate_red_periplasmic':'NapAB nitrate reduction',
 'TCA_step5_succinylCoA_synthetase_panelSucCD':'Succinyl-CoA synthetase',
 'N_nitrate_red_membrane':'NarGHI-like module', 'TCA_step2_aconitase':'Aconitase',
 'Salt_hydroxyectoine_ext':'Hydroxyectoine (ectD)',
 'TCA_step1_citrate_synthase':'Citrate synthase', 'Salt_betaine_transport':'Choline/betaine transport'}
top=read('display_top12_modules').set_index('feature'); modules=top.index.tolist()
sm=read('GMP_station7',0)
diag=read('module_direction_diagnostics');diag=diag[diag.response=='GMP'].set_index('module_id')
sens=read('module_catalogue_breadth_sensitivity');qual=read('module_MAG_quality_sensitivity')
def heading(fig,letter,title,x,y,subtitle):
    fig.text(x,y,letter,fontsize=12,weight='bold',va='top')
    fig.text(x+.035,y-.002,title,fontsize=8.5,weight='bold',va='top')
    fig.text(x+.035,y-.024,subtitle,fontsize=6.7,color=MID,va='top')
def clean(ax):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2.6,width=.5)
    ax.set_axisbelow(True);ax.grid(axis='x',color=P['col_light'],lw=.4)

fig=plt.figure(figsize=(7.48,9.05))
heading(fig,'a','Station-level distributions and module contrasts',.015,.982,
        'NZ: four stations; SZ: three stations. Each point averages spring and autumn.')
ax=fig.add_axes([.295,.566,.258,.34]);ae=fig.add_axes([.604,.566,.245,.34],sharey=ax)
at=fig.add_axes([.875,.566,.109,.34],sharey=ax)
for i,m in enumerate(modules):
    for z in [ax,ae,at]:
        if i%2==0:z.axhspan(i-.5,i+.5,color=PALE)
    for stations,off,color,marker in [(['N1','N2','N3','N4'],-.14,N,'o'),(['S1','S2','S3'],.14,S,'s')]:
        ax.scatter(sm.loc[stations,m],i+off+np.linspace(-.04,.04,len(stations)),s=14,c=color,marker=marker,ec='white',lw=.3,zorder=3)
    row=top.loc[m];color=N if row.delta_NZ_minus_SZ>0 else S
    ae.plot([row.CI95_low,row.CI95_high],[i,i],color=color,lw=1.2)
    ae.scatter(row.delta_NZ_minus_SZ,i,s=22,c=color,zorder=4)
    at.text(0,i,f'{row.P_exact:.3f}',va='center',fontsize=6.7)
    at.text(.56,i,f'{row.q_BH:.3f}',va='center',fontsize=6.7)
ax.set_yticks(range(12),[LABELS[m] for m in modules],fontsize=6.8)
ax.tick_params(axis='y',length=0,pad=5);ax.set_xlim(-.02,1.02)
ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('Graded potential')
ax.legend(handles=[Line2D([],[],color=N,marker='o',ls='',label='NZ'),Line2D([],[],color=S,marker='s',ls='',label='SZ')],loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,fontsize=6.8)
ae.set_xlim(-.65,.65);ae.set_xticks([-.6,-.3,0,.3,.6]);ae.axvline(0,color=MID,ls='--',lw=.6)
ae.tick_params(axis='y',labelleft=False,left=False);ae.set_xlabel('NZ − SZ (95% interval)')
at.axis('off');at.set_xlim(0,1);at.text(0,-1,'P',weight='bold',fontsize=7);at.text(.56,-1,'BH q',weight='bold',fontsize=7)
for z in [ax,ae,at]:z.set_ylim(11.6,-.6)
clean(ax);clean(ae)
fig.text(.295,.517,'Pointwise bootstrap intervals are not multiplicity-adjusted or inversions of the permutation test.',fontsize=6.5,color=MID)
heading(fig,'b','Direction and denominator sensitivity',.015,.465,
        'All cells show NZ − SZ; these are sensitivity contrasts, not additional independent tests.')
columns=['Primary','Spring','Autumn','LOSO min','LOSO max','414 MAGs','Breadth .10','Breadth .25','HQ fixed','HQ renorm.']
vv=[]
for m in modules:
    d=diag.loc[m]
    def sval(cat,breadth):
        z=sens[(sens.feature==m)&(sens.catalogue==cat)&np.isclose(sens.breadth,breadth)&(sens.family=='GMP')]
        assert len(z)==1
        return z.iloc[0].delta_NZ_minus_SZ
    def qval(den):
        z=qual[(qual.module_id==m)&(qual.quality_subset=='HQ_only')&(qual.denominator==den)&(qual.response=='GMP')]
        assert len(z)==1
        return z.iloc[0].delta_NZ_minus_SZ
    vv.append([top.loc[m].delta_NZ_minus_SZ,d.spring_delta,d.autumn_delta,d.LOSO_min_delta,d.LOSO_max_delta,
               sval(414,.5),sval(269,.1),sval(269,.25),qval('fixed269_denominator'),qval('renormalized_retained_MAGs')])
vv=np.array(vv)
np.testing.assert_allclose(vv,read('FigureS18b_direction_sensitivity',0).loc[modules,columns],rtol=1e-12,atol=1e-14)
vv*=100
ax=fig.add_axes([.295,.105,.685,.30])
norm=TwoSlopeNorm(vmin=-max(abs(vv).max(),1),vcenter=0,vmax=max(abs(vv).max(),1))
cmap=LinearSegmentedColormap.from_list('SZ_white_NZ',[S,'white',N])
ax.pcolormesh(np.arange(11)-.5,np.arange(13)-.5,vv,cmap=cmap,norm=norm,rasterized=False)
ax.set_ylim(11.5,-.5);ax.set_xticks(range(10),columns,rotation=40,ha='right',fontsize=6.5)
ax.set_yticks(range(12),[LABELS[m] for m in modules],fontsize=6.8);ax.tick_params(length=0,pad=5)
for i in range(12):
    for j in range(10):
        ax.text(j,i,f'{vv[i,j]:+.1f}',ha='center',va='center',fontsize=6.2,color='white' if abs(vv[i,j])>.68*norm.vmax else INK)
ax.set_xticks(np.arange(-.5,10,1),minor=True);ax.set_yticks(np.arange(-.5,12,1),minor=True)
ax.grid(which='minor',color='white',lw=.7);ax.tick_params(which='minor',length=0)
fig.text(.295,.020,'Blue: NZ higher; orange: SZ higher. Values are percentage-point differences in graded potential.',fontsize=6.6,color=MID)
export(fig,'Figure_S18')
