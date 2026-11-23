"""Fixed-reference MAG and genomic-potential matrices for assembly sensitivity."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path('@REVISION@')
group=sys.argv[1];cat=ROOT/'catalogues'/group;out=cat/'responses';out.mkdir(exist_ok=True)
meta=pd.read_csv(ROOT/'inputs/paired14.tsv',sep='\t')
score=pd.read_csv(cat/'annotation/MAG_module_scores.tsv',sep='\t',index_col=0)
complete=pd.read_csv(cat/'annotation/MAG_module_complete.tsv',sep='\t',index_col=0)
assert set(score.index)==set(complete.index) and score.shape[1]==49
mags=score.index;complete=complete.loc[mags,score.columns]
depth=pd.DataFrame(index=meta['sample'],columns=mags,dtype=float)
breadth=depth.copy();abundance=depth.copy();mapping={}
for sample in meta['sample']:
    p=cat/'coverm'/f'{sample}.tsv';assert p.with_suffix('.complete').exists(),sample
    df=pd.read_csv(p,sep='\t',index_col=0)
    assert set(mags).issubset(df.index),sample
    for dest,keyword in [(depth,' Mean'),(breadth,' Covered Fraction'),(abundance,' Relative Abundance (%)')]:
        cols=[c for c in df if c.endswith(keyword)];assert len(cols)==1,cols
        dest.loc[sample]=df.loc[mags,cols[0]].astype(float)
    mapping[sample]=100-float(df.loc['unmapped',[c for c in df if c.endswith(' Relative Abundance (%)')][0]])
assert np.isfinite(depth.values).all() and np.isfinite(breadth.values).all() and np.isfinite(abundance.values).all()
qc=[]
for cutoff in [.10,.25,.50]:
    detected=(depth>=1)&(breadth>=cutoff)
    a=abundance.where(detected,0)
    denominator=a.sum(axis=1)
    assert (denominator>0).all(),f'No detected abundance at cutoff {cutoff}; do not silently treat as ecological zero'
    a=a.div(denominator,axis=0)
    n=detected.sum(axis=1)
    matrices={'MAG':a,'GMP':pd.DataFrame(a.values@score.values,index=a.index,columns=score.columns),
      'SCA':pd.DataFrame(a.values@complete.values,index=a.index,columns=score.columns),
      'DCP':pd.DataFrame(detected.values.astype(float)@complete.values,index=a.index,columns=score.columns).div(n,axis=0)}
    for response,mat in matrices.items():
        assert np.isfinite(mat.values).all()
        mat.to_csv(out/f'{response}_breadth{int(100*cutoff)}.tsv',sep='\t',index_label='sample_id')
    for sample in a.index:
        qc.append(dict(sample_id=sample,catalogue=group,breadth_cutoff=cutoff,catalogue_MAGs=len(mags),
         n_detected_mags=int(n[sample]),detected_MAG_fraction=float(n[sample]/len(mags)),
         thresholded_abundance_sum_pct=float(denominator[sample]),
         prethreshold_catalogue_capture_pct=float(abundance.loc[sample].sum()),
         mapping_rate_pct=mapping[sample]))
pd.DataFrame(qc).to_csv(out/'sample_QC.tsv',sep='\t',index=False)
for name,df in [('mean_depth',depth),('covered_fraction',breadth),('unthresholded_abundance_pct',abundance)]:
    df.to_csv(out/f'{name}.tsv',sep='\t',index_label='sample_id')
(out/'definitions.json').write_text(json.dumps({'genome_depth_min':1,'breadth_cutoffs':[.10,.25,.50],
 'denominator':'Each fixed catalogue shared by all 14 samples; retained CoverM abundance sums to 1',
 'GMP':'abundance-weighted graded module score',
 'SCA':'abundance-weighted strict complete-carrier indicator',
 'DCP':'detected strict complete carriers / detected MAGs',
 'zero_detected_samples':'Raises error rather than interpreting non-recovery as ecological absence'},indent=2))
