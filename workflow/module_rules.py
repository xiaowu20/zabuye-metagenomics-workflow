from __future__ import annotations
import numpy as np
import pandas as pd

def classify(n,row):
    def val(k):
        v=row[k]; return None if pd.isna(v) or v=='' else int(v)
    complete,near,pmin,pmax,fmin,fmax=[val(k) for k in ('complete_min','near_complete_min','partial_min','partial_max','fragment_min','fragment_max')]
    if complete is not None and n>=complete:return 'complete'
    if near is not None and n>=near:return 'near_complete'
    if pmin is not None and pmax is not None and pmin<=n<=pmax:return 'partial'
    if fmin is not None and fmax is not None and fmin<=n<=fmax:return 'fragment'
    return 'absent' if n==0 else 'fragment'

def calculate(presence,defs,members,catalog,empty_member_fix=False):
    # Exactly the historical dictionary/order and missing-value conversion.
    gene_to_ko=dict(zip(members.gene,members.KO))
    rows=[]
    leaf=defs[defs.scoring_function!='aggregate_step_count']
    agg=defs[defs.scoring_function=='aggregate_step_count']
    for mag in presence.index:
        def has(g):
            ko=gene_to_ko.get(g)
            return bool(presence.loc[mag,ko]) if ko in presence.columns else False
        bymod={}
        for _,r in leaf.iterrows():
            mod=r.module_id; sf=r.scoring_function
            def parse(key):
                value=r.get(key,'')
                if empty_member_fix and pd.isna(value):return []
                return [g for g in str(value or '').split(',') if g]
            req=parse('required_members');opt=parse('optional_members')
            policy=r.validation_policy
            tent=not (policy=='none' or pd.isna(policy))
            if sf=='single_gene':
                p=has(req[0] if req else opt[0]); status='complete' if p else 'absent'; n=int(p);d=1;trigger=tent and p
            elif sf=='single_gene_or':
                genes=opt if opt else req;p=any(has(g) for g in genes);status='complete' if p else 'absent';n=sum(has(g) for g in genes);d=None;trigger=tent and p
            elif sf in ('subunit_count_2','subunit_count_3','subunit_count_5','subunit_count_7_special','subunit_count_with_optional'):
                n=sum(has(g) for g in req);d=len(req);status=classify(n,r)
                if mod=='Salt_Mrp_complex' and status=='near_complete' and not(has('mrpA') and has('mrpD')):status='partial'
                trigger=tent and status!='absent'
            elif sf=='nested_and_or':
                g1=has('sucA') and has('sucB');g2=has('korA') and has('korB')
                status='complete' if g1 or g2 else 'absent';n=int(g1)+int(g2);d=2;trigger=False
            elif sf=='logic_expression':
                mxa=has('mxaF') and has('mxaI');xox=has('xoxF');p=mxa or xox
                status='complete' if p else 'absent';n=int(p);d=None;trigger=xox
            elif sf in ('directional_marker','context_only'):
                genes=req+opt;n=sum(has(g) for g in genes);d=len(genes) if genes else None;status='present' if n else 'absent';trigger=False
            else:raise ValueError(sf)
            valid='pending_validation' if tent and trigger else status
            item=dict(mag_id=mag,module_id=mod,n_present=n,denominator=d,potential_status=status,validated_status=valid,validation_policy=policy,output_field=r.output_field)
            rows.append(item);bymod[mod]=item
        for _,r in agg.iterrows():
            children=[c for c in str(r.required_members).split(',') if c]
            npot=sum(bymod[c]['potential_status']=='complete' for c in children)
            nval=sum(bymod[c]['validated_status']=='complete' for c in children)
            rows.append(dict(mag_id=mag,module_id=r.module_id,n_present=npot,denominator=len(children),potential_status=classify(npot,r),validated_status=classify(nval,r),validation_policy=r.validation_policy,output_field=r.output_field))
    df=pd.DataFrame(rows)
    meta=['mag_id','source_assembly','quality_tier','gtdb_phylum','gtdb_genus']
    return df.merge(catalog[meta],on='mag_id',how='left')

def scores(long):
    d=pd.to_numeric(long.denominator,errors='coerce')
    score=np.where(d.notna() & (d>0),long.n_present/d,long.potential_status.isin(['complete','present']).astype(float))
    result=long.assign(score=score,complete=(long.potential_status=='complete').astype(int))
    return (result.pivot(index='mag_id',columns='module_id',values='score'),result.pivot(index='mag_id',columns='module_id',values='complete'))
