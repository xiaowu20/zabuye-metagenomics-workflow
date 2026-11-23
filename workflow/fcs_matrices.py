"""FCS stress tests using frozen original mapping and documented denominators."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess
import numpy as np
import pandas as pd

import argparse
from stage import REPO, SUBDIR

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project", type=Path, required=True, help="Project staged with S5 and S6")
parser.add_argument("--out", type=Path, required=True, help="New scenario matrix and model directory")
parser.add_argument("--fit", action="store_true", help="Also run the exact tests and station bootstraps")
parser.add_argument("--rscript", default="Rscript")
args = parser.parse_args()
ROOT = args.project.resolve() / SUBDIR
FOLLOW = ROOT / "deposition_qc/followup_20261112"
REVIEW = FOLLOW / "scientific_review"
if not (REVIEW / "scientific_review_complete.json").is_file():
    parser.error("Stage the S6 reviewed inputs before running the scenario analysis")
OUT = args.out.resolve()
OUT.mkdir(parents=True, exist_ok=False)
for name in ['inputs','code','catalogues']:(OUT/name).mkdir()
hashes={}
def record(p):hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
for rel in ['inputs/paired14.tsv','inputs/original_feature_tests.tsv','code/paired_station_method_source.R','code/module_feature_tests.R']:
    p=ROOT/rel;record(p);shutil.copy2(p,OUT/rel)
meta=pd.read_csv(ROOT/'inputs/paired14.tsv',sep='\t')
original=ROOT/'catalogues/archival_original269'
score=pd.read_csv(original/'annotation/MAG_module_scores.tsv',sep='\t',index_col=0)
complete=pd.read_csv(original/'annotation/MAG_module_complete.tsv',sep='\t',index_col=0)
actions=pd.read_csv(FOLLOW/'tables/FCS_action_calls.tsv',sep='\t')
all_action=set(actions.mag_id)&set(score.index)
exclude=set(actions.loc[actions.action=='EXCLUDE','mag_id'])&set(score.index)
assert score.shape==(269,49) and len(exclude)==4 and len(all_action)==11
candidate_score=pd.read_csv(REVIEW/'EXCLUDE_contigs_module_scores_414.tsv',sep='\t',index_col=0).loc[score.index,score.columns]
candidate_complete=pd.read_csv(REVIEW/'EXCLUDE_contigs_module_complete_414.tsv',sep='\t',index_col=0).loc[score.index,score.columns]
qc=pd.read_csv(original/'responses/sample_QC.tsv',sep='\t')
scenarios={'baseline_original269':(set(),score,complete),
    'without_EXCLUDE_MAGs':(exclude,score,complete),
    'without_all_action_MAGs':(all_action,score,complete),
    'EXCLUDE_contigs_function_only':(set(),candidate_score,candidate_complete)}
qa=[];impact=[]
for group,(drop,gs,gc) in scenarios.items():
    print('Preparing matrices: '+group,flush=True)
    dest=OUT/'catalogues'/group
    for rel in ['responses','statistics','annotation']:(dest/rel).mkdir(parents=True,exist_ok=False)
    keep=[m for m in score.index if m not in drop];gs=gs.loc[keep];gc=gc.loc[keep]
    gs.to_csv(dest/'annotation/MAG_module_scores.tsv',sep='\t');gc.to_csv(dest/'annotation/MAG_module_complete.tsv',sep='\t')
    new_qc=[]
    for breadth in [10,25,50]:
        path=original/'responses'/('MAG_breadth'+str(breadth)+'.tsv');record(path)
        base=pd.read_csv(path,sep='\t',index_col=0).loc[meta['sample'],score.index]
        assert np.allclose(base.sum(axis=1),1,atol=1e-12,rtol=0)
        q=qc.loc[qc.breadth_cutoff==breadth/100].set_index('sample_id').loc[base.index]
        assert ((base>0).sum(axis=1).values==q.n_detected_mags.values).all()
        retained=base[keep];retained_sum=retained.sum(axis=1);assert (retained_sum>0).all()
        a=retained.div(retained_sum,axis=0);detected=(a>0);n=detected.sum(axis=1)
        matrices={'MAG':a,'GMP':pd.DataFrame(a.values@gs.values,index=a.index,columns=gs.columns),
            'SCA':pd.DataFrame(a.values@gc.values,index=a.index,columns=gs.columns),
            'DCP':pd.DataFrame(detected.values.astype(float)@gc.values,index=a.index,columns=gs.columns).div(n,axis=0)}
        for response,matrix in matrices.items():
            assert np.isfinite(matrix.values).all()
            matrix.to_csv(dest/'responses'/(response+'_breadth'+str(breadth)+'.tsv'),sep='\t',index_label='sample_id')
            if group=='baseline_original269':
                ref=original/'responses'/(response+'_breadth'+str(breadth)+'.tsv');record(ref)
                reference=pd.read_csv(ref,sep='\t',index_col=0).loc[matrix.index,matrix.columns]
                delta=float(np.max(np.abs(matrix.values-reference.values)));assert delta<1e-12,(response,breadth,delta)
                qa.append(dict(response=response,breadth=breadth/100,maximum_absolute_difference=delta))
        for sample in a.index:
            new_qc.append(dict(sample_id=sample,catalogue=group,breadth_cutoff=breadth/100,
                catalogue_MAGs=len(keep),n_detected_mags=int(n[sample]),
                retained_fraction_of_original_detected_abundance=float(retained_sum[sample]),
                removed_fraction_of_original_detected_abundance=float(1-retained_sum[sample])))
            if breadth==50:impact.append(dict(scenario=group,sample_id=sample,
                removed_detected_MAGs=int(q.loc[sample,'n_detected_mags']-n[sample]),
                removed_detected_abundance_fraction=float(1-retained_sum[sample])))
    pd.DataFrame(new_qc).to_csv(dest/'responses/sample_QC.tsv',sep='\t',index=False)
    (dest/'definition.json').write_text(json.dumps(dict(removed_MAGs=sorted(drop),retained_MAGs=len(keep),
        source_mapping='Frozen original 269-MAG mapping; reads are not reassigned',
        interpretation='All-action exclusion is a stress test, not a curation instruction; function-only keeps the original abundances and detection values.'),indent=2))
pd.DataFrame(qa).to_csv(OUT/'baseline_matrix_reproduction_QA.tsv',sep='\t',index=False)
pd.DataFrame(impact).to_csv(OUT/'primary_mapping_mass_and_detection_impact.tsv',sep='\t',index=False)
record(original/'annotation/MAG_module_scores.tsv');record(original/'annotation/MAG_module_complete.tsv')
record(FOLLOW/'tables/FCS_action_calls.tsv')
record(REVIEW/'EXCLUDE_contigs_module_scores_414.tsv');record(REVIEW/'EXCLUDE_contigs_module_complete_414.tsv')
(OUT/'input_receipt.json').write_text(json.dumps(dict(created_at=datetime.datetime.now().astimezone().isoformat(),
    job_id=os.environ.get('SLURM_JOB_ID'),source_sha256=hashes,scenarios=list(scenarios),
    source_mapping_and_sequences_changed=False),indent=2))
print("Wrote four fixed-reference scenario matrices and baseline checks.", flush=True)
if args.fit:
    subprocess.run([args.rscript, "--vanilla", str(REPO/"native/revision/fcs_fixed_reference_statistics.R"),
                    str(OUT), str(ROOT)], check=True)
    if not all((OUT/"catalogues"/g/"statistics/COMPLETE.txt").is_file() for g in scenarios):
        raise RuntimeError("One or more scenario fits did not complete")
    (OUT/"COMPLETE.json").write_text(json.dumps({"scenarios":list(scenarios)}, indent=2))
