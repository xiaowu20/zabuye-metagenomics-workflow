"""Summarise competitive CoverM recruitment in one fixed MAG catalogue."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from primary import read_samples


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("work","samples","scores","complete","out"):p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--breadth",type=float,default=.50)
    p.add_argument("--catalogue",type=Path)
    p.add_argument("--source",default="all",choices=("all","ZBY_water"))
    a=p.parse_args()
    if not 0<=a.breadth<=1:p.error("Breadth is a fraction in [0, 1].")
    score=pd.read_csv(a.scores,sep="\t",index_col=0);complete=pd.read_csv(a.complete,sep="\t",index_col=0)
    if not all(frame.index.is_unique and frame.columns.is_unique for frame in (score,complete)) or set(score.index)!=set(complete.index) or set(score.columns)!=set(complete.columns):
        raise ValueError("Module matrices must use identical unique MAG and module IDs.")
    complete=complete.loc[score.index,score.columns]
    if not np.isfinite(score.to_numpy()).all() or not complete.isin([0,1]).all().all():
        raise ValueError("Scores must be finite and strict carrier states must be binary.")
    if a.source!="all":
        if not a.catalogue:p.error("--source requires --catalogue")
        catalogue=pd.read_csv(a.catalogue,sep="\t").set_index("mag_id",verify_integrity=True)
        chosen=catalogue.index[catalogue.source_assembly==a.source]
        score,complete=score.loc[chosen],complete.loc[chosen]
    matrices={k:{} for k in ("MAG","GMP","SCA","DCP")}; qc=[]
    for row in read_samples(a.samples):
        if a.source=="ZBY_water" and row["group"]!="ZBY_water":continue
        sample=row["sample_id"]
        tab=pd.read_csv(a.work/"recruit"/sample/"coverm.tsv",sep="\t",index_col=0)
        def column(suffix):
            names=[c for c in tab.columns if c.lower().endswith(suffix)]
            if len(names)!=1:raise ValueError(f"{sample}: expected one {suffix} column, found {names}")
            return pd.to_numeric(tab.loc[score.index,names[0]],errors="raise")
        means=column(" mean");breadth=column(" covered fraction");raw=column(" relative abundance (%)")
        if not (np.isfinite(means).all() and np.isfinite(raw).all() and (means>=0).all()
                and raw.between(0,100).all() and breadth.between(0,1).all()):
            raise ValueError(f"Invalid CoverM input for {sample}")
        detect=(means>=1)&(breadth>=a.breadth)
        abundance=raw.where(detect,0)
        if not detect.any() or abundance.sum()<=0:raise ValueError(f"No supported within-catalogue denominator for {sample}")
        abundance/=abundance.sum()
        matrices["MAG"][sample]=abundance;matrices["GMP"][sample]=abundance@score
        matrices["SCA"][sample]=abundance@complete
        matrices["DCP"][sample]=detect.astype(float)@complete/detect.sum()
        qc.append((sample,int(detect.sum()),raw.sum()))
    a.out.mkdir(parents=True,exist_ok=True)
    for key,records in matrices.items():pd.DataFrame(records).T.rename_axis("sample").to_csv(a.out/(key+".tsv"),sep="\t")
    pd.DataFrame(qc,columns=["sample","detected_MAGs","prethreshold_catalogue_capture_pct"]).to_csv(a.out/"QC.tsv",sep="\t",index=False)


if __name__=="__main__":main()
