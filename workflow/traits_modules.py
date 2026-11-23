"""Build MAG-level traits and module matrices from primary workflow products."""
from pathlib import Path
import argparse
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from module_rules import calculate, scores


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--definitions",type=Path,required=True)
    parser.add_argument("--members",type=Path,required=True)
    parser.add_argument("--catalogue",type=Path,required=True,help="TSV with mag_id,source_assembly,quality_tier,gtdb_phylum,gtdb_genus")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    catalogue=pd.read_csv(args.catalogue,sep="\t")
    defs=pd.read_csv(args.definitions);members=pd.read_csv(args.members)
    if not catalogue.mag_id.is_unique:raise ValueError("Duplicate MAG identifier")
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    proteins=args.work/"proteins/project"
    ids=pd.read_csv(proteins/"protein_to_mag.tsv",sep="\t").set_index("query",verify_integrity=True).mag_id
    counts=pd.DataFrame(0,index=catalogue.mag_id,columns=sorted(members.KO.unique()),dtype=int)
    with (args.work/"annotate/project/catalogue.emapper.annotations").open() as handle:
        for line in handle:
            if line.startswith("#"):continue
            fields=line.rstrip("\n").split("\t")
            mag=ids.at[fields[0]]
            if mag not in counts.index:continue
            for ko in set(fields[11].replace("ko:","").split(",")).intersection(counts.columns):
                counts.at[mag,ko]+=1
    long=calculate(counts>0,defs,members,catalogue,empty_member_fix=True)
    graded,strict=scores(long)
    for name,data in [("panel_KO_gene_counts",counts),("MAG_module_scores",graded),("MAG_module_complete",strict)]:
        data.to_csv(out/(name+".tsv"),sep="\t")
    long.to_csv(out/"module_status_long.tsv",sep="\t",index=False)
    values=defaultdict(list); valid=set("ACDEFGHIKLMNPQRSTVWY")
    with (out/"protein_pI.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(["mag_id","protein_id","length","pI","acidic_fraction","basic_fraction"])
        for rec in SeqIO.parse(proteins/"proteins.faa","fasta"):
            mag=ids.at[rec.id]
            if mag not in counts.index:continue
            seq="".join(c for c in str(rec.seq).upper() if c in valid)
            if len(seq)<20:continue
            protein=ProteinAnalysis(seq)
            pi=float(f"{protein.isoelectric_point():.3f}")
            # Preserve the historical column names. Biopython reports these residues on a 0–100 scale.
            pct=protein.amino_acids_percent
            w.writerow([mag,rec.id,len(seq),f"{pi:.3f}",f"{pct['D']+pct['E']:.4f}",f"{pct['K']+pct['R']+pct['H']:.4f}"])
            values[mag].append(pi)
    pd.DataFrame([(mag,len(v),np.median(v)) for mag,v in sorted(values.items())],
                 columns=["mag_id","n_proteins","median_pI"]).to_csv(out/"MAG_median_pI.tsv",sep="\t",index=False)


if __name__=="__main__":main()
