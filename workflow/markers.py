"""Confirm potassium markers, audit neighbourhoods, and reconstruct the K03498 marker tree."""
import argparse
import csv
from itertools import product
from pathlib import Path
import subprocess
from Bio import SeqIO


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("proteins","gff-dir","kofam-hits","out"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--pfam-library",type=Path,help="Full Pfam-A HMM library")
    p.add_argument("--pfam-candidates",type=Path,help="Protein FASTA for the archived Pfam assessment subset")
    p.add_argument("--threads",type=int,default=16)
    p.add_argument("--tree",action="store_true")
    a=p.parse_args()
    if bool(a.pfam_library)!=bool(a.pfam_candidates):p.error("Supply both Pfam inputs together")
    hits={}
    with a.kofam_hits.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            r=line.split()
            if len(r)>=4 and r[0]=="*" and r[1].startswith("ZBY_water__") and r[2] in {"K03498","K03499"}:
                hits.setdefault(r[1],set()).add(r[2])
    if not hits:raise ValueError("No confirmed ZBY-water markers in KofamScan detail-tsv")
    features={}
    for path in sorted(a.gff_dir.glob("ZBY_water__*.gff")):
        with path.open() as f:
            for line in f:
                if line.startswith("#"):continue
                fields=line.rstrip().split("\t")
                if len(fields)!=9 or fields[2]!="CDS":continue
                attributes=dict(pair.split("=",1) for pair in fields[8].split(";") if "=" in pair)
                ordinal=int(attributes["ID"].rsplit("_",1)[1])
                identifier=f"{fields[0]}_{ordinal}"
                if identifier in hits:
                    features[identifier]=(path.stem,fields[0],ordinal,int(fields[3]),int(fields[4]),fields[6])
    if set(features)!=set(hits):raise ValueError("KOfam IDs do not match the source-qualified Prodigal GFFs")
    a.out.mkdir(parents=True,exist_ok=False)
    by_contig={}
    best={}
    for rec in SeqIO.parse(a.proteins,"fasta"):
        if rec.id not in hits:continue
        mag,contig,*_=features[rec.id]
        sets=by_contig.setdefault((mag,contig),{"K03498":[],"K03499":[]})
        for k in hits[rec.id]:sets[k].append(rec.id)
        if "K03498" in hits[rec.id] and (mag not in best or len(rec.seq)>len(best[mag].seq)):
            best[mag]=rec
    if {x for v in by_contig.values() for ids in v.values() for x in ids}!=set(hits):
        raise ValueError("Missing marker sequences")
    with (a.out/"neighbourhood_pairs.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["mag","contig","K03498","K03499","intervening_CDS","codirectional","gap_bp"])
        for (mag,contig),sets in sorted(by_contig.items()):
            for left,right in product(sets["K03498"],sets["K03499"]):
                if left==right:continue
                first,second=sorted((features[left],features[right]),key=lambda row:row[3])
                w.writerow([mag,contig,left,right,abs(first[2]-second[2])-1,first[5]==second[5],second[3]-first[4]-1])
    with (a.out/"marker_carriers.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(["mag","protein","KO"])
        for protein,ks in sorted(hits.items()):
            for ko in sorted(ks):w.writerow([features[protein][0],protein,ko])
    fasta=a.out/"K03498_longest_per_MAG.faa"
    with fasta.open("w") as f:
        for mag,rec in sorted(best.items()):f.write(f">{mag}\n{rec.seq}\n")
    if a.pfam_library:
        subprocess.run(["hmmsearch","--cpu",str(a.threads),"-E","1e-5","--tblout",str(a.out/"Pfam_full_sequence.tbl"),
                        "-o",str(a.out/"Pfam.log"),str(a.pfam_library),str(a.pfam_candidates)],check=True)
    if a.tree:
        alignment=a.out/"K03498.aligned.faa"
        with alignment.open("w") as f:subprocess.run(["mafft","--auto","--thread",str(a.threads),str(fasta)],stdout=f,check=True)
        subprocess.run(["iqtree","-s",str(alignment),"-m","MFP","-bb","1000","-alrt","1000",
                        "-nt",str(a.threads),"-pre",str(a.out/"K03498_iqtree")],check=True)


if __name__=="__main__":
    main()
