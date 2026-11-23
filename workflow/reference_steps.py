"""Run reference-dependent steps on the original input alignments, candidate proteins or external genomes."""
import argparse
import gzip
from pathlib import Path
import shutil
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    stages = p.add_subparsers(dest="stage",required=True)
    archival = stages.add_parser("archival-tree", help="Earlier GTDB r220 archaeal alignment before tip pruning")
    archival.add_argument("--msa",type=Path,required=True,help="gtdbtk.ar53.user_msa.fasta.gz from the earlier collection")
    kofam = stages.add_parser("kofam", help="Score the source-qualified candidate proteins with KOfamScan")
    for key in ("query","profiles","ko-list"):
        kofam.add_argument("--"+key,type=Path,required=True)
    external = stages.add_parser("external-taxonomy", help="Classify all external MAGs with GTDB-Tk 2.5.1 and r226")
    external.add_argument("--genomes",type=Path,required=True)
    for step in (archival,kofam,external):
        step.add_argument("--out",type=Path,required=True)
        step.add_argument("--threads",type=int,default=32)
    args=p.parse_args()
    if args.threads<1:p.error("Expected a positive thread count")
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    if args.stage=="archival-tree":
        alignment=out/"ar53_user_msa.fasta"
        with gzip.open(args.msa,"rb") as source,alignment.open("wb") as dest:
            shutil.copyfileobj(source,dest)
        command=["iqtree","-s",alignment,"-m","LG+G4+F","-bb","1000","-alrt","1000",
                 "-nt",str(args.threads),"-seed","388326","-pre",out/"ar53_iqtree"]
    elif args.stage=="kofam":
        command=["exec_annotation","-f","detail-tsv","-o",out/"kofam_detail_output.tsv",
                 "--cpu",str(args.threads),"-p",args.profiles,"-k",args.ko_list,
                 "--tmp-dir",out/"tmp_kofam","--create-alignment",args.query]
    else:
        command=["gtdbtk","classify_wf","--genome_dir",args.genomes,"--out_dir",out,
                 "--extension","fa","--cpus",str(args.threads),"--skip_ani_screen",
                 "--pplacer_cpus",str(min(16,args.threads))]
    subprocess.run(list(map(str,command)),check=True)


if __name__=="__main__":
    main()
