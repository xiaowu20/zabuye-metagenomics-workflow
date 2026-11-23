"""Run FCS-GX on supplied analytical MAGs with their assigned NCBI taxids."""
import argparse
import csv
import hashlib
from pathlib import Path
import subprocess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,required=True,help="TSV with mag_id, fasta, taxid and sha256")
    p.add_argument("--database",type=Path,required=True)
    p.add_argument("--runner",type=Path,required=True,help="NCBI fcs.py from the matching release")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    with a.manifest.open(encoding="utf-8-sig") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    if not rows or len({r["mag_id"] for r in rows})!=len(rows):raise ValueError("Expected unique MAG identifiers")
    for r in rows:
        ident=r["mag_id"]
        if not ident or ident in (".","..") or Path(ident).name!=ident:raise ValueError("Invalid identifier")
        genome=Path(r["fasta"])
        if not genome.is_absolute():genome=a.manifest.resolve().parent/genome
        with genome.open("rb") as f:actual=hashlib.file_digest(f,"sha256").hexdigest()
        if actual!=r["sha256"] or int(r["taxid"])<=0:raise ValueError(f"Check identity or taxid for {ident}")
        out=a.out/ident
        if out.exists():raise FileExistsError(f"Use a new report directory: {out}")
        out.mkdir(parents=True)
        subprocess.run([str(a.runner.resolve()),"screen","genome","--fasta",str(genome.resolve()),"--tax-id",r["taxid"],
                        "--gx-db",str(a.database.resolve()),"--out-dir",str(out.resolve())],check=True)


if __name__=="__main__":
    main()
