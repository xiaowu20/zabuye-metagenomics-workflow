"""Run the archived Nonpareil settings on one paired library's R1 mate."""
import argparse
import csv
from pathlib import Path
import subprocess
from stage import SUBDIR


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("sample")
    p.add_argument("--project",type=Path,required=True)
    p.add_argument("--threads",type=int,default=8)
    a=p.parse_args()
    root=a.project.resolve()/SUBDIR
    with (root/"inputs/paired14.tsv").open() as f:
        reads={r["sample"]:Path(r["r1"]) for r in csv.DictReader(f,delimiter="\t")}
    if a.sample not in reads or not reads[a.sample].is_file():p.error("Missing paired R1 input")
    out=root/"nonpareil355"/a.sample
    if (out/"COMPLETE.txt").exists():return
    out.mkdir(parents=True,exist_ok=True)
    subprocess.run(["nonpareil","-s",str(reads[a.sample]),"-T","kmer","-f","fastq","-b",a.sample,
                    "-t",str(a.threads),"-R","24000"],cwd=out,check=True)
    if not (out/(a.sample+".npo")).stat().st_size:raise ValueError("Empty Nonpareil curve")
    (out/"COMPLETE.txt").write_text("ok\n")


if __name__=="__main__":
    main()
