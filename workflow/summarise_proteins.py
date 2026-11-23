"""Summarise rounded protein records by MAG with the archived pI definitions."""
import argparse
import csv
from itertools import groupby
from pathlib import Path
from statistics import mean, median


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    seen=set()
    with a.proteins.open() as src,a.out.open("x",newline="") as out:
        writer=csv.writer(out,delimiter="\t")
        writer.writerow(["mag_id","n_proteins","median_pI","mean_pI","mean_acidic_fraction","mean_basic_fraction","acidic_minus_basic"])
        for mag,records in groupby(csv.DictReader(src,delimiter="\t"),key=lambda r:r["mag_id"]):
            if mag in seen:raise ValueError("Sort the protein rows by mag_id before summarising")
            seen.add(mag)
            pis,acid,basic=zip(*[(float(r["pI"]),float(r["acidic_fraction"]),float(r["basic_fraction"])) for r in records])
            av,bv=mean(acid),mean(basic)
            writer.writerow([mag,len(pis),median(pis),mean(pis),av,bv,av-bv])


if __name__=="__main__":
    main()
