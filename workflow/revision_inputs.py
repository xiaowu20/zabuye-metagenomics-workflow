"""Materialise the archived baseline and the identical-genome restricted remap."""
import argparse
import csv
from pathlib import Path
from stage import SUBDIR


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project",type=Path,required=True)
    base=p.parse_args().project.resolve()
    root=base/SUBDIR
    pairs=list(csv.DictReader((root/"inputs/paired14.tsv").read_text().splitlines(),delimiter="\t"))
    for row in pairs:
        for mate in ("r1","r2"):
            if not Path(row[mate]).is_file():raise FileNotFoundError(row[mate])
    for name in ("archival_original269","remapped_original269"):
        cat=root/"catalogues"/name
        quality=list(csv.DictReader((cat/"representative_quality.tsv").read_text().splitlines(),delimiter="\t"))
        genomes=[base/"09_drep/dereplicated_genomes"/(row["mag_id"]+".fa") for row in quality]
        if len(set(genomes))!=269 or any(not p.is_file() for p in genomes):
            raise ValueError("Restore the archived water-derived representatives with analytical MAG identifiers")
        (cat/"genomes.txt").write_text("".join(str(p)+"\n" for p in genomes))
        required={r["mag_id"] for r in quality}
        if name.startswith("archival"):
            for row in pairs:
                path=cat/"coverm"/(row["sample"]+".tsv")
                with path.open() as src:
                    records=csv.reader(src,delimiter="\t");next(records)
                    labels={r[0] for r in records}
                if not required|{"unmapped"} <= labels:raise ValueError(f"Incomplete archived recruitment {path}")
                path.with_suffix(".complete").write_text("Verified frozen S5 input\n")
        (cat/"COMPLETE.txt").write_text("Verified reference IDs\n")


if __name__=="__main__":
    main()
