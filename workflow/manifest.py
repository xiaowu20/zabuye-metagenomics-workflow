"""Join a local SRA run table to the physical library order in Supplementary Data S5."""
import argparse
import csv
import io
import json
import tempfile
from pathlib import Path, PurePosixPath
from zipfile import ZipFile
from primary import GROUPS, read_samples


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s5", required=True, type=Path)
    p.add_argument("--runs", required=True, type=Path, help="TSV or CSV export of the library ID and run accession")
    p.add_argument("--sample-column", default="sample_id")
    p.add_argument("--run-column", default="run")
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    with a.runs.open(encoding="utf-8-sig", newline="") as f:
        header = f.readline()
        f.seek(0)
        records = list(csv.DictReader(f, delimiter="\t" if "\t" in header else ","))
    if not records or not {a.sample_column,a.run_column}.issubset(records[0]):
        p.error("Check the library and run column headers in --runs")
    runs = {r[a.sample_column].strip(): r[a.run_column].strip() for r in records}
    if len(runs) != len(records) or len(set(runs.values())) != len(records):
        raise ValueError("Expected a one-to-one sample-to-run mapping")
    rows = []
    with ZipFile(a.s5) as source:
        paired = source.read("Supplementary_Data_S5/inputs/paired14.tsv").decode("utf-8-sig")
        ids = {row["sample"] for row in csv.DictReader(io.StringIO(paired), delimiter="\t")}
        for group in GROUPS:
            path = f"Supplementary_Data_S5/provenance/original_commands/05_assembly/{group}/options.json"
            options = json.loads(source.read(path))
            for position, (r1, r2) in enumerate(zip(options["pe1"], options["pe2"], strict=True), 1):
                sample, suffix = PurePosixPath(r1).name.split(".R1.", 1)
                if PurePosixPath(r2).name != f"{sample}.R2.{suffix}":
                    raise ValueError(f"Mate mismatch in {group} at position {position}")
                rows.append([sample, runs[sample], group, position, int(sample in ids)])
    if set(r[0] for r in rows) != set(runs) or len(rows) != 35 or sum(r[-1] for r in rows) != 14:
        raise ValueError("Sample IDs do not reconcile with the study design")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.out.parent) as temp:
        candidate=Path(temp)/"samples.tsv"
        with candidate.open("w",encoding="utf-8",newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "run", "group", "assembly_order", "paired"])
            writer.writerows(rows)
        read_samples(candidate)
        with a.out.open("xb") as f:
            f.write(candidate.read_bytes())
    print(f"Wrote {len(rows)} libraries to {a.out}")


if __name__ == "__main__":
    main()
