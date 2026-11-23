#!/usr/bin/env python3
"""Stream regional MAG FASTA files into one Kraken2-labelled library."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--missing", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    with args.mapping.open("r", encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    required = {"source", "reference_id", "species_taxid", "fasta_path_hpc"}
    mapping_columns = set(mapping[0]) if mapping else set()
    if missing_cols := required.difference(mapping_columns):
        raise ValueError(f"Missing mapping columns: {sorted(missing_cols)}")

    missing_paths = [r["fasta_path_hpc"] for r in mapping if not Path(r["fasta_path_hpc"]).is_file()]
    args.missing.parent.mkdir(parents=True, exist_ok=True)
    args.missing.write_text("\n".join(missing_paths) + ("\n" if missing_paths else ""), encoding="utf-8")
    if missing_paths:
        raise FileNotFoundError(f"{len(missing_paths)} mapped FASTA files are missing; see {args.missing}")
    if args.validate_only:
        print(f"validated_reference_paths={len(mapping)}")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    genomes = 0
    sequences = 0
    bases = 0
    with args.output.open("w", encoding="utf-8", newline="\n") as out:
        for rec in mapping:
            source = rec["source"]
            reference = rec["reference_id"]
            taxid = rec["species_taxid"]
            seq_index = 0
            with open_text(Path(rec["fasta_path_hpc"])) as handle:
                for line in handle:
                    if line.startswith(">"):
                        seq_index += 1
                        sequences += 1
                        out.write(f">{source}|{reference}|seq{seq_index}|kraken:taxid|{taxid}\n")
                    else:
                        seq = line.strip()
                        bases += len(seq)
                        out.write(seq + "\n")
            if seq_index == 0:
                raise ValueError(f"No sequences found in {rec['fasta_path_hpc']}")
            genomes += 1

    print(f"genomes={genomes}")
    print(f"sequences={sequences}")
    print(f"bases={bases}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
