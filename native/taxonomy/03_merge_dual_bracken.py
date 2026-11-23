#!/usr/bin/env python3
"""Merge PlusPF and residual-read regional Bracken counts without double counting."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


COUNT_CANDIDATES = ("new_est_reads", "new_est_reads ", "estimated_reads")


def read_bracken(path: Path, source: str, rank: str) -> list[dict[str, object]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        count_col = next((c for c in COUNT_CANDIDATES if c in fields), None)
        if count_col is None:
            raise ValueError(f"No Bracken estimated-read column in {path}: {fields}")
        name_col = "name" if "name" in fields else fields[0]
        taxid_col = "taxonomy_id" if "taxonomy_id" in fields else fields[1]
        rows = []
        for row in reader:
            count = float(row[count_col])
            if count <= 0:
                continue
            rows.append({"taxon": row[name_col].strip(), "taxid": row[taxid_col].strip(),
                         "count": count, "source_database": source, "rank": rank})
        return rows


def kraken_root_count(path: Path) -> tuple[int, int]:
    classified = 0
    unclassified = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            rank, taxid = fields[3].strip(), fields[4].strip()
            clade = int(fields[1])
            if rank == "U" or taxid == "0":
                unclassified = clade
            elif rank == "R" or taxid == "1":
                classified = clade
    return classified, unclassified


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--sample-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    samples = [s.strip() for s in args.sample_list.read_text(encoding="utf-8").splitlines() if s.strip()]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    long_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []

    for sample in samples:
        pf_dir = args.run_root / "classification" / "pluspf"
        rg_dir = args.run_root / "classification" / "regional_residual"
        pf_class, pf_unclass = kraken_root_count(pf_dir / f"{sample}.report")
        rg_class, rg_unclass = kraken_root_count(rg_dir / f"{sample}.report")
        if rg_class + rg_unclass != pf_unclass:
            raise ValueError(f"Residual accounting failed for {sample}: {rg_class + rg_unclass} != {pf_unclass}")

        bracken_totals: dict[str, float] = {}
        for rank, suffix in (("species", "S"), ("genus", "G")):
            pluspf_rows = read_bracken(pf_dir / f"{sample}.{suffix}.bracken.tsv", "PlusPF", rank)
            regional_rows = read_bracken(
                rg_dir / f"{sample}.{suffix}.bracken.tsv", "Regional residual", rank
            )
            bracken_totals[f"pluspf_bracken_{rank}_total"] = sum(float(row["count"]) for row in pluspf_rows)
            bracken_totals[f"regional_bracken_{rank}_total"] = sum(float(row["count"]) for row in regional_rows)
            rows = pluspf_rows + regional_rows
            grouped: dict[str, dict[str, object]] = {}
            for row in rows:
                taxon = str(row["taxon"])
                if rank == "species" and row["source_database"] == "Regional residual":
                    taxon = f"REGIONAL::{taxon}"
                if taxon not in grouped:
                    grouped[taxon] = {"taxon": taxon, "taxid": set(), "count": 0.0,
                                      "source_database": set(), "rank": rank}
                grouped[taxon]["taxid"].add(str(row["taxid"]))
                grouped[taxon]["source_database"].add(str(row["source_database"]))
                grouped[taxon]["count"] += float(row["count"])
            total = sum(float(row["count"]) for row in grouped.values())
            for row in grouped.values():
                long_rows.append({"sample_id": sample, "rank": rank, "taxon": row["taxon"],
                                  "taxid": ";".join(sorted(row["taxid"])),
                                  "source_database": "+".join(sorted(row["source_database"])),
                                  "count": row["count"],
                                  "relative_abundance": float(row["count"]) / total if total else 0.0})

        total_fragments = pf_class + pf_unclass
        qc_rows.append({"sample_id": sample, "pluspf_classified": pf_class,
                        "pluspf_unclassified": pf_unclass, "regional_classified_from_residual": rg_class,
                        "unclassified_after_both": rg_unclass, "total_fragments": total_fragments,
                        "pluspf_classified_fraction": pf_class / total_fragments,
                        "regional_increment_fraction_total": rg_class / total_fragments,
                        "dual_classified_fraction": (pf_class + rg_class) / total_fragments,
                        **bracken_totals})

    long_fields = ["sample_id", "rank", "taxon", "taxid", "source_database", "count", "relative_abundance"]
    write_tsv(args.output_dir / "dual_bracken_long.tsv", long_rows, long_fields)
    write_tsv(args.output_dir / "dual_classification_qc.tsv", qc_rows, list(qc_rows[0]))
    for rank in ("species", "genus"):
        rank_rows = [r for r in long_rows if r["rank"] == rank]
        taxa = sorted({str(r["taxon"]) for r in rank_rows})
        counts: dict[tuple[str, str], float] = defaultdict(float)
        rel: dict[tuple[str, str], float] = defaultdict(float)
        for row in rank_rows:
            key = (str(row["taxon"]), str(row["sample_id"]))
            counts[key] += float(row["count"])
            rel[key] += float(row["relative_abundance"])
        for label, values in (("counts", counts), ("relative_abundance", rel)):
            with (args.output_dir / f"dual_{rank}_{label}.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([rank] + samples)
                for taxon in taxa:
                    writer.writerow([taxon] + [values[(taxon, sample)] for sample in samples])
    print(f"samples={len(samples)}")
    print(f"species={len({r['taxon'] for r in long_rows if r['rank'] == 'species'})}")
    print(f"genera={len({r['taxon'] for r in long_rows if r['rank'] == 'genus'})}")


if __name__ == "__main__":
    main()
