#!/usr/bin/env python3
"""Build a reproducible GTDB-derived taxonomy for the regional MAG library.

TPMC genomes are assigned to their published 95% ANI genome-species OTU.
QXLSG contributes only the 2,742 published dereplicated species representatives,
because the public workbook does not map the other MAGs to representatives.
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd


RANKS = ("phylum", "class", "order", "family", "genus")
UNKNOWN = {"", "unknown", "unclassified", "na", "nan", "none"}


def clean(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def parse_gtdb(lineage: object) -> dict[str, str]:
    prefixes = {
        "d__": "domain", "p__": "phylum", "c__": "class",
        "o__": "order", "f__": "family", "g__": "genus", "s__": "species",
    }
    out = {v: "" for v in prefixes.values()}
    for token in clean(lineage).split(";"):
        token = token.strip()
        for prefix, rank in prefixes.items():
            if token.startswith(prefix):
                out[rank] = clean(token[len(prefix):])
                break
    return out


class Taxonomy:
    def __init__(self) -> None:
        self.next_taxid = 10_000_000
        self.nodes: dict[int, tuple[int, str]] = {
            1: (1, "no rank"),
            131567: (1, "no rank"),
            2: (131567, "superkingdom"),
            2157: (131567, "superkingdom"),
        }
        self.names: dict[int, str] = {
            1: "root", 131567: "cellular organisms",
            2: "Bacteria", 2157: "Archaea",
        }
        self.cache: dict[tuple[int, str, str], int] = {}

    def add(self, parent: int, rank: str, name: str) -> int:
        key = (parent, rank, name)
        if key in self.cache:
            return self.cache[key]
        taxid = self.next_taxid
        self.next_taxid += 1
        self.cache[key] = taxid
        self.nodes[taxid] = (parent, rank)
        self.names[taxid] = name
        return taxid

    def lineage(self, values: dict[str, str]) -> tuple[int, str, dict[str, int]]:
        domain = values.get("domain", "")
        domain_taxid = 2157 if domain.lower().startswith("archaea") else 2
        parent = domain_taxid
        parent_name = self.names[parent]
        rank_taxids: dict[str, int] = {}
        for rank in RANKS:
            value = clean(values.get(rank, ""))
            if value.lower() in UNKNOWN:
                value = f"unclassified {parent_name} {rank}"
            parent = self.add(parent, rank, value)
            parent_name = value
            rank_taxids[rank] = parent
        return parent, parent_name, rank_taxids

    def species(self, genus_taxid: int, name: str) -> int:
        return self.add(genus_taxid, "species", name)

    def write(self, out_dir: Path) -> None:
        out_dir.mkdir(parents=True, exist_ok=True)
        with (out_dir / "nodes.dmp").open("w", encoding="utf-8", newline="\n") as handle:
            for taxid in sorted(self.nodes):
                parent, rank = self.nodes[taxid]
                fields = [taxid, parent, rank, "", 0, 0, 1, 0, 0, 0, 0, 0, ""]
                handle.write("\t|\t".join(map(str, fields)) + "\t|\n")
        with (out_dir / "names.dmp").open("w", encoding="utf-8", newline="\n") as handle:
            for taxid in sorted(self.names):
                handle.write(f"{taxid}\t|\t{self.names[taxid]}\t|\t\t|\tscientific name\t|\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--tpmc-genomes",type=Path,required=True)
    parser.add_argument("--qxlsg-genomes",type=Path,required=True)
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_root = args.output_root.resolve()
    taxonomy_dir = output_root / "taxonomy"
    data_dir = output_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    tpmc_file = input_dir / "TPMC_MAG_summary.xlsx"
    qxlsg_file = input_dir / "QXLSG_mmc2.xlsx"
    tpmc = pd.read_excel(tpmc_file, sheet_name="Table S3")
    genome_col = "Genome ID (ANI 0.95 of OTU representative sequences are marked with *)"
    otu_col = "Genomespecies (dRep OTU_ID, ANI 0.95)"
    lineage_col = "GTDB_classification"
    tpmc["mag_id"] = tpmc[genome_col].map(clean).str.rstrip("*")
    tpmc["is_representative"] = tpmc[genome_col].map(clean).str.endswith("*")
    tpmc["otu_id"] = tpmc[otu_col].map(clean)
    if tpmc["mag_id"].duplicated().any():
        raise ValueError("TPMC MAG identifiers are not unique")

    rep_counts = tpmc.loc[tpmc["is_representative"]].groupby("otu_id").size()
    if not (rep_counts == 1).all() or set(rep_counts.index) != set(tpmc["otu_id"]):
        raise ValueError("TPMC does not contain exactly one marked representative per 95% ANI OTU")
    reps = tpmc.loc[tpmc["is_representative"]].set_index("otu_id")

    qxlsg = pd.read_excel(qxlsg_file, sheet_name="Sheet4", header=1)
    qxlsg.columns = [clean(c) for c in qxlsg.columns]
    required_q = {"Genome", "RawName", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "lineage"}
    if missing := required_q.difference(qxlsg.columns):
        raise ValueError(f"QXLSG Sheet4 missing columns: {sorted(missing)}")
    qxlsg["RawName"] = qxlsg["RawName"].map(clean)
    qxlsg = qxlsg[qxlsg["RawName"] != ""].copy()
    if qxlsg["RawName"].duplicated().any():
        raise ValueError("QXLSG species-representative RawName values are not unique")

    tax = Taxonomy()
    records: list[dict[str, object]] = []

    otu_to_taxid: dict[str, int] = {}
    otu_to_genus: dict[str, str] = {}
    for otu_id in sorted(reps.index):
        rep = reps.loc[otu_id]
        parsed = parse_gtdb(rep[lineage_col])
        genus_taxid, genus_name, _ = tax.lineage(parsed)
        species_name = f"TPMC_OTU_{safe_id(otu_id)}"
        species_taxid = tax.species(genus_taxid, species_name)
        otu_to_taxid[otu_id] = species_taxid
        otu_to_genus[otu_id] = genus_name

    for _, row in tpmc.sort_values(["otu_id", "mag_id"]).iterrows():
        otu_id = row["otu_id"]
        mag_id = row["mag_id"]
        records.append({
            "source": "TPMC",
            "reference_id": mag_id,
            "source_species_cluster": otu_id,
            "species_taxid": otu_to_taxid[otu_id],
            "operational_species": tax.names[otu_to_taxid[otu_id]],
            "genus": otu_to_genus[otu_id],
            "gtdb_lineage": clean(row[lineage_col]),
            "fasta_path_hpc": str(args.tpmc_genomes.resolve()/f"{mag_id}.fa"),
            "is_species_representative": bool(row["is_representative"]),
        })

    for row in qxlsg.sort_values("RawName").to_dict("records"):
        parsed = parse_gtdb(row["lineage"])
        if not parsed["domain"]:
            parsed["domain"] = clean(row["Kingdom"])
        for rank, column in zip(RANKS, ("Phylum", "Class", "Order", "Family", "Genus")):
            if not parsed[rank]:
                parsed[rank] = clean(row[column])
        genus_taxid, genus_name, _ = tax.lineage(parsed)
        species_name = f"QXLSG_{safe_id(clean(row['RawName']))}"
        species_taxid = tax.species(genus_taxid, species_name)
        raw_name = clean(row["RawName"])
        # The deposited co-assembly FASTA files omit the leading "all_" used
        # by Dataset S1; no Sheet3 genome uses a native bare "bin.*" name.
        fasta_name = re.sub(r"^all_bin\.", "bin.", raw_name)
        records.append({
            "source": "QXLSG",
            "reference_id": raw_name,
            "source_species_cluster": clean(row["Genome"]),
            "species_taxid": species_taxid,
            "operational_species": species_name,
            "genus": genus_name,
            "gtdb_lineage": clean(row["lineage"]),
            "fasta_path_hpc": str(args.qxlsg_genomes.resolve()/f"{fasta_name}.fa"),
            "is_species_representative": True,
        })

    mapping = pd.DataFrame.from_records(records)
    mapping.to_csv(data_dir / "regional_reference_map.tsv", sep="\t", index=False)
    species = mapping[["source", "species_taxid", "operational_species", "genus", "source_species_cluster"]].drop_duplicates()
    species.to_csv(data_dir / "regional_operational_species.tsv", sep="\t", index=False)
    tax.write(taxonomy_dir)

    source_counts = Counter(mapping["source"])
    species_counts = species.groupby("source")["species_taxid"].nunique().to_dict()
    summary = pd.DataFrame([
        {"metric": "TPMC reference MAGs", "value": source_counts["TPMC"]},
        {"metric": "TPMC 95% ANI operational species", "value": species_counts.get("TPMC", 0)},
        {"metric": "QXLSG published species representatives", "value": source_counts["QXLSG"]},
        {"metric": "QXLSG operational species", "value": species_counts.get("QXLSG", 0)},
        {"metric": "Regional reference genomes", "value": len(mapping)},
        {"metric": "Regional operational species", "value": mapping["species_taxid"].nunique()},
        {"metric": "Taxonomy nodes", "value": len(tax.nodes)},
    ])
    summary.to_csv(data_dir / "regional_taxonomy_summary.tsv", sep="\t", index=False)

    audit = output_root / "regional_taxonomy_audit.md"
    audit.write_text(
        "# Regional operational-species taxonomy audit\n\n"
        f"- TPMC input: `{tpmc_file}`.\n"
        f"- TPMC MAGs: {source_counts['TPMC']:,}; published 95% ANI OTUs: {species_counts.get('TPMC', 0):,}.\n"
        f"- QXLSG input: `{qxlsg_file}`, Sheet4.\n"
        f"- QXLSG species representatives retained: {source_counts['QXLSG']:,}.\n"
        "- QXLSG nonrepresentative MAGs were not included because the workbook does not provide their representative mapping.\n"
        f"- Total regional operational species: {mapping['species_taxid'].nunique():,}.\n"
        "- Taxonomy below domain follows the supplied GTDB lineage. Missing ranks receive parent-specific unclassified labels.\n"
        "- Operational species labels are source-specific and must not be interpreted as formally named species.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
