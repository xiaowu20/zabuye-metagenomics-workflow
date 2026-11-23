"""Connect per-library profiles and the labelled genome catalogue to study-level analyses."""
import argparse
import csv
from pathlib import Path
import shutil
import sys
from primary import GROUPS, read_samples, require, run

HERE = Path(__file__).resolve().parents[1]


def link(src, dst):
    src = src.resolve(strict=True)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.resolve() == src:
            return
        raise FileExistsError(dst)
    dst.symlink_to(src, target_is_directory=src.is_dir())


def catalogue(work, target):
    import pandas as pd
    ids = [p.stem for p in sorted((work/"dereplicate/project/run/dereplicated_genomes").glob("*.fa"))]
    if not ids:
        raise ValueError("Missing dRep representatives")
    q = pd.read_csv(work/"quality/project/checkm2/quality_report.tsv", sep="\t").set_index("Name", verify_integrity=True).loc[ids]
    files = sorted((work/"gtdb/project/r226/classify").glob("*.summary.tsv"))
    tax = pd.concat([pd.read_csv(p, sep="\t") for p in files]).set_index("user_genome", verify_integrity=True)
    rows = []
    for mag, row in q.iterrows():
        group = mag.split("__", 1)[0]
        if group not in GROUPS:
            raise ValueError("Unknown assembly of " + mag)
        lineage = dict(item.split("__", 1) for item in tax.at[mag, "classification"].split(";") if "__" in item)
        rows.append((mag, group, row.Completeness, row.Contamination,
                     "HQ" if row.Completeness >= 90 and row.Contamination < 5 else "MQ",
                     lineage.get("p") or "Unclassified", lineage.get("g") or "Unclassified"))
    columns = ["mag_id","source_assembly","completeness","contamination","quality_tier","gtdb_phylum","gtdb_genus"]
    pd.DataFrame(rows, columns=columns).to_csv(target, sep="\t", index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["taxonomy","humann","catalogue","reads","reference"])
    for key in ("samples","work","project"):
        p.add_argument("--"+key, type=Path, required=True)
    a=p.parse_args()
    rows=read_samples(a.samples)
    w,base=a.work.resolve(),a.project.resolve()
    if a.stage=="reads":
        for r in rows:
            for mate in (1,2):
                sample=r["sample_id"]
                link(w/"trim"/sample/f"R{mate}.fastq.gz",base/"01_fastp"/f"{sample}.R{mate}.clean.fastq.gz")
    elif a.stage=="taxonomy":
        tax=base/"dual_species_reanalysis_20260819"
        for r in rows:
            sample=r["sample_id"]
            for rank,folder in (("pluspf","pluspf"),("regional","regional_residual")):
                for ending in ("report","S.tsv","G.tsv"):
                    dst=ending.replace(".tsv",".bracken.tsv")
                    link(w/"classify"/sample/f"{rank}.{ending}",tax/"classification"/folder/f"{sample}.{dst}")
        sample_list=tax/"input/sample_list.txt"
        sample_list.parent.mkdir(parents=True,exist_ok=True)
        sample_list.write_text("".join(r["sample_id"]+"\n" for r in rows))
        run(sys.executable,HERE/"native/taxonomy/03_merge_dual_bracken.py","--run-root",tax,
            "--sample-list",sample_list,"--output-dir",tax/"data")
    elif a.stage=="humann":
        root=w/"community_ko"
        inputs=root/"by_library"
        for r in rows:
            sample=r["sample_id"]
            for feature in ("genefamilies","pathabundance","pathcoverage"):
                name=f"{sample}_{feature}.tsv"
                link(w/"humann"/sample/"profile"/name,inputs/name)
        for feature in ("genefamilies","pathabundance","pathcoverage"):
            run("humann_join_tables","-s","--input",inputs,"--output",root/f"{feature}_joined.tsv","--file_name",feature)
        run("humann_regroup_table","--input",root/"genefamilies_joined.tsv","--groups","uniref90_ko","--output",root/"KO.tsv")
        run("humann_renorm_table","--input",root/"KO.tsv","--output",root/"KO_cpm.tsv","--units","cpm")
        run("humann_split_stratified_table","--input",root/"KO_cpm.tsv","--output",root/"strata")
    elif a.stage=="catalogue":
        out=w/"catalogue"
        out.mkdir(parents=True,exist_ok=True)
        catalogue(w,out/"MAG_catalogue.tsv")
    elif a.stage=="reference":
        link(w/"catalogue/MAG_catalogue.tsv",base/"packages/pkg_D_mag_qc/mag_catalogue_414_unified.tsv")
        link(w/"dereplicate/project/run/dereplicated_genomes",base/"09_drep/dereplicated_genomes")
        link(w/"quality/project/bins",base/"06_bins/all_dastool_bins")
        link(w/"quality/project/checkm2",base/"07_checkm2")
        for group in GROUPS:
            link(w/"coassemble"/group/"assembly/final.contigs.fa",base/"05_assembly"/group/f"{group}.contigs.fa")
            for algorithm in ("metabat2","maxbin2"):
                link(w/"bin"/group/algorithm,base/"06_bins"/group/algorithm)
            for r in (r for r in rows if r["group"]==group):
                sample=r["sample_id"]
                for end in ("bam","bam.bai"):
                    link(w/"bin"/group/"bam"/f"{sample}.{end}",base/"06_bins"/group/"bam"/f"{sample}.sorted.{end}")


if __name__=="__main__":
    main()
