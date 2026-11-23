"""Stage the specific paper supplements for a replay or fresh revision run."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
from zipfile import ZipFile

REPO=Path(__file__).resolve().parents[1]
SUBDIR=Path("Major_revised_FM/R1_20261110")
EXPECTED={"s5":"7cdb44a05f6cbe159146e619f0a3660f9c27d1a0f0e7d327ab41163769df52b5",
          "s6":"0fa127d2ef4737c0dca0e936b32939e8508bffb608d3b589ec3b910dd5208115"}
S5_CODE=set("""analyze_nonpareil.R assembly_task.py assembly_local_scratch.py bin_assembly.py
build_catalogue_responses.py build_codon_alignments.py catalogue_ani.py catalogue_annotate.py
catalogue_coverm.py catalogue_pipeline.py catalogue_statistics.R compare_catalogues.py
compare_phylogenies.R controlled_catalogue_ani.py coverage_sensitivity.R eggnog_local_database.py
module_feature_tests.R paired_station_method_source.R phylogeny_initial.slurm codon_tree_stage.slurm
codon_trees.slurm concoct_original.slurm nonpareil355_array.slurm gtdb251_check.slurm
replay_phylogenetic_models.R summarize_assembly_contigs.py summarize_controlled_catalogues.py
verify_genus_replay.R""".split())
S6_CODE=set("""audit_exact_module_inference.py fcs_candidate_rescreen.py fcs_contig_scientific_review.py
fcs_fixed_reference_sensitivity.py fcs_fixed_reference_statistics.R fcs_inference_figures.R
fcs_paired_checkm2.py fcs_phylogenetic_sensitivity.R fcs_phylogenetic_sensitivity.py
fcs_quality_figures.R fcs_remaining.slurm fcs_revision_screen.slurm prepare_fcs_followup_inputs.py""".split())


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):h.update(data)
    return h.hexdigest()


def portable(text, project, conda):
    home_pattern=r"/dssg/home/[^/\s'\";]+/[^/\s'\";]+"
    root=(project/SUBDIR).as_posix()
    project=project.as_posix();conda=conda.as_posix()
    for template,target in [(home_pattern+"/Metagenome/Major_revised_FM/R1_20261110",root),
                            (home_pattern+"/Metagenome",project),(home_pattern+"/miniconda3",conda)]:
        text=re.sub(template,lambda _:target,text)
    for label,value in [("@REVISION@",root),("@PROJECT@",project),("@CONDA@",conda)]:
        text=text.replace(label,value)
    return "\n".join(l for l in text.splitlines() if not re.match(r"^#SBATCH\s+(?:--(?:account|partition|mail-user|mail-type)|-[Ap]\b)",l))+"\n"


def put(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    if isinstance(data,str):data=data.encode()
    path.write_bytes(data)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ["s5","project"]:p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--s6",type=Path,help="Add the FCS-GX follow-up and Figures S27–S28")
    p.add_argument("--conda",type=Path,default=Path(sys.prefix))
    p.add_argument("--mode",choices=["replay","fresh"],default="replay")
    a=p.parse_args();project=a.project.resolve();root=project/SUBDIR;conda=a.conda.resolve()
    supplied={"s5":a.s5}
    if a.s6:supplied["s6"]=a.s6
    for key in supplied:
        if digest(getattr(a,key))!=EXPECTED[key]:p.error("File differs from the frozen "+key.upper()+" revision archive.")
    if root.exists():p.error("Choose a new workspace. Staging does not overlay an existing run.")
    fresh=a.mode=="fresh"
    if any(c in str(project)+str(conda) for c in ("\n","'","\"")):p.error("Use paths without quotes or newlines for template rendering.")
    if fresh and any(c.isspace() for c in str(project)+str(conda)):p.error("Fresh runs use native shell templates; choose project and conda paths without spaces.")
    follow=root/"deposition_qc/followup_20261112"
    root.mkdir(parents=True)
    for label,archive in supplied.items():
        with ZipFile(archive) as z:
            for member in z.infolist():
                if member.is_dir():continue
                parts=PurePosixPath(member.filename).parts
                if ".." in parts or not parts or PurePosixPath(member.filename).is_absolute():raise ValueError(member.filename)
                if label=="s5":
                    if parts[0]!="Supplementary_Data_S5":continue
                    parts=parts[1:]
                    if not parts:continue
                    dest=root.joinpath(*parts)
                    if parts[0]=="code":
                        if parts[-1] not in S5_CODE:continue
                    elif parts[0]=="inputs":pass
                    elif fresh and parts[0]=="catalogues" and len(parts)>2 and parts[1] in ("archival_original269","remapped_original269"):
                        # Historical all-reference mapping is the comparison baseline. Rerun the restricted remap.
                        if parts[2] not in ("annotation","representative_quality.tsv") and not (parts[1]=="archival_original269" and parts[2]=="coverm"):continue
                    elif not fresh and parts[0] in ("catalogues","results","phylogeny","nonpareil355"):
                        if any(n in parts for n in ("figures","plots")) or dest.suffix.lower() in (".png",".svg",".pdf",".tif"):continue
                    else:continue
                else:
                    dest=follow.joinpath(*parts)
                    if parts[0]=="code":
                        if parts[-1] not in S6_CODE:continue
                        dest=root/"code"/parts[-1]
                    elif not fresh and parts[0] in ("tables","paired_quality_results","reference_context","scientific_review","fixed_reference_sensitivity","phylogenetic_sensitivity"):pass
                    else:continue
                data=z.read(member)
                if label=="s5" and parts==("code","compare_catalogues.py"):
                    match=re.search(rb"(?m)^\s*target='([^'\r\n]+)'",data)
                    if not match:raise ValueError("Missing study focal-carrier selector")
                    put(root/"inputs/focal_carrier.txt",match[1]+b"\n")
                if parts[0]=="code" or (label=="s5" and parts[0]=="inputs" and dest.suffix in (".json",".tsv")):
                    data=portable(data.decode("utf-8-sig"),project,conda)
                put(dest,data)
    for file in (REPO/"native/s5").glob("*.in"):
        name=file.name.removesuffix(".in")
        target=project/"scripts"/name if name=="fasta_to_s2b.py" else root/"code"/name
        put(target,portable(file.read_text(encoding="utf8"),project,conda))
    for file in (REPO/"native/revision").iterdir():
        if file.is_file():put(root/"code"/file.name,portable(file.read_text(encoding="utf-8"),project,conda))
    shutil.copy2(REPO/"workflow/module_rules.py",root/"code/module_scoring_source.py")
    for file in (REPO/"native/taxonomy").iterdir():
        rel=Path("hpc" if file.name.endswith(".in") else "code")/file.name.removesuffix(".in")
        put(project/"dual_species_reanalysis_20260819"/rel,portable(file.read_text(encoding="utf-8"),project,conda))
    for file in (REPO/"native/external").glob("*.in"):
        put(project/"external"/file.name.removesuffix(".in"),portable(file.read_text(encoding="utf-8"),project,conda))
    (root/"figures").mkdir(exist_ok=True)
    for name in ("logs","provenance","results"):
        (root/name).mkdir(exist_ok=True)
    if not fresh and a.s6:
        shutil.copy2(follow/"paired_quality_results/paired_CheckM2_comparison_6.tsv",follow/"tables/paired_CheckM2_comparison_6.tsv")
    put(root/"workspace.json",json.dumps({"mode":a.mode,"sources":{k:EXPECTED[k] for k in supplied}},indent=2))
    print(root)


if __name__=="__main__":main()
