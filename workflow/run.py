"""Plan or execute the study workflow, one dependent operation at a time."""
import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from primary import GROUPS, read_samples
from stage import REPO, SUBDIR

BRANCHES = ("reads", "taxonomy", "community-ko", "historical", "primary-traits",
            "sensitivity", "trees", "replay", "study-models", "figures")
NEW = ("single_megahit14", "single_metaspades14", "pooled_megahit14",
       "original_ZBY_water_das2_matched", "original_ZBY_water_das3")
ARCHIVE, REMAP = "archival_original269", "remapped_original269"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("branch", choices=BRANCHES)
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--samples", type=Path)
    p.add_argument("--work", type=Path)
    for name in ("books", "local-inputs", "out"):
        p.add_argument("--"+name, type=Path)
    p.add_argument("--assembly-reads", "--provider", type=Path, help="Archived 00_qc_raw FASTQ directory")
    p.add_argument("--assembly-threads", type=int, help="Override --threads for the four historical MEGAHIT jobs")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--execute", action="store_true", help="Default prints commands without launching any compute")
    a = p.parse_args()
    if a.threads < 1 or (a.assembly_threads is not None and a.assembly_threads < 1):
        p.error("Choose a positive thread count")
    base = a.project.resolve()
    if a.branch not in ("replay","study-models","figures"):
        if not a.samples or not a.work:p.error("This branch requires --samples and --work")
        rows=read_samples(a.samples)
        counts = {g:sum(r["group"] == g for r in rows) for g in GROUPS}
        if counts != dict(zip(GROUPS, (21,4,5,5))) or sum(r["paired"]=="1" for r in rows)!=14:
            p.error("The manifest must describe the 35 physical libraries and paired subset")
        work=a.work.resolve()
        common=["--samples",a.samples.resolve(),"--work",work]
        library_ids=[r["sample_id"] for r in rows]
    root = base/SUBDIR
    code=root/"code"
    mode_path=root/"workspace.json"
    if not mode_path.is_file():
        p.error("Stage S5/S6 into --project before planning")
    mode=json.loads(mode_path.read_text())["mode"]
    if a.branch in ("sensitivity","trees") and mode!="fresh":
        p.error("Use a workspace staged with --mode fresh for genome reconstruction")
    if a.branch=="replay" and mode!="replay":
        p.error("Replay uses the frozen catalogue and mapping files staged with --mode replay")
    allocation = os.environ.copy()
    allocation["SLURM_CPUS_PER_TASK"]=str(a.threads)
    ops=0
    def cmd(*argv, index=None):
        nonlocal ops
        command=list(map(str,argv))
        env=allocation.copy()
        prefix=""
        if index is not None:
            env["SLURM_ARRAY_TASK_ID"]=str(index)
            prefix=f"SLURM_ARRAY_TASK_ID={index} "
        print(prefix+shlex.join(command),flush=True)
        if a.execute:
            subprocess.run(command,env=env,check=True)
        ops+=1
    def python(name,*args):
        cmd(sys.executable,REPO/"workflow"/name,*args)
    def primary(stage,units):
        for unit in units:
            threads = (a.assembly_threads or a.threads) if stage == "coassemble" else a.threads
            args=[stage,unit,*common,"--threads",threads]
            if stage=="coassemble":
                args.extend(["--assembly-reads",a.assembly_reads.resolve()])
            python("primary.py",*args)
    def join(stage):
        python("joins.py",stage,*common,"--project",base)
    def r(file,*args):
        cmd(a.rscript,"--vanilla",file,*args)

    if a.branch=="reads":
        for step in ("fetch","trim"):
            primary(step,library_ids)
        join("reads")
    elif a.branch=="taxonomy":
        if not a.local_inputs:p.error("Pass --local-inputs with taxonomy/manifest_with_site.tsv")
        primary("classify",library_ids)
        join("taxonomy")
        tax=base/"dual_species_reanalysis_20260819"
        if a.execute:
            from joins import link
            link(a.local_inputs/"taxonomy/manifest_with_site.tsv",tax/"data/manifest_with_site.tsv")
        r(tax/"code/04_run_dual_species_analysis.R")
    elif a.branch=="community-ko":
        primary("humann",library_ids)
        join("humann")
        if a.local_inputs and a.out:
            r(REPO/"workflow/community_ko.R",a.local_inputs/"community_ko",a.out/"community_ko",work/"community_ko/KO_cpm.tsv")
    elif a.branch=="historical":
        if not a.assembly_reads:p.error("Historical assemblies require --assembly-reads with the archived input FASTQs")
        primary("coassemble",GROUPS)
        primary("bin",GROUPS)
        for step in ("quality","dereplicate","gtdb"):
            primary(step,["project"])
        join("catalogue")
        join("reference")
    elif a.branch=="primary-traits":
        for step in ("proteins","annotate"):
            primary(step,["project"])
        primary("recruit",library_ids)
        traits=work/"traits"
        python("traits_modules.py","--work",work,"--definitions",root/"inputs/module_definitions_corrected.csv",
               "--members",root/"inputs/zby_ko_panel_members_v3.csv","--catalogue",work/"catalogue/MAG_catalogue.tsv","--out",traits)
        for name in ("all","ZBY_water"):
            python("responses.py",*common,"--scores",traits/"MAG_module_scores.tsv","--complete",traits/"MAG_module_complete.tsv",
                   "--catalogue",work/"catalogue/MAG_catalogue.tsv","--source",name,"--out",work/"responses"/name)
    elif a.branch=="sensitivity":
        python("revision_inputs.py","--project",base)
        with (root/"inputs/paired14.tsv").open() as f:
            paired=list(csv.DictReader(f,delimiter="\t"))
        tasks=json.loads((root/"inputs/assembly_tasks.json").read_text())
        if len(tasks)!=29 or {r["sample"] for r in paired}!={r["sample_id"] for r in rows if r["paired"]=="1"}:
            p.error("Task and library indices do not reconcile with S5")
        for index,_ in enumerate(tasks):
            cmd(sys.executable,code/"assembly_task.py",index)
            cmd(sys.executable,code/"bin_assembly.py",index)
        for index in range(len(GROUPS)):
            cmd("bash",code/"concoct_original.slurm",index=index)
        for item in paired:
            python("nonpareil.py",item["sample"],"--project",base,"--threads",a.threads)
        r(code/"analyze_nonpareil.R",root)
        for group in (*NEW,REMAP,ARCHIVE):
            if group in NEW:
                cmd(sys.executable,code/"catalogue_pipeline.py",group)
                cmd(sys.executable,code/"catalogue_annotate.py",group)
                cmd(sys.executable,code/"catalogue_ani.py",group)
            if group!=ARCHIVE:
                for index,_ in enumerate(paired):
                    cmd(sys.executable,code/"catalogue_coverm.py",group,index=index)
            cmd(sys.executable,code/"build_catalogue_responses.py",group)
            r(code/"catalogue_statistics.R",root,group)
            r(code/"coverage_sensitivity.R",root,group)
        r(code/"coverage_sensitivity.R",root,"read_species")
        cmd(sys.executable,code/"controlled_catalogue_ani.py")
        cmd(sys.executable,code/"compare_catalogues.py")
    elif a.branch=="trees":
        for script in ("phylogeny_initial.slurm","gtdb251_check.slurm"):
            cmd("bash",code/script)
        cmd(sys.executable,code/"build_codon_alignments.py")
        for kind in ("AA","NT"):
            cmd("bash",code/"codon_tree_stage.slurm",kind)
        r(code/"replay_phylogenetic_models.R",root)
        r(code/"compare_phylogenies.R",root)
    elif a.branch=="replay":
        if bool(a.books) != bool(a.out):
            p.error("Supply both --books and --out to add the S4 external-lake model, or omit both for S5 alone")
        if a.books and not (a.books/"Supplementary_Data_S4.xlsx").is_file():
            p.error("The requested external-lake model needs Supplementary_Data_S4.xlsx in --books")
        r(code/"replay_phylogenetic_models.R",root)
        r(code/"compare_phylogenies.R",root)
        for group in (ARCHIVE,REMAP,*NEW):
            r(code/"catalogue_statistics.R",root,group)
            r(code/"coverage_sensitivity.R",root,group)
        r(code/"coverage_sensitivity.R",root,"read_species")
        if a.books:
            r(REPO/"workflow/external_models.R",a.books/"Supplementary_Data_S4.xlsx",a.out/"external")
    elif a.branch=="study-models":
        if not a.books or not a.local_inputs or not a.out:p.error("Supply --books, --local-inputs and a new --out")
        src=a.local_inputs
        for script,folder in (("paired_models","ions"),("incidence","ions"),("ion_models","ions"),
                              ("dom_models","ions"),("crosslake","ions"),("community_ko","community_ko"),
                              ("proteome_models","proteome")):
            r(REPO/"workflow"/(script+".R"),src/folder,a.out/script)
        r(REPO/"workflow/dom_engines.R",src/"dom_engines",a.out/"dom_engines")
        r(REPO/"workflow/external_models.R",a.books/"Supplementary_Data_S4.xlsx",a.out/"external")
        python("carrier_attribution.py",src/"carrier",a.out/"carrier")
        python("figure3_tables.py","--carrier",a.out/"carrier","--models",src/"figure3/models",
               "--palette",src/"si_display/palette.json","--out",a.out/"figure3/source_data")
    elif a.branch=="figures":
        if not a.books or not a.local_inputs or not a.out:p.error("Supply --books, --local-inputs and --out")
        python("figures.py","--books",a.books,"--local-inputs",a.local_inputs,"--project",base,"--out",a.out,"--rscript",a.rscript)
    print(f"{'Executed' if a.execute else 'Planned'} {ops} operations.",file=sys.stderr)


if __name__=="__main__":
    main()
