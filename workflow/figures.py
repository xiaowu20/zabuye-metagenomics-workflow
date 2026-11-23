"""Redraw figures from local source tables. Existing artwork is never used."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from stage import REPO,SUBDIR,digest

GROUPS={"workbooks":["1","5","S1","S2","S3","S4","S5","S11","S14","S16"],
        "traits":["4","S6","S7","S9","S10","S12","S13","S15","S19"],
        "community":["2"],"carriers":["3"],"display":["S8","S17","S18"],
        "coverage":["S20"],"assembly":["S21","S22","S23"],"tree":["S24"],
        "covariates":["S25"],"matches":["S26"],"fcs_quality":["S27"],"fcs_inference":["S28"]}
R_PLOTS={"coverage":("nonpareil.R",["R1_Nonpareil14"]),
         "assembly":("plot_assembly_sensitivity.R",["R1_Catalogue_recovery","R1_Assembly_ecological_sensitivity","R1_Module_direction_sensitivity"]),
         "tree":("plot_phylogenetic_sensitivity.R",["R1_Phylogenetic_sensitivity"]),
         "covariates":("plot_coverage_sensitivity.R",["R1_Coverage_model_sensitivity"]),
         "matches":("plot_controlled_population_comparisons.R",["R1_Controlled_population_comparisons"]),
         "fcs_quality":("fcs_quality_figures.R",["R1_FCS_quality_followup"]),
         "fcs_inference":("fcs_inference_figures.R",["R1_FCS_inference_sensitivity"])}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("books","local-inputs","project","out"):p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--group",choices=["all","paper",*GROUPS],default="all",
                   help="paper uses the supplements and skips the two author-held density displays")
    p.add_argument("--rscript",default="Rscript")
    a=p.parse_args();books=a.books.resolve();data=a.local_inputs.resolve();out=a.out.resolve()
    root=a.project.resolve()/SUBDIR
    if a.group in ("all","paper") or a.group in R_PLOTS:
        if not (root/"workspace.json").is_file():p.error("Stage S5 before rendering the revision figures.")
        if (a.group=="all" or a.group.startswith("fcs_")) and "s6" not in json.loads((root/"workspace.json").read_text())["sources"]:
            p.error("Figures S27 and S28 require a workspace staged with --s6.")
    out.mkdir(parents=True,exist_ok=True)
    groups=list(GROUPS) if a.group in ("all","paper") else [a.group]
    if a.group=="paper" and "s6" not in json.loads((root/"workspace.json").read_text())["sources"]:
        groups=[g for g in groups if not g.startswith("fcs_")]
        print("S6 was not supplied. Figures S27 and S28 are omitted.")
    for g in groups:
        ids=GROUPS[g];job=out/"runs"/g
        if job.exists():p.error(f"Existing plot run {job}; choose a new output or move the previous run.")
        job.mkdir(parents=True)
        if g=="workbooks":cmd=[sys.executable,REPO/"plots/workbook_panels.py",books,job];stems=["Figure_"+i for i in ids];loc=job
        elif g=="traits":
            cmd=[sys.executable,REPO/"plots/additional.py",books,data,job]
            if a.group=="paper":
                ids=[i for i in ids if i not in ("4","S9")]
                cmd.append("--paper")
            stems=["Figure_"+i for i in ids];loc=job
        elif g=="community":
            cmd=[sys.executable,REPO/"plots/figure2.py",data/"figure2",job];stems=["Figure_2_top40_occurrence_reordered"];loc=job
        elif g=="carriers":
            shutil.copytree(data/"figure3/source_data",job/"source_data")
            cmd=[a.rscript,"--vanilla",REPO/"plots/figure3.R",job];stems=["Figure_3_full12"];loc=job/"figures"
        elif g=="display":
            cmd=[sys.executable,REPO/"plots/display_panels.py",data/"si_display",job]
            stems=["Figure_"+i+"_revised" for i in ids];loc=job
        else:
            script,stems=R_PLOTS[g]
            body=(root/"code"/script if g.startswith("fcs_") else REPO/("plots" if g=="coverage" else "native/revision")/script).read_text(encoding="utf-8")
            # Redirect plot outputs only. Keep the original source tables read-only.
            import re
            body=re.sub(r'out\s*<-\s*file.path\(root,\s*["\x27]figures["\x27]\)',"out<-commandArgs(trailingOnly=TRUE)[2]",body)
            body=body.replace('fig<-file.path(root,"figures")','fig<-commandArgs(trailingOnly=TRUE)[2]')
            body=body.replace('file.path(root,"figures",paste0(', 'file.path(commandArgs(trailingOnly=TRUE)[2],paste0(')
            body=re.sub(r'out<-file.path\(follow,"(?:inference|quality)_figures"\);stopifnot\(!dir.exists\(out\)\);dir.create\(out\)',"out<-commandArgs(trailingOnly=TRUE)[2]",body)
            scriptpath=job/script;scriptpath.write_text(body,encoding="utf-8")
            cmd=[a.rscript,"--vanilla",scriptpath,root,job];loc=job
        subprocess.run([str(x) for x in cmd],check=True)
        for fig,stem in zip(ids,stems,strict=True):
            copied=[]
            for ext in (".pdf",".svg",".png"):
                src=loc/(stem+ext)
                if g=="carriers" and ext==".png" and not src.exists():src=job/"qa/Figure_3_preview.png"
                if src.is_file() and src.stat().st_size>100:
                    dest=out/("Figure_"+fig+ext);shutil.copy2(src,dest);copied.append(dest)
            if len(copied)!=3:raise RuntimeError("Expected PDF, SVG and PNG for Figure "+fig)
            print("Redrew",fig,flush=True)
        manifest={p.name:digest(p) for p in sorted(out.glob("Figure_*")) if p.is_file()}
        (out/"checksums.json").write_text(json.dumps(manifest,indent=2))


if __name__=="__main__":main()
