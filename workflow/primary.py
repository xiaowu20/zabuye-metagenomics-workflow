"""Run individual Linux compute stages for the primary Zabuye analysis."""
import argparse
import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


GROUPS = ("ZBY_water", "ZBY_sed", "DXC_water", "DXC_sed")
SAMPLE_STEPS = ("fetch", "trim", "classify", "humann", "recruit")
GROUP_STEPS = ("coassemble", "bin")
PROJECT_STEPS = ("quality", "dereplicate", "gtdb", "proteins", "annotate")


def run(*command, **kwargs):
    cmd = [str(x) for x in command]
    print("+", subprocess.list2cmdline(cmd), flush=True)
    subprocess.run(cmd, check=True, **kwargs)


def require(paths):
    for p in paths:
        if not p.is_file() or not p.stat().st_size:
            raise FileNotFoundError(p)


def database(variable):
    value = os.environ.get(variable)
    if not value or not Path(value).exists():
        raise ValueError("Set "+variable+" to the installed database path.")
    return Path(value).resolve()


def read_samples(path):
    import re
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows or len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Expected a nonempty manifest of unique physical libraries.")
    for row in rows:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row["sample_id"]) or row["group"] not in GROUPS:
            raise ValueError("Invalid sample or assembly group")
        if not re.fullmatch(r"(SRR|ERR|DRR)\d+", row["run"]):
            raise ValueError("Expected a public run accession for each library.")
    if len({(r["group"],int(r["assembly_order"])) for r in rows})!=len(rows):
        raise ValueError("Duplicate assembly input position within a group.")
    for group in GROUPS:
        order=sorted(int(row["assembly_order"]) for row in rows if row["group"]==group)
        if order!=list(range(1,len(order)+1)):
            raise ValueError("Assembly order must start at one and be contiguous in "+group)
    if len({r["run"] for r in rows})!=len(rows) or any(r.get("paired") not in ("0","1") for r in rows):
        raise ValueError("Each library needs one unique run accession and a 0/1 paired flag.")
    return rows


def fastqs(root, step, sample):
    pair = [root/step/sample/f"R{mate}.fastq.gz" for mate in (1,2)]
    require(pair)
    return pair


def fasta_ids(path):
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].split()[0]


def fasta_map(directory, extension, out):
    with out.open("w") as f:
        for genome in sorted(directory.glob("*."+extension)):
            for contig in fasta_ids(genome):
                print(contig,genome.stem,sep="\t",file=f)
    require([out])


def run_step(a, rows, out):
    w, key, threads = a.work, a.unit, a.threads
    stage = a.stage
    genomes = w/"dereplicate"/"project"/"run"/"dereplicated_genomes"
    if stage=="fetch":
        acc=next(r["run"] for r in rows if r["sample_id"]==key)
        with tempfile.TemporaryDirectory(dir=out) as temp:
            tmp=Path(temp)
            run("prefetch",acc,"--max-size","u","--output-directory",tmp)
            run("fasterq-dump","--split-files","-e",threads,"-t",tmp,"-O",tmp,tmp/acc)
            for mate in (1,2):
                raw=tmp/f"{acc}_{mate}.fastq"
                require([raw])
                with raw.open("rb") as src,gzip.open(out/f"R{mate}.fastq.gz","wb") as dst:
                    shutil.copyfileobj(src,dst,1024*1024)
    elif stage=="trim":
        r1,r2=fastqs(w,"fetch",key)
        run("fastp","-i",r1,"-I",r2,"-o",out/"R1.fastq.gz","-O",out/"R2.fastq.gz",
            "--qualified_quality_phred","20","--length_required","50",
            "--detect_adapter_for_pe","--thread",min(threads,16),"--json",out/"fastp.json","--html",out/"fastp.html")
    elif stage=="coassemble":
        ordered=sorted((r for r in rows if r["group"]==key),key=lambda r:int(r["assembly_order"]))
        if not ordered: raise ValueError(f"No libraries for {key}")
        if a.provider is None: raise ValueError("--provider is required for the historical assemblies.")
        pairs=[[a.provider/f"{r['sample_id']}.R{m}.clean.fastq.gz" for m in (1,2)] for r in ordered]
        require(p for pair in pairs for p in pair)
        run("megahit","-1",",".join(str(x[0]) for x in pairs),"-2",",".join(str(x[1]) for x in pairs),
            "--k-list","21,29,39,59,79,99,119,141","--min-count","2","--min-contig-len","500",
            "-m","0.85","-t",threads,"-o",out/"assembly")
    elif stage=="bin":
        contigs=w/"coassemble"/key/"assembly"/"final.contigs.fa"
        require([contigs])
        for name in ["bam","metabat2","maxbin2","das"]: (out/name).mkdir()
        idx=out/"index"
        run("bowtie2-build","--threads",threads,contigs,idx)
        bams=[]
        for r in sorted((r for r in rows if r["group"]==key),key=lambda r:int(r["assembly_order"])):
            sample=r["sample_id"]; r1,r2=fastqs(w,"trim",sample); bam=out/"bam"/(sample+".bam")
            with subprocess.Popen(["bowtie2","-x",str(idx),"-1",str(r1),"-2",str(r2),"-p",str(max(1,threads//2)),"--no-unal"],stdout=subprocess.PIPE) as align:
                try: run("samtools","sort","-@",max(1,threads//2),"-o",bam,"-",stdin=align.stdout)
                finally: align.stdout.close()
                if align.wait(): raise RuntimeError("Bowtie2 failed for "+sample)
            run("samtools","index",bam);bams.append(bam)
        run("jgi_summarize_bam_contig_depths","--outputDepth",out/"depth.tsv",*bams)
        run("metabat2","-i",contigs,"-a",out/"depth.tsv","-m","1500","-t",threads,"-o",out/"metabat2/bin")
        run("jgi_summarize_bam_contig_depths","--outputDepth",out/"maxbin_depth.tsv","--noIntraDepthVariance",*bams)
        with (out/"maxbin_depth.tsv").open() as src,(out/"mean.tsv").open("w") as dest:
            depth=csv.reader(src,delimiter="\t");next(depth)
            for row in depth: print(row[0],row[2],sep="\t",file=dest)
        run("run_MaxBin.pl","-contig",contigs,"-abund",out/"mean.tsv","-out",out/"maxbin2/bin","-thread",threads)
        files=[out/"metabat2.tsv",out/"maxbin2.tsv"]
        fasta_map(out/"metabat2","fa",files[0]);fasta_map(out/"maxbin2","fasta",files[1])
        run("DAS_Tool","-i",",".join(map(str,files)),"-l","metabat2,maxbin2","-c",contigs,
            "-o",out/"das/dastool","-t",threads,"--score_threshold","0.5","--write_bins","--write_bin_evals")
    elif stage=="quality":
        bins=out/"bins";bins.mkdir()
        for group in GROUPS:
            found=sorted((w/"bin"/group/"das/dastool_DASTool_bins").glob("*.fa"))
            if not found:raise FileNotFoundError("No DAS Tool genomes for "+group)
            for src in found:os.link(src,bins/f"{group}__{src.name}")
        run("checkm2","predict","--input",bins,"--output-directory",out/"checkm2","--extension","fa",
            "--threads",threads,"--database_path",database("CHECKM2_DB"))
    elif stage=="dereplicate":
        q=w/"quality/project"
        with (q/"checkm2/quality_report.tsv").open() as handle,(out/"genome_info.csv").open("w",newline="") as dest:
            writer=csv.writer(dest);writer.writerow(["genome","completeness","contamination"]);files=[]
            for r in csv.DictReader(handle,delimiter="\t"):
                c,x=float(r["Completeness"]),float(r["Contamination"])
                if c>=50 and x<10:
                    path=q/"bins"/(r["Name"]+".fa");require([path])
                    files.append(path);writer.writerow([path.name,c,x])
        if not files: raise ValueError("No quality-eligible bins.")
        run("dRep","dereplicate",out/"run","--genomes",*sorted(files),"--genomeInfo",out/"genome_info.csv",
            "--ignoreGenomeQuality","--processors",threads,"-pa","0.9","-sa","0.95","-nc","0.30")
    elif stage=="gtdb":
        database("GTDBTK_DATA_PATH")
        run("gtdbtk","classify_wf","--genome_dir",genomes,"--extension","fa","--out_dir",out/"r226","--cpus",threads,"--skip_ani_screen")
    elif stage=="proteins":
        refs=sorted(genomes.glob("*.fa"))
        if not refs:raise FileNotFoundError(genomes)
        with (out/"proteins.faa").open("w") as merged,(out/"protein_to_mag.tsv").open("w") as ids:
            print("query","mag_id",sep="\t",file=ids)
            for fa in refs:
                faa=out/(fa.stem+".faa")
                with tempfile.TemporaryDirectory(dir=out) as temp:
                    prefixed=Path(temp)/"genome.fa"
                    with fa.open() as original,prefixed.open("w") as dst:
                        for line in original:
                            if line.startswith(">"): print(">"+fa.stem+"__"+line[1:].split()[0],file=dst)
                            else: dst.write(line)
                    run("prodigal","-i",prefixed,"-a",faa,"-d",out/(fa.stem+".fna"),"-o",out/(fa.stem+".gff"),"-f","gff","-p","meta","-q")
                with faa.open() as f:
                    for line in f:
                        if line.startswith(">"):
                            query=line[1:].split()[0]
                            print(">"+query,file=merged);print(query,fa.stem,sep="\t",file=ids)
                        else:merged.write(line)
    elif stage=="annotate":
        run("emapper.py","-i",w/"proteins/project/proteins.faa","--itype","proteins","--data_dir",database("EGGNOG_DATA_DIR"),
            "-m","diamond","--cpu",threads,"--target_orthologs","all","--tax_scope","auto","--output","catalogue","--output_dir",out)
    elif stage=="recruit":
        r1,r2=fastqs(w,"trim",key);refs=sorted(genomes.glob("*.fa"))
        if not refs:raise FileNotFoundError(genomes)
        run("coverm","genome","--coupled",r1,r2,"--genome-fasta-files",*refs,"--methods","mean","covered_fraction","relative_abundance",
            "--min-covered-fraction","0","--min-read-percent-identity","0.95","--min-read-aligned-percent","0.75",
            "--threads",threads,"--output-file",out/"coverm.tsv")
    elif stage=="classify":
        r1,r2=fastqs(w,"trim",key)
        for rank,db in [("pluspf",database("PLUSPF_DB")),("regional",database("REGIONAL_DB"))]:
            args=["--gzip-compressed","--unclassified-out",str(out/"residual#.fq")] if rank=="pluspf" else []
            mates=[r1,r2] if rank=="pluspf" else [out/"residual_1.fq",out/"residual_2.fq"]
            report=out/(rank+".report")
            run("kraken2","--db",db,"--threads",threads,"--paired","--confidence","0.10",
                "--report",report,"--output",os.devnull,*args,*mates)
            for level in ["S","G"]:
                run("bracken","-d",db,"-i",report,"-o",out/f"{rank}.{level}.tsv","-r","150","-l",level,"-t","10")
    elif stage=="humann":
        r1,r2=fastqs(w,"trim",key)
        with tempfile.TemporaryDirectory(dir=out) as temp:
            merged=Path(temp)/"paired.fastq.gz"
            with merged.open("wb") as dest:
                for mate in (r1,r2):
                    with mate.open("rb") as src:shutil.copyfileobj(src,dest,1024*1024)
            run("humann","--input",merged,"--output",out/"profile","--threads",threads,
                "--nucleotide-database",database("CHOCO_DB"),"--protein-database",database("UNIREF_DB"),
                "--metaphlan-options",f"--db_dir {database('METAPHLAN_DB')} -x mpa_vJun23_CHOCOPhlAnSGB_202307",
                "--output-basename",key,"--remove-temp-output")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=SAMPLE_STEPS+GROUP_STEPS+PROJECT_STEPS)
    p.add_argument("unit",help="Library ID, assembly group, or 'project' for a joint operation")
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--samples",type=Path,required=True)
    p.add_argument("--provider",type=Path,help="Archived assembly-input FASTQs")
    p.add_argument("--threads",type=int,default=8)
    a=p.parse_args();a.work=a.work.resolve()
    rows=read_samples(a.samples)
    valid={r["sample_id"] for r in rows} if a.stage in SAMPLE_STEPS else GROUPS if a.stage in GROUP_STEPS else ("project",)
    if a.unit not in valid or a.threads<1:p.error("Invalid stage unit or thread count")
    target=a.work/a.stage/a.unit
    if (target/"COMPLETE").is_file():print("Completed:",target);return
    target.parent.mkdir(parents=True,exist_ok=True)
    lock=target.with_name(target.name+".lock")
    lock.mkdir()
    try:
        if target.exists():raise FileExistsError("Inspect existing incomplete output: "+str(target))
        # Stable paths are important for tools that record input/output filenames.
        target.mkdir()
        run_step(a,rows,target)
        (target/"COMPLETE").write_text("ok\n")
    finally:lock.rmdir()


if __name__=="__main__": main()
