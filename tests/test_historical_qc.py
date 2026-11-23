"""Regression checks for the two study QC recipes and historical assembly inputs."""
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

WORKFLOW = Path(__file__).resolve().parents[1] / "workflow"
sys.path.insert(0, str(WORKFLOW))
import primary
from manifest import assembly_libraries
from stage import SUBDIR


class HistoricalQCTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        self.work = self.root / "work"
        self.manifest = self.root / "samples.tsv"
        self.header = "sample_id\trun\tgroup\tassembly_order\tpaired\thistorical_input_reads\n"
        self.manifest.write_text(self.header + "B\tSRR2\tZBY_water\t2\t0\t8\nA\tSRR1\tZBY_water\t1\t1\t8\n")

    def step(self, stage, unit, *options):
        args = ["primary.py", stage, unit, "--work", str(self.work), "--samples", str(self.manifest), *map(str, options)]
        with patch.object(sys, "argv", args):
            primary.main()

    def test_log_counts_have_matching_mates_and_order(self):
        options = {"pe1": ["/reads/A.R1.clean.fq"], "pe2": ["/reads/A.R2.clean.fq"]}
        log = "Lib 0 (/reads/A.R1.clean.fq,/reads/A.R2.clean.fq): pe, 8 reads, 151 max length"
        self.assertEqual(list(assembly_libraries("group", options, log)), [("A", 1, 8)])
        for bad in (log.replace("A.R2", "B.R2"), log.replace("Lib 0", "Lib 1"), "", log + "\n" + log):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                list(assembly_libraries("group", options, bad))

    def test_qc_sources_flags_versions_and_count_checks(self):
        fetched = self.work / "fetch/A"
        fetched.mkdir(parents=True)
        for name in ("R1.fastq.gz", "R2.fastq.gz", "COMPLETE"):
            (fetched / name).write_bytes(b"input")
        raw = self.root / "raw"
        raw.mkdir()
        for mate in (1, 2):
            (raw / f"A.R{mate}.raw.fastq.gz").write_bytes(b"original")

        def write_products(*cmd):
            for key in ("-o", "-O", "--html"):
                cmd[cmd.index(key) + 1].write_bytes(b"output")
            report = {"summary": {"after_filtering": {"total_reads": 8}}}
            cmd[cmd.index("--json") + 1].write_text(json.dumps(report))

        options = ["--historical-fastp", "fastp-old", "--historical-raw", raw, "--threads", 32]
        with patch("primary.subprocess.check_output", return_value="fastp v0.22.0"), patch("primary.run", side_effect=write_products) as run:
            self.step("historical-trim", "A", *options)
            cmd = run.call_args.args
            self.assertEqual(cmd[:5], ("fastp-old", "-i", raw/"A.R1.raw.fastq.gz", "-I", raw/"A.R2.raw.fastq.gz"))
            self.assertEqual(cmd[cmd.index("--detect_adapter_for_pe"):cmd.index("--json")],
                             ("--detect_adapter_for_pe", "--qualified_quality_phred", "20", "--length_required", "50", "--correction", "--thread", 8))
            self.step("historical-trim", "A", *options)
            run.assert_called_once()
            with patch("sys.stderr", new_callable=io.StringIO) as error, self.assertRaises(SystemExit):
                self.step("historical-trim", "A")  # cached execution had a different raw source
            self.assertIn("different raw input path", error.getvalue())
        with patch("primary.subprocess.check_output", return_value="fastp 1.3.2"), patch("primary.run", side_effect=write_products) as run:
            self.step("trim", "A")
            cmd = run.call_args.args
            self.assertEqual(cmd[:5], ("fastp", "-i", fetched/"R1.fastq.gz", "-I", fetched/"R2.fastq.gz"))
            self.assertNotIn("--correction", cmd)
            self.assertEqual(cmd[cmd.index("-o") + 1], self.work/"trim/A/R1.fastq.gz")
        self.manifest.write_text(self.manifest.read_text().replace("\t1\t8\n", "\t1\t10\n"))
        with self.assertRaisesRegex(ValueError, "fastp retained 8.*records 10"):
            self.step("historical-trim", "A", *options)

    def test_failed_count_check_prevents_use_of_partial_reads(self):
        self.manifest.write_text(self.header + "A\tSRR1\tZBY_water\t1\t1\t8\n")
        fetched = self.work / "fetch/A"
        fetched.mkdir(parents=True)
        for name in ("R1.fastq.gz", "R2.fastq.gz", "COMPLETE"):
            (fetched / name).write_bytes(b"raw")

        def mismatched(program, version, pair, out, threads, correction):
            self.assertEqual(pair, [fetched / f"R{m}.fastq.gz" for m in (1, 2)])
            for mate in (1, 2):
                (out / f"R{mate}.fastq.gz").write_bytes(b"filtered")
            (out/"fastp.json").write_text('{"summary":{"after_filtering":{"total_reads":6}}}')

        with patch("primary.trim_reads", side_effect=mismatched), self.assertRaisesRegex(ValueError, "retained 6.*records 8"):
            self.step("historical-trim", "A")
        self.assertFalse((self.work/"historical-trim/A/COMPLETE").exists())
        with patch("primary.run") as run, self.assertRaises(FileNotFoundError):
            self.step("coassemble", "ZBY_water")
        run.assert_not_called()

    def test_wrong_version_stops_before_filtering(self):
        with patch("primary.run") as run:
            for wanted, found in (("0.22.0", "1.3.2"), ("1.3.2", "0.22.0")):
                with self.subTest(wanted=wanted), patch("primary.subprocess.check_output", return_value=f"fastp {found}"), self.assertRaisesRegex(ValueError, "Expected fastp"):
                    primary.trim_reads("fastp", wanted, [], self.work, 8, wanted == "0.22.0")
            run.assert_not_called()

    def test_coassembly_selects_completed_pairs_in_archived_order(self):
        archived = self.root / "archived"
        archived.mkdir()
        for sample in ("A", "B"):
            folder = self.work / "historical-trim" / sample
            folder.mkdir(parents=True)
            for m in (1, 2):
                (folder / f"R{m}.fastq.gz").write_bytes(b"reconstructed")
                (archived / f"{sample}.R{m}.clean.fastq.gz").write_bytes(b"archived")
            (folder/"fastp.json").write_text('{"summary":{"after_filtering":{"total_reads":8}}}')
            (folder/"COMPLETE").write_text("ok")

        def megahit(*cmd):
            dest = cmd[cmd.index("-o") + 1]
            dest.mkdir()
            (dest / "options.json").write_text(json.dumps({"pe1": cmd[2].split(","), "pe2": cmd[4].split(",")}))

        with patch("primary.run", side_effect=megahit) as run:
            self.step("coassemble", "ZBY_water")
            self.assertEqual(run.call_args.args[2].split(","), [str(self.work/f"historical-trim/{s}/R1.fastq.gz") for s in "AB"])
            self.assertEqual(run.call_args.args[4].split(","), [str(self.work/f"historical-trim/{s}/R2.fastq.gz") for s in "AB"])
            self.step("coassemble", "ZBY_water")
            with patch("sys.stderr", new_callable=io.StringIO) as error, self.assertRaises(SystemExit):
                self.step("coassemble", "ZBY_water", "--assembly-reads", archived)
            self.assertIn("different input paths or library order", error.getvalue())
            self.work = self.root / "archive-run"
            self.step("coassemble", "ZBY_water", "--provider", archived)
            self.assertEqual(run.call_args.args[2].split(","), [str(archived/f"{s}.R1.clean.fastq.gz") for s in "AB"])
            self.assertEqual(run.call_args.args[4].split(","), [str(archived/f"{s}.R2.clean.fastq.gz") for s in "AB"])
            self.assertEqual(run.call_count, 2)

    def test_historical_plans_keep_archive_bypass_and_old_manifests(self):
        meta = self.root / SUBDIR
        meta.mkdir(parents=True)
        (meta / "workspace.json").write_text('{"mode":"fresh"}')
        lines = [self.header.rsplit("\t", 1)[0] + "\n"]
        n = 0
        for group, size in zip(primary.GROUPS, (21, 4, 5, 5)):
            for order in range(1, size + 1):
                n += 1
                lines.append(f"L{n}\tSRR{n}\t{group}\t{order}\t{int(n <= 14)}\n")
        self.manifest.write_text("".join(lines))
        cmd = [sys.executable, str(WORKFLOW/"run.py"), "historical", "--project", str(self.root),
               "--work", str(self.work), "--samples", str(self.manifest)]
        plan = subprocess.check_output(cmd + ["--historical-fastp", "fastp-old"], text=True).splitlines()
        self.assertEqual(sum(" historical-trim " in line and "--historical-fastp fastp-old" in line for line in plan), 35)
        self.assertTrue(all(" historical-trim " in line for line in plan[:35]))
        bypass = subprocess.check_output(cmd + ["--provider", str(self.root/"archive")], text=True).splitlines()
        self.assertFalse(any(" historical-trim " in line for line in bypass))
        self.assertEqual(sum(" coassemble " in line and "--assembly-reads" in line for line in bypass), 4)


if __name__ == "__main__":
    unittest.main()
