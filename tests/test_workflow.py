"""Synthetic examples for sample identity, KO logic and response denominators."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

WORKFLOW = Path(__file__).resolve().parents[1] / "workflow"
sys.path.insert(0, str(WORKFLOW))
from module_rules import calculate, scores
from primary import read_samples


class WorkflowTests(unittest.TestCase):
    def test_single_gene_with_empty_required_members(self):
        members = pd.DataFrame({"gene": ["example"], "KO": ["test_ko"]})
        definitions = pd.DataFrame([dict(module_id="choice", scoring_function="single_gene_or",
                                         required_members=float("nan"), optional_members="example",
                                         validation_policy="none", output_field="score")])
        catalog = pd.DataFrame({"mag_id":["a","b"], "quality_tier":"MQ", "source_assembly":"unit",
                                "gtdb_phylum":"test", "gtdb_genus":"test"})
        ko = pd.DataFrame({"test_ko":[True,False]}, index=["a","b"])
        graded, strict = scores(calculate(ko, definitions, members, catalog, empty_member_fix=True))
        self.assertEqual(graded.choice.tolist(), [1,0])
        self.assertEqual(strict.choice.tolist(), [1,0])

    def test_manifest_and_coverm_denominators(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root/"samples.tsv"
            manifest.write_text("sample_id\trun\tgroup\tassembly_order\tpaired\nL\tSRR1\tZBY_water\t1\t0\n")
            reads = root/"recruit/L"
            reads.mkdir(parents=True)
            table = pd.DataFrame({"Genome":["a","b","c","unmapped"],"L Mean":[3,1,4,0],
                                  "L Covered Fraction":[.9,.5,.49,0],
                                  "L Relative Abundance (%)":[30,10,20,40]})
            table.to_csv(reads/"coverm.tsv", sep="\t", index=False)
            for name, values in (("scores",[1,.4,1]),("complete",[1,0,1])):
                pd.DataFrame({"id":["a","b","c"],"m":values}).to_csv(root/(name+".tsv"),sep="\t",index=False)
            out = root/"output"
            subprocess.run([sys.executable,str(WORKFLOW/"responses.py"),"--samples",str(manifest),"--work",str(root),
                            "--out",str(out),"--scores",str(root/"scores.tsv"),"--complete",str(root/"complete.tsv")],
                           check=True)
            for key, expected in (("GMP",.85),("SCA",.75),("DCP",.5)):
                self.assertAlmostEqual(pd.read_csv(out/(key+".tsv"),sep="\t").m.iloc[0],expected)
            qc = pd.read_csv(out/"QC.tsv",sep="\t").iloc[0]
            self.assertEqual(qc.detected_MAGs, 2)
            self.assertEqual(qc.prethreshold_catalogue_capture_pct, 60)
            with manifest.open("a") as f:
                f.write("M\tSRR1\tZBY_water\t2\t0\n")
            with self.assertRaises(ValueError):
                read_samples(manifest)


if __name__ == "__main__":
    unittest.main()
