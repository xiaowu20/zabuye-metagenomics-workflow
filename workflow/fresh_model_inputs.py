"""Join fresh recruitment, response matrices and measured study metadata for paired fits."""
import argparse
from pathlib import Path
import shutil
import tempfile

import pandas as pd
from primary import read_samples


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "samples", "design", "s4", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error("--out must be a new directory")
    rows = read_samples(args.samples)
    zby = [r["sample_id"] for r in rows if r["group"] == "ZBY_water"]
    water = [r["sample_id"] for r in rows if r["group"].endswith("_water")]
    if len(zby) != 21 or len(water) != 26:
        p.error("Expected the study's 21 Zabuye and five comparison water libraries")
    code = args.work / "responses"
    def matrix(name, response):
        return pd.read_csv(code / name / (response + ".tsv"), sep="\t").set_index("sample", verify_integrity=True)
    study=pd.read_excel(args.s4,sheet_name="MAG_catalogue_414").set_index("mag_id",verify_integrity=True)
    for scope, expected in (("all",study),("ZBY_water",study.loc[study.source_assembly=="ZBY_water"])):
        if set(matrix(scope,"MAG").columns)!=set(expected.index):
            p.error(f"Genome identifiers for {scope} differ from S4. The fixed-catalogue model labels would be invalid.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out.parent) as tmp:
        out = Path(tmp) / "models"
        out.mkdir()
        for filename in ("ZBYinternal_manifest_with_site.tsv", "sample_membership_corrected.tsv",
                         "site_chemistry.xlsx", "hydrology.xlsx"):
            shutil.copyfile(args.design / filename, out / filename)
        matrix("ZBY_water", "MAG").loc[zby].T.to_csv(out / "within_zby_relative_abundance.tsv",
                                                      sep="\t", index_label="mag_id")
        for prefix, name, libraries in (("ZBYinternal","ZBY_water",zby),
                                        ("zbyonly21_cat414","all",zby),
                                        ("crosslake26_cat414","all",water)):
            for response, label in (("GMP","graded_weighted"),("SCA","strict_carrier_abundance"),
                                    ("DCP","detected_carrier_fraction")):
                matrix(name,response).loc[libraries].to_csv(out / f"{prefix}_response_{label}.tsv",
                                                             sep="\t", index_label="sample_id")
        qc = matrix("all", "QC").rename(columns={"detected_MAGs":"n_detected_mags",
                                                "prethreshold_catalogue_capture_pct":"catalogue_capture_pct"})
        for sample in water:
            coverm = pd.read_csv(args.work/"recruit"/sample/"coverm.tsv",sep="\t",index_col=0)
            abundance = [c for c in coverm if c.lower().endswith(" relative abundance (%)")]
            if len(abundance)!=1 or not coverm.index.is_unique:
                raise ValueError(f"Unrecognised CoverM header or duplicate row in {sample}")
            qc.loc[sample,"mapping_rate_pct"] = 100-float(coverm.at["unmapped",abundance[0]])
        qc["lake"] = {r["sample_id"]:r["group"].split("_")[0] for r in rows}
        for prefix, libraries in (("zbyonly21_cat414",zby),("crosslake26_cat414",water)):
            qc.loc[libraries].to_csv(out / f"{prefix}_sample_qc.tsv",sep="\t",index_label="sample_id")
        out.rename(args.out)
    print("Prepared fresh inputs for paired_models.R, incidence.R, ion_models.R and crosslake.R")


if __name__ == "__main__":
    main()
