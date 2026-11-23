"""Prepare figure inputs locally from the study supplements."""
import argparse
from pathlib import Path
import shutil
import tempfile
from zipfile import ZipFile

import pandas as pd
from stage import SUBDIR


def prepare(books, staged, destination, addendum):
    if destination.exists():
        raise FileExistsError(f"Choose a new figure-input directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
        inputs = Path(temp) / "inputs"
        fig2 = inputs / "figure2"
        fig3 = inputs / "figure3/source_data"
        fig2.mkdir(parents=True)
        fig3.mkdir(parents=True)
        with pd.ExcelFile(books / "Supplementary_Data_S1.xlsx") as s1:
            heat = pd.read_excel(s1, sheet_name="Fig2a_genus_abundance")
            order = ["domain", "mean_abundance_rank", "genus", "mean_relative_abundance_pct", "display_order"]
            panels = {
                "Figure_2a_top40_domain_abundance_source": heat,
                "Figure_2a_top40_domain_order": heat[order].drop_duplicates().sort_values("display_order"),
                "Figure_2b_n21_species_PCoA": pd.read_excel(s1, sheet_name="Fig2c_species_PCoA"),
                "Figure_2c_n14_PERMANOVA_effects": pd.read_excel(s1, sheet_name="Species_PERMANOVA").query(
                    "sample_scope == 'paired14_primary' and occurrence_filter == 2"),
                "Figure_S18bc_overlap_thresholds": pd.read_excel(s1, sheet_name="Overlap_thresholds"),
                "Figure_S19_domain_overlap_thresholds": pd.read_excel(s1, sheet_name="Domain_overlap"),
            }
        panels["Figure_2d_n14_robust_genera"] = pd.read_excel(
            books / "Supplementary_Data_S2.xlsx", sheet_name="Fig2e_displayed_genera")
        for name, frame in panels.items():
            frame.to_csv(fig2 / (name + ".tsv"), sep="\t", index=False)
        for source in sorted((staged / "inputs/si_display_sources").glob("*")):
            if source.suffix == ".tsv" or source.name == "palette.json":
                (inputs / "si_display").mkdir(exist_ok=True)
                shutil.copyfile(source, inputs / "si_display" / source.name)
        main = staged / "inputs/main_figure_inputs"
        for name in ("Figure3a_response_effects", "Figure3b_module_means", "Figure3c_genus_attribution",
                     "Figure3c_lineage_labels", "Figure3d_salinity_attenuation"):
            shutil.copyfile(main / (name + ".tsv"), fig3 / (name + ".tsv"))
        (inputs / "main_figure_inputs").mkdir()
        shutil.copyfile(main / "archival_midpoint_halo84.nwk", inputs / "main_figure_inputs/archival_midpoint_halo84.nwk")
        extra = inputs / "extra"
        extra.mkdir()
        names = ("environment_n14", "ion_PCA_loadings", "FigureS6c_ion_effects", "FigureS6d_basin_R2",
                 "FigureS6d_basin_P", "display_top16_carriers", "paired_manifest_14", "MAG_abundance_sample14",
                 "ion_PCA_variance", "predictor_VIF_n14", "MAG_differential")
        with ZipFile(addendum) as bundle:
            for name in names:
                filename = name + ".tsv"
                (extra / filename).write_bytes(bundle.read("source_data/" + filename))
        if not (inputs / "si_display/palette.json").is_file():
            raise ValueError("S5 display inputs are missing")
        inputs.rename(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("books", "project", "out"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--addendum", type=Path, help="Defaults to Supplementary_Data_S3_Fig3_Addendum.zip in --books")
    args = parser.parse_args()
    prepare(args.books, args.project / SUBDIR, args.out,
            args.addendum or args.books / "Supplementary_Data_S3_Fig3_Addendum.zip")
    print(f"Prepared inputs in {args.out}. Use figures.py --group paper with the same books and project.")


if __name__ == "__main__":
    main()
