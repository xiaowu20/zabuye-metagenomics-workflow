"""Build the Figure 3 plotting tables from the attribution and statistical outputs."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__)
for key in ("carrier","models","palette","out"):
    p.add_argument("--"+key,type=Path,required=True)
a=p.parse_args()
SRC=a.out
SRC.mkdir(parents=True,exist_ok=False)
def read(name):return pd.read_csv(a.carrier/(name+".tsv"),sep="\t")
global_path=a.models/"global_paired_models.tsv"
ko_path=a.models/"community_KO_primary_models.tsv"
env_path=a.models/"continuous_environment_PERMANOVA_n14.tsv"
P=json.loads(a.palette.read_text())
N,CORAL=P["col_north"],P["col_coral"]
top = read("display_top12_modules").set_index("feature")
modules = top.index.tolist()
global_all = pd.read_csv(global_path, sep="\t")
ko_all = pd.read_csv(ko_path, sep="\t")
env = pd.read_csv(env_path, sep="\t")
rows = []
response_specs = [
    ("dual_species_occurrence_ge2", "none", "Read-based species", N),
    ("mag_relative_abundance_269", "n_detected_mags", "MAG composition", P["col_gold"]),
    ("graded_weighted_modules_269", "n_detected_mags", "Graded potential (GMP)", CORAL),
    ("strict_carrier_abundance_269", "n_detected_mags", "Strict carriers (SCA)", CORAL),
    ("detected_carrier_fraction_269", "n_detected_mags", "Carrier prevalence", P["col_green"]),
]
for key, adjustment, label, color in response_specs:
    x = global_all[(global_all.response == key) & (global_all.term == "basin") &
                   (global_all.sample_scope == "paired14_primary") &
                   (global_all.technical_adjustment == adjustment)].drop_duplicates()
    assert len(x) == 1, (key, x.shape)
    r = x.iloc[0]
    rows.append(dict(response=label, R2=r.R2, low=r.R2_CI_low, high=r.R2_CI_high,
                     P=r.P_exact, color=color, source_response=key))
for model, label in [("M1_MAG_catalogue_detection", "Community KO\n(MAG detection)"),
                     ("M2_CoverM_mapping_rate", "Community KO\n(mapping rate)")]:
    r = ko_all[(ko_all.model == model) & (ko_all.term == "basin")].iloc[0]
    rows.append(dict(response=label, R2=r.R2, low=r.CI95_low, high=r.CI95_high,
                     P=r.P_exact, color=P["col_autumn"], source_response=model))
panel_a = pd.DataFrame(rows)
panel_a.to_csv(SRC / "Figure3a_response_effects.tsv", sep="\t", index=False)

RESP = [
    ("MAG_relative_abundance_269", "MAG composition · 269"),
    ("graded_modules_269", "GMP · 269"),
    ("strict_modules_269", "SCA · 269"),
    ("carrier_fraction_269", "Prevalence · 269"),
    ("graded_modules_414", "GMP · 414"),
    ("strict_modules_414", "SCA · 414"),
    ("carrier_fraction_414", "Prevalence · 414"),
]
# Response identifiers are checked rather than inferred from plot labels.
available = set(env.response)
for i, (key, label) in enumerate(RESP):
    if key not in available and "carrier_prevalence" in key:
        alt = key.replace("carrier_prevalence", "detected_carrier_fraction")
        if alt in available:
            RESP[i] = (alt, label)
assert set(k for k, _ in RESP).issubset(available), available
sal_rows = []
for cat in [269, 414]:
    for key, short in [(f"graded_modules_{cat}", "GMP"),
                       (f"strict_modules_{cat}", "SCA")]:
        base = env[(env.response == key) & (env.model == "baseline") & (env.term == "basin")].iloc[0]
        cond = env[(env.response == key) & (env.model == "salinity_conditioned") & (env.term == "basin")].iloc[0]
        sal_rows.append(dict(catalogue=cat, response=short, baseline_R2=base.R2,
                             conditioned_R2=cond.R2, baseline_P=base.P_restricted,
                             conditioned_P=cond.P_restricted,
                             attenuation_pct=100*(1-cond.R2/base.R2)))
sal = pd.DataFrame(sal_rows)
sal.to_csv(SRC / "Figure3d_salinity_attenuation.tsv", sep="\t", index=False)

gt = read("genus_module_attribution")
gt = gt[(gt.response == "GMP") & gt.module_id.isin(modules)]
rank = gt.groupby("genus_key").delta_contribution.apply(lambda x: np.abs(x).sum()).sort_values(ascending=False)
keys = rank.head(10).index.tolist()
lookup = gt.groupby("genus_key").label.first()
mat = gt.pivot(index="genus_key", columns="module_id", values="delta_contribution").reindex(index=keys, columns=modules)
mat.loc["Other lineages (net)"] = gt[~gt.genus_key.isin(keys)].groupby("module_id").delta_contribution.sum().reindex(modules)
lineage_labels = [lookup[k] for k in keys] + ["Other lineages (net)"]
net = top.delta_NZ_minus_SZ.to_numpy()
assert np.allclose(mat.sum(axis=0), net, atol=1e-12)
mat.to_csv(SRC / "Figure3c_genus_attribution.tsv", sep="\t")
pd.DataFrame({"lineage_key": keys + ["Other lineages (net)"],
              "display_label": lineage_labels}).to_csv(SRC / "Figure3c_lineage_labels.tsv", sep="\t", index=False)
top.to_csv(SRC/"Figure3b_module_means.tsv",sep="\t")
