"""Basin contrasts and exact additive carrier attribution; no sample pseudoreplication."""
from pathlib import Path
import hashlib
import itertools
import json
import platform
import time
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", type=Path)
parser.add_argument("out", type=Path)
args = parser.parse_args()
INP, OUT = args.inputs, args.out
OUT.mkdir(parents=True, exist_ok=False)
START = time.time()
QA = {}

def read(name, index=True):
    x = pd.read_csv(INP / name, sep="\t", index_col=0 if index else None)
    assert not x.index.duplicated().any()
    assert not x.columns.duplicated().any()
    return x

def save(x, name, index=False):
    x.to_csv(OUT / f"{name}.tsv", sep="\t", index=index, float_format="%.17g")

def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    s = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate(s[::-1])[::-1])
    return q

catalog = read("catalogue_414.tsv").fillna("")
scores = read("scores_414.tsv")
complete = read("complete_414.tsv")
zbyids = read("scores_269.tsv").index
depth, breadth, raw = [read(f"coverm_{x}.tsv") for x in [
    "mean_depth", "covered_fraction", "relative_abundance"]]
meta = read("paired_manifest_14.tsv")
roles = read("sample_roles_35.tsv")
mods = list(scores.columns)
assert scores.shape == (414, 49) and len(zbyids) == 269
assert set(catalog.index) == set(scores.index) == set(raw.index) == set(complete.index)
assert set(raw.columns) == set(depth.columns) == set(breadth.columns) == set(roles.index)
for x in [scores, complete, depth, breadth, raw]:
    assert np.isfinite(x.to_numpy()).all() and (x.to_numpy() >= 0).all()
assert scores.le(1).all().all() and breadth.le(1).all().all()
assert set(np.unique(complete)) <= {0, 1}
depth, breadth, raw, complete = [x.loc[scores.index] for x in [depth, breadth, raw, complete]]
QA.update({"finite_nonnegative_matrices": True, "unique_identifiers": True,
           "MAGs": 414, "primary_MAGs": 269, "modules": 49,
           "samples": 35, "primary_libraries": 14, "station_counts_NZ_SZ": [4, 3]})

stations = ["N1", "N2", "N3", "N4", "S1", "S2", "S3"]
group = np.array([True] * 4 + [False] * 3)
siteidx = np.array([stations.index(s) for s in meta.station])
pair_ids = meta.index.tolist()
assert meta.groupby("station").season.nunique().eq(2).all()
assignments = []
for north in itertools.combinations(range(7), 4):
    b = np.zeros(7, dtype=bool)
    b[list(north)] = True
    assignments.append(b)
assignments = np.array(assignments)
save(pd.DataFrame(assignments.astype(int), columns=stations).assign(
    assignment=np.arange(35)), "basin_assignments_35")
bootN = np.array(list(itertools.product(range(4), repeat=4)))
bootS = np.array(list(itertools.product(range(4, 7), repeat=3)))

def station_means(x):
    return x.loc[pair_ids].groupby(meta.station).mean().loc[stations]

def welch_t(a, b):
    d = a.mean(axis=0) - b.mean(axis=0)
    se2 = a.var(axis=0, ddof=1) / len(a) + b.var(axis=0, ddof=1) / len(b)
    result = np.zeros_like(d)
    ok = se2 > 1e-28
    result[ok] = d[ok] / np.sqrt(se2[ok])
    result[~ok & (np.abs(d) > 1e-14)] = np.sign(d[~ok & (np.abs(d) > 1e-14)]) * 1e100
    return result

def test_features(mat, family):
    """mat is 7 stations by features; include invariant features in BH family."""
    v = mat.to_numpy(dtype=float)
    obs = welch_t(v[group], v[~group])
    perm = np.array([welch_t(v[g], v[~g]) for g in assignments])
    p = (np.abs(perm) >= np.abs(obs)[None, :] - 1e-10).mean(axis=0)
    invariant = np.ptp(v, axis=0) < 1e-14
    p[invariant] = 1
    n = v[group].mean(axis=0)
    s = v[~group].mean(axis=0)
    bn, bs = v[bootN].mean(axis=1), v[bootS].mean(axis=1)
    boot = (bn[:, None, :] - bs[None, :, :]).reshape(-1, v.shape[1])
    lo, hi = np.quantile(boot, [.025, .975], axis=0)
    lfc = np.full(v.shape[1], np.nan)
    both = (n > 0) & (s > 0)
    lfc[both] = np.log2(n[both] / s[both])
    out = pd.DataFrame(dict(feature=mat.columns, family=family,
        NZ_mean=n, SZ_mean=s, delta_NZ_minus_SZ=n-s,
        CI95_low=lo, CI95_high=hi, log2_mean_ratio=lfc,
        zero_mean_basin=np.where(n == 0, "NZ", np.where(s == 0, "SZ", "")),
        welch_t=obs, P_exact=p, q_BH=bh(p), invariant=invariant,
        n_stations_NZ=4, n_stations_SZ=3, n_assignments=35))
    out["FDR05"] = out.q_BH < .05
    out["nominal05"] = out.P_exact < .05
    return out

defs = pd.read_csv(INP / "module_definitions.csv").fillna("").set_index("module_id")
members = pd.read_csv(INP / "KO_panel_members.csv").fillna("")
long = pd.read_csv(INP / "module_status_long.tsv", sep="\t").fillna("")
def category(m):
    return "Osmoadaptation" if m.startswith("Salt_") else (
        "Nitrogen" if m.startswith("N_") else "Sulfur" if m.startswith("S_") else "Carbon")
def module_kos(m, visited=None):
    visited = set() if visited is None else visited
    if m in visited: return set()
    visited.add(m)
    row = defs.loc[m]
    genes = set(filter(None, (str(row.required_members) + "," + str(row.optional_members)).split(",")))
    k = set(members.loc[members.gene.isin(genes), "KO"])
    for child in genes.intersection(set(defs.index)):
        k |= module_kos(child, visited)
    return k
md = pd.DataFrame({"module_id": mods, "category": [category(m) for m in mods],
                   "KO_members": [";".join(sorted(module_kos(m))) for m in mods]})
md = md.merge(defs.reset_index(), on="module_id", validate="one_to_one")
save(md, "module_dictionary")

def qualified_labels(rank):
    names = []
    for mag, row in catalog.iterrows():
        parts = [row.domain_gtdb, row.gtdb_phylum, row.gtdb_class,
                 row.gtdb_order, row.gtdb_family]
        if rank == "genus": parts.append(row.gtdb_genus)
        cleaned = [str(x) if str(x) not in ["", "nan", "Unclassified"] else "unclassified" for x in parts]
        names.append(";".join(cleaned))
    return pd.Series(names, index=catalog.index, name=rank + "_key")

for rank in ["genus", "family"]:
    catalog[rank + "_key"] = qualified_labels(rank)
    col = "gtdb_" + rank
    labels = catalog[col].astype(str)
    missing = labels.isin(["", "Unclassified", "nan"])
    ancestor_cols = (["gtdb_family"] if rank == "genus" else []) + [
        "gtdb_order", "gtdb_class", "gtdb_phylum", "domain_gtdb"]
    parent = catalog[ancestor_cols].replace(
        {"": np.nan, "Unclassified": np.nan, "nan": np.nan}).bfill(axis=1).iloc[:, 0].fillna("lineage")
    catalog[rank + "_label"] = np.where(missing, "Unclassified " + parent, labels)
save(catalog.reset_index(), "MAG_taxonomy_quality")

responses = {}
abundances = {}
qc_rows = []
reproduction = []
families = {}
sensitivity_rows = []
primary_adj = {}
primary_contrasts = {}
for catn, ids in [(269, zbyids), (414, scores.index)]:
    for threshold in [.10, .25, .50]:
        D = ((depth.loc[ids] >= 1) & (breadth.loc[ids] >= threshold)).T
        r = raw.loc[ids].T.where(D, 0.0)
        den = r.sum(axis=1)
        # No primary library has a zero denominator; non-primary zeros stay explicit.
        assert (den.loc[pair_ids] > 0).all()
        A = r.div(den.replace(0, np.nan), axis=0).fillna(0)
        qc = pd.DataFrame({"sample_id": A.index, "catalogue": catn, "breadth": threshold,
            "n_detected": D.sum(axis=1).values, "capture_before_filter_pct": raw.loc[ids].sum().values,
            "capture_after_filter_pct": den.values, "normalized_sum": A.sum(axis=1).values})
        qc_rows.append(qc)
        mats = {
            "GMP": A @ scores.loc[ids],
            "SCA": A @ complete.loc[ids],
            "DCP": (D.astype(float) @ complete.loc[ids]).div(D.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        }
        abundances[(catn, threshold)] = A
        for response, Y in mats.items():
            responses[(catn, threshold, response)] = Y
            sm = station_means(Y)
            tab = test_features(sm, response)
            tab["catalogue"] = catn
            tab["breadth"] = threshold
            sensitivity_rows.append(tab)
            if catn == 269 and threshold == .5:
                families[response] = tab.copy()
                primary_contrasts[response] = tab.set_index("feature").delta_NZ_minus_SZ
                save(Y.loc[pair_ids].rename_axis("sample_id"), f"{response}_sample14", True)
                save(sm.rename_axis("station"), f"{response}_station7", True)
                ref = read(f"reference_{response}_269_n21.tsv")
                err = float(np.max(np.abs(Y.loc[ref.index, ref.columns] - ref)))
                reproduction.append({"matrix": response+"_269_n21", "max_abs_error": err})
                assert err < 1e-12
            if catn == 414 and threshold == .5 and response == "GMP":
                ref = read("reference_GMP_414_n35.tsv")
                err = float(np.max(np.abs(Y.loc[ref.index, ref.columns] - ref)))
                reproduction.append({"matrix": "GMP_414_n35", "max_abs_error": err})
                assert err < 1e-12
save(pd.concat(qc_rows, ignore_index=True), "sample_QC_all_catalogues_thresholds")
save(pd.DataFrame(reproduction), "response_reproduction")
save(pd.concat(sensitivity_rows, ignore_index=True), "module_catalogue_breadth_sensitivity")

A = abundances[(269, .5)].loc[pair_ids]
AS = station_means(A)
save(A.rename_axis("sample_id"), "MAG_abundance_sample14", True)
save(AS.rename_axis("station"), "MAG_abundance_station7", True)
D = ((depth.loc[zbyids, pair_ids] >= 1) & (breadth.loc[zbyids, pair_ids] >= .5)).T
ndet = D.sum(axis=1).groupby(meta.station).mean().loc[stations].to_numpy()
zdet = (ndet - ndet.mean()) / ndet.std(ddof=1)
X = np.column_stack([np.ones(7), group.astype(float), zdet])
assert np.linalg.matrix_rank(X) == 3
beta_A = np.linalg.lstsq(X, AS.to_numpy(), rcond=None)[0][1]
save(pd.DataFrame({"station": stations, "basin": np.where(group, "NZ", "SZ"),
    "mean_detected_MAGs": ndet}), "station_technical_diagnostics")
QA["station_detection_basin_correlation"] = float(np.corrcoef(group.astype(float), zdet)[0, 1])

allfunc = pd.concat([families[r] for r in ["GMP", "SCA", "DCP"]], ignore_index=True)
allfunc["q_BH_joint147"] = bh(allfunc.P_exact)
allfunc = allfunc.merge(md[["module_id", "category", "validation_policy"]], left_on="feature", right_on="module_id")
save(allfunc, "functional_differential_all147")

tax_tables = []
for rank in ["MAG", "genus", "family"]:
    if rank == "MAG":
        mat = AS
    else:
        keys = catalog.loc[AS.columns, rank + "_key"]
        mat = AS.T.groupby(keys).sum().T
        save(mat.rename_axis("station"), f"{rank}_abundance_station7", True)
    tb = test_features(mat, rank)
    if rank == "MAG":
        tb = tb.merge(catalog.reset_index(), left_on="feature", right_on="mag_id", validate="one_to_one")
        tb["detected_libraries"] = D.sum().reindex(tb.feature).values
        tb["detected_stations"] = D.groupby(meta.station).any().sum().reindex(tb.feature).values
    else:
        lookup = catalog.groupby(rank+"_key")[rank+"_label"].first()
        tb["label"] = tb.feature.map(lookup)
        tb["n_MAGs"] = tb.feature.map(catalog.loc[zbyids].groupby(rank+"_key").size())
    tax_tables.append(tb)
taxonomy_test = pd.concat(tax_tables, ignore_index=True)
taxonomy_test["q_BH_joint_taxonomy"] = bh(taxonomy_test.P_exact)
save(taxonomy_test, "taxonomy_differential_all")
for rank in ["MAG", "genus", "family"]:
    save(taxonomy_test[taxonomy_test.family.eq(rank)], f"{rank}_differential")

# Feature effects with technical adjustment and direction checks.
diagnostics = []
for response in ["GMP", "SCA", "DCP"]:
    Y = responses[(269, .5, response)].loc[pair_ids]
    S = station_means(Y)
    adj = np.linalg.lstsq(X, S.to_numpy(), rcond=None)[0][1]
    primary_adj[response] = pd.Series(adj, index=mods)
    b = families[response].set_index("feature")
    spring = Y.loc[meta.season.eq("spring")]
    autumn = Y.loc[meta.season.eq("autumn")]
    sd = spring.loc[meta.loc[spring.index].basin.eq("north")].mean() - spring.loc[meta.loc[spring.index].basin.eq("south")].mean()
    ad = autumn.loc[meta.loc[autumn.index].basin.eq("north")].mean() - autumn.loc[meta.loc[autumn.index].basin.eq("south")].mean()
    loos = []
    for omitted in range(7):
        g = group.copy()
        keep = np.arange(7) != omitted
        v = S.to_numpy()
        loos.append(v[keep & g].mean(axis=0) - v[keep & ~g].mean(axis=0))
    loo = np.array(loos)
    out = pd.DataFrame({"module_id": mods, "response": response,
        "delta": b.loc[mods].delta_NZ_minus_SZ.values, "detected_MAG_adjusted_beta": adj,
        "spring_delta": sd.reindex(mods).values, "autumn_delta": ad.reindex(mods).values,
        "LOSO_min_delta": loo.min(axis=0), "LOSO_max_delta": loo.max(axis=0),
        "LOSO_same_direction_count": (np.sign(loo) == np.sign(b.loc[mods].delta_NZ_minus_SZ.values)).sum(axis=0)})
    out["same_direction_both_seasons"] = (np.sign(out.spring_delta) == np.sign(out.delta)) & (np.sign(out.autumn_delta) == np.sign(out.delta))
    out["same_direction_detection_adjusted"] = np.sign(out.detected_MAG_adjusted_beta) == np.sign(out.delta)
    diagnostics.append(out)
save(pd.concat(diagnostics), "module_direction_diagnostics")

# Exact arithmetic attribution, separate from inferential testing.
attrib = []
delta_A = AS.loc[np.array(stations)[group]].mean() - AS.loc[np.array(stations)[~group]].mean()
mean_N, mean_S = AS.iloc[:4].mean(), AS.iloc[4:].mean()
for response, C in [("GMP", scores.loc[zbyids]), ("SCA", complete.loc[zbyids])]:
    C = C.loc[AS.columns, mods]
    for j, m in enumerate(mods):
        c = C[m]
        term = delta_A * c
        adjusted = beta_A * c.to_numpy()
        target_delta = primary_contrasts[response][m]
        target_beta = primary_adj[response][m]
        assert abs(term.sum() - target_delta) < 1e-12
        assert abs(adjusted.sum() - target_beta) < 1e-12
        absden = np.abs(term).sum()
        block = pd.DataFrame({"mag_id": AS.columns, "module_id": m, "response": response,
            "score_or_complete": c.values,
            "NZ_contribution": (mean_N * c).values, "SZ_contribution": (mean_S * c).values,
            "delta_contribution": term.values, "adjusted_contribution": adjusted,
            "share_absolute_MAG_difference": np.abs(term.values)/absden if absden > 0 else 0,
            "sign_matches_module": np.sign(term.values) == np.sign(target_delta)})
        attrib.append(block)
attr = pd.concat(attrib, ignore_index=True).merge(
    catalog.reset_index(), on="mag_id", validate="many_to_one")
status_cols = ["mag_id", "module_id", "potential_status", "validated_status", "validation_policy"]
attr = attr.merge(long[status_cols], on=["mag_id","module_id"], validate="many_to_one")
save(attr, "MAG_module_attribution_all")
for rank in ["genus", "family", "gtdb_phylum"]:
    key = rank+"_key" if rank in ["genus", "family"] else rank
    cols = ["NZ_contribution","SZ_contribution","delta_contribution","adjusted_contribution"]
    agg = attr.groupby(["response", "module_id", key])[cols].sum().reset_index()
    if rank in ["genus", "family"]:
        lookup = catalog.groupby(key)[rank+"_label"].first()
        agg["label"] = agg[key].map(lookup)
    else:
        agg["label"] = agg[key]
    den = attr.groupby(["response","module_id"]).delta_contribution.apply(lambda v: np.abs(v).sum())
    agg["share_abs_MAG_denominator"] = [
        abs(r.delta_contribution) / den.loc[(r.response,r.module_id)]
        if den.loc[(r.response,r.module_id)] > 0 else 0 for r in agg.itertuples()]
    save(agg, f"{rank}_module_attribution")
QA["attribution_identity_max_error_bound"] = 1e-12
QA["attribution_interpretation"] = "algebraic decomposition, not causation or independent correlation"

# Predeclared display ranking: top absolute raw contrasts, not top significant P.
top = families["GMP"].loc[lambda x: ~x.invariant].copy()
top["abs_delta"] = top.delta_NZ_minus_SZ.abs()
top = top.sort_values(["abs_delta", "feature"], ascending=[False, True]).head(12)
save(top, "display_top12_modules")
chosen = top.feature.tolist()
importance = attr[(attr.response == "GMP") & attr.module_id.isin(chosen)].groupby("mag_id").delta_contribution.apply(lambda v: np.abs(v).sum()).sort_values(ascending=False)
topmag = importance.head(16).index
mt = taxonomy_test[taxonomy_test.family.eq("MAG")].set_index("feature").loc[topmag].copy()
mt["sum_absolute_contributions_top12"] = importance.loc[topmag]
save(mt.rename_axis("feature").reset_index(), "display_top16_carriers")
save(attr[(attr.response == "GMP") & attr.module_id.isin(chosen) & attr.mag_id.isin(topmag)], "display_MAG_module_links")

# Reproduce original aggregate Fig3 statistics, without endorsing conditional inference.
aggregate = []
for catn in [269, 414]:
    nd = ((depth.loc[abundances[(catn,.5)].columns, pair_ids] >= 1) &
          (breadth.loc[abundances[(catn,.5)].columns, pair_ids] >= .5)).sum().values
    Z = np.column_stack([np.ones(14), (meta.season == "autumn").astype(float), nd])
    H0 = Z @ np.linalg.pinv(Z)
    def stat(G, b):
        X1 = np.column_stack([Z, b])
        H1 = X1 @ np.linalg.pinv(X1)
        term = np.trace((H1-H0) @ G)
        resid = np.trace((np.eye(14)-H1) @ G)
        return term / np.trace(G), term/(resid/(14-np.linalg.matrix_rank(X1)))
    for resp in ["GMP", "SCA", "DCP"]:
        Y = responses[(catn,.5,resp)].loc[pair_ids]
        dist = squareform(pdist(Y, metric="braycurtis"))
        J = np.eye(14) - np.ones((14,14))/14
        G = -.5 * J @ (dist**2) @ J
        r2, f = stat(G, group[siteidx].astype(float))
        fs = [stat(G,b[siteidx].astype(float))[1] for b in assignments]
        aggregate.append(dict(catalogue=catn,response=resp,R2=r2,F=f,
            P_original_label_scheme=float(np.mean(np.array(fs)>=f-1.5e-8)),
            role="version-alignment reproduction; conditional exchangeability not newly certified"))
save(pd.DataFrame(aggregate), "aggregate_Fig3_reproduction")
ref = read("reference_aggregate_statistics.tsv", index=False)
for row in aggregate:
    rn = {"GMP":"graded_weighted_modules", "SCA":"strict_carrier_abundance", "DCP":"detected_carrier_fraction"}[row["response"]] + "_" + str(row["catalogue"])
    selected = ref[(ref.response==rn)&(ref.term=="basin")&(ref.scope=="paired14_primary")&(ref.model=="n_detected_mags")]
    assert len(selected)==1
    assert abs(row["R2"]-selected.R2.iloc[0])<1e-10
QA["aggregate_version_reproduction"] = True

summary = {
    "counts": {"MAG":269,"modules":49,"stations":7,"libraries":14,"permutations":35,"bootstrap_draws":6912},
    "functional": {r:{"nominal05":int(families[r].nominal05.sum()),
        "FDR05":int(families[r].FDR05.sum()), "min_q":float(families[r].q_BH.min()),
        "noninvariant":int((~families[r].invariant).sum())} for r in families},
    "taxonomy": {r:{"tested":int((taxonomy_test.family==r).sum()),
        "nominal05":int(taxonomy_test.loc[taxonomy_test.family==r,"nominal05"].sum()),
        "FDR05":int(taxonomy_test.loc[taxonomy_test.family==r,"FDR05"].sum()),
        "min_q":float(taxonomy_test.loc[taxonomy_test.family==r,"q_BH"].min())} for r in ["MAG","genus","family"]},
    "top_modules": chosen, "top_MAGs": list(topmag),
    "runtime_seconds": time.time()-START,
    "software": {"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__},
    "HPC": "Not required: 414 x 35 abundance and 414 x 49 annotation matrices analysed locally."
}
(OUT/"analysis_QA.json").write_text(json.dumps(QA,indent=2),encoding="utf-8")
(OUT/"result_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
