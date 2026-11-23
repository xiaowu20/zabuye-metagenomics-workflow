#!/usr/bin/env python3
"""Compose the reordered five-panel Figure 2 from checked numerical source tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Ellipse, Patch, Rectangle
import numpy as np
import pandas as pd


ZSCORE_CLIP = 2.5
ABUNDANCE_SCALE = 0.05
MIN_POSITIVE_AREA = 7.0
MAX_AREA_INCREMENT = 68.0
FIG_WIDTH_IN = 8.35
FIG_HEIGHT_IN = 10.35
FIG_DPI = 600

COL = {
    "north": "#3F6B99", "south": "#E5A442",
    "spring": "#55A47B", "summer": "#D8A53A",
    "autumn": "#7A63A8", "winter": "#7897B4",
    "Archaea": "#7561A8", "Bacteria": "#2A9D8F",
    "teal": "#34988F", "orange": "#E5A442",
    "ink": "#25282B", "muted": "#68747D",
    "grid": "#E6EBED", "edge": "#626A70",
}


def abundance_area(values: np.ndarray, maximum: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.zeros(values.shape, dtype=float)
    positive = values > 0
    denominator = np.log1p(maximum / ABUNDANCE_SCALE)
    out[positive] = MIN_POSITIVE_AREA + MAX_AREA_INCREMENT * (
        np.log1p(values[positive] / ABUNDANCE_SCALE) / denominator
    )
    return out


def add_covariance_ellipse(ax, x: np.ndarray, y: np.ndarray, color: str) -> None:
    """Add an 80% bivariate-normal covariance ellipse as a visual summary."""
    points = np.column_stack([x, y])
    covariance = np.cov(points, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    chi_square_80_df2 = -2.0 * np.log(1.0 - 0.80)
    width, height = 2.0 * np.sqrt(eigenvalues * chi_square_80_df2)
    ax.add_patch(Ellipse(
        xy=points.mean(axis=0), width=width, height=height, angle=angle,
        facecolor="none", edgecolor=color, linewidth=1.25,
        linestyle=(0, (4.0, 3.0)), zorder=1,
    ))


def style_numeric_axis(ax) -> None:
    ax.grid(True, color=COL["grid"], linewidth=0.65, zorder=0)
    for spine in ax.spines.values():
        spine.set_color("#34383C")
        spine.set_linewidth(0.72)
    ax.tick_params(length=3.0, width=0.65, color="#454C51")


def save_all(fig: mpl.figure.Figure, output_dir: Path) -> None:
    stem = output_dir / "Figure_2_top40_occurrence_reordered"
    fig.savefig(stem.with_suffix(".png"), dpi=FIG_DPI, facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"), dpi=FIG_DPI, facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    output_dir = args.output
    data_out = output_dir / 'Data'
    data_out.mkdir(parents=True, exist_ok=True)
    a_source = args.inputs / "Figure_2a_top40_domain_abundance_source.tsv"
    a_order = args.inputs / "Figure_2a_top40_domain_order.tsv"
    pcoa_source = args.inputs / "Figure_2b_n21_species_PCoA.tsv"
    permanova_source = args.inputs / "Figure_2c_n14_PERMANOVA_effects.tsv"
    genera_source = args.inputs / "Figure_2d_n14_robust_genera.tsv"
    pooled_overlap_source = args.inputs / "Figure_S18bc_overlap_thresholds.tsv"
    domain_overlap_source = args.inputs / "Figure_S19_domain_overlap_thresholds.tsv"
    required = [
        a_source, a_order, pcoa_source, permanova_source, genera_source,
        pooled_overlap_source, domain_overlap_source,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing source tables:\n" + "\n".join(missing))

    heat = pd.read_csv(a_source, sep="\t")
    genus_info = pd.read_csv(a_order, sep="\t").sort_values("display_order")
    pcoa = pd.read_csv(pcoa_source, sep="\t")
    permanova = pd.read_csv(permanova_source, sep="\t")
    genera = pd.read_csv(genera_source, sep="\t")
    pooled_overlap_all = pd.read_csv(pooled_overlap_source, sep="\t")
    domain_overlap_all = pd.read_csv(domain_overlap_source, sep="\t")
    pooled_overlap = pooled_overlap_all.loc[
        pooled_overlap_all["minimum_occurrence"] == 2
    ].copy()
    domain_overlap = domain_overlap_all.loc[
        (domain_overlap_all["minimum_occurrence"] == 2)
        & domain_overlap_all["domain"].isin(["Archaea", "Bacteria"])
    ].copy()
    domain_overlap["domain"] = pd.Categorical(
        domain_overlap["domain"], categories=["Bacteria", "Archaea"], ordered=True
    )
    domain_overlap = domain_overlap.sort_values("domain").reset_index(drop=True)
    occurrence = pd.concat([
        pd.DataFrame({
            "group": ["All features"],
            "retained_features": pooled_overlap["retained_features"].astype(int).to_numpy(),
            "spring_only": pooled_overlap["spring_only"].astype(int).to_numpy(),
            "spring_only_percent": pooled_overlap["spring_only_percent"].to_numpy(),
            "shared": pooled_overlap["shared"].astype(int).to_numpy(),
            "shared_percent": pooled_overlap["shared_percent"].to_numpy(),
            "autumn_only": pooled_overlap["autumn_only"].astype(int).to_numpy(),
            "autumn_only_percent": pooled_overlap["autumn_only_percent"].to_numpy(),
        }),
        domain_overlap.rename(columns={"domain": "group"})[
            ["group", "retained_features", "spring_only", "spring_only_percent",
             "shared", "shared_percent", "autumn_only", "autumn_only_percent"]
        ],
    ], ignore_index=True)

    if len(heat) != 840 or heat["genus"].nunique() != 40 or heat["sample_id"].nunique() != 21:
        raise ValueError("Panel a must contain 40 genera x 21 samples")
    if len(occurrence) != 3 or occurrence["group"].tolist() != ["All features", "Bacteria", "Archaea"]:
        raise ValueError("Panel b must contain All features, Bacteria and Archaea")
    if len(pcoa) != 21 or set(pcoa["basin"]) != {"north", "south"}:
        raise ValueError("Panel c must contain 21 north/south samples")
    if len(permanova) != 2 or set(permanova["term"]) != {"season", "basin"}:
        raise ValueError("Panel d must contain season and basin effects")
    if len(genera) != 12:
        raise ValueError("Panel e must contain 12 displayed genera")

    # Archive the exact input tables used by this composite.
    heat.to_csv(data_out / "Figure_2a_source.tsv", sep="\t", index=False)
    genus_info.to_csv(data_out / "Figure_2a_domain_order.tsv", sep="\t", index=False)
    occurrence.to_csv(data_out / "Figure_2b_source.tsv", sep="\t", index=False)
    pcoa.to_csv(data_out / "Figure_2c_source.tsv", sep="\t", index=False)
    permanova.to_csv(data_out / "Figure_2d_source.tsv", sep="\t", index=False)
    genera.to_csv(data_out / "Figure_2e_source.tsv", sep="\t", index=False)

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7.5,
        "axes.titlesize": 9.2, "axes.titleweight": "bold",
        "axes.labelsize": 8.0, "xtick.labelsize": 7.0, "ytick.labelsize": 7.1,
        "legend.fontsize": 6.7, "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.linewidth": 0.7,
    })

    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), facecolor="white")
    main_grid = fig.add_gridspec(
        1, 2, width_ratios=[1.30, 1.0],
        left=0.025, right=0.985, top=0.972, bottom=0.045,
        wspace=0.255,
    )

    # ------------------------------------------------------------------ panel a
    left = main_grid[0, 0].subgridspec(
        3, 3, height_ratios=[1.12, 18.8, 2.30],
        width_ratios=[0.17, 1.28, 3.35], hspace=0.025, wspace=0.035,
    )
    ax_ann = fig.add_subplot(left[0, 2])
    ax_group = fig.add_subplot(left[1, 0])
    ax_labels = fig.add_subplot(left[1, 1], sharey=ax_group)
    ax_heat = fig.add_subplot(left[1, 2], sharey=ax_group)
    ax_a_legend = fig.add_subplot(left[2, :])

    sample_info = (
        heat[["sample_id", "sample_order", "season", "basin", "site_corrected"]]
        .drop_duplicates("sample_id").sort_values("sample_order").reset_index(drop=True)
    )
    sample_order = sample_info["sample_id"].tolist()
    genus_order = genus_info["genus"].tolist()
    sample_to_x = {sample: i for i, sample in enumerate(sample_order)}
    genus_to_y = {genus: i for i, genus in enumerate(genus_order)}
    heat["x"] = heat["sample_id"].map(sample_to_x)
    heat["y"] = heat["genus"].map(genus_to_y)
    max_abundance = float(heat["relative_abundance_pct"].max())
    heat["area_composite_pt2"] = abundance_area(heat["relative_abundance_pct"], max_abundance)

    ax_ann.set_xlim(-0.5, 20.5)
    ax_ann.set_ylim(-0.5, 1.5)
    for x, row in sample_info.iterrows():
        ax_ann.add_patch(Rectangle((x - 0.48, 0.52), 0.96, 0.96,
                                   facecolor=COL[row["basin"]], edgecolor="white", linewidth=0.4))
        ax_ann.add_patch(Rectangle((x - 0.48, -0.48), 0.96, 0.96,
                                   facecolor=COL[row["season"]], edgecolor="white", linewidth=0.4))
    ax_ann.text(-0.88, 1.0, "Basin", ha="right", va="center", fontsize=7.2, fontweight="bold")
    ax_ann.text(-0.88, 0.0, "Season", ha="right", va="center", fontsize=7.2, fontweight="bold")
    ax_ann.text(-0.5, 1.87, "Dominant genera", ha="left", va="bottom",
                fontsize=9.5, fontweight="bold", color=COL["ink"], clip_on=False)
    captured = genus_info.mean_relative_abundance_pct.sum()
    ax_ann.text(-0.5, 1.56, f"Top {len(genus_info)} across {heat.sample_id.nunique()} ZBY water metagenomes ({captured:.2f}% of abundance)",
                ha="left", va="bottom", fontsize=6.5, color=COL["muted"], clip_on=False)
    ax_ann.axis("off")

    ax_group.set_xlim(0, 1)
    ax_group.set_ylim(39.5, -0.5)
    for domain in ["Archaea", "Bacteria"]:
        rows = genus_info.index[genus_info["domain"] == domain].to_numpy()
        y0 = float(rows.min()) - 0.48
        height = float(rows.max() - rows.min() + 1) - 0.04
        ax_group.add_patch(Rectangle((0.50, y0), 0.40, height,
                                     facecolor=COL[domain], edgecolor="none", alpha=0.92))
        ax_group.text(0.27, (rows.min() + rows.max()) / 2, domain, rotation=90,
                      ha="center", va="center", color=COL[domain],
                      fontsize=7.0, fontweight="bold")
    ax_group.axis("off")

    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(39.5, -0.5)
    for y, row in genus_info.iterrows():
        ax_labels.text(0.99, y, row["genus"], ha="right", va="center",
                       fontsize=6.75, fontstyle="italic", color=COL[row["domain"]])
    ax_labels.axis("off")

    cmap = LinearSegmentedColormap.from_list(
        "coral_white_teal", ["#E98770", "#FAFAF8", "#4FA99A"], N=256
    )
    norm = TwoSlopeNorm(vmin=-ZSCORE_CLIP, vcenter=0, vmax=ZSCORE_CLIP)
    positive = heat["relative_abundance_pct"] > 0
    ax_heat.scatter(
        heat.loc[positive, "x"], heat.loc[positive, "y"],
        s=heat.loc[positive, "area_composite_pt2"],
        c=heat.loc[positive, "displayed_z"], cmap=cmap, norm=norm,
        edgecolors=COL["edge"], linewidths=0.34, zorder=3,
    )
    ax_heat.set_xlim(-0.5, 20.5)
    ax_heat.set_ylim(39.5, -0.5)
    ax_heat.set_xticks(np.arange(21))
    ax_heat.set_xticklabels(sample_order, rotation=58, ha="right", rotation_mode="anchor", fontsize=5.9)
    ax_heat.set_yticks([])
    ax_heat.set_xticks(np.arange(-0.5, 21, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, 40, 1), minor=True)
    ax_heat.grid(which="minor", color=COL["grid"], linewidth=0.33)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    ax_heat.tick_params(axis="x", length=2.4, width=0.55, pad=1.0)
    for spine in ax_heat.spines.values():
        spine.set_color("#34383C")
        spine.set_linewidth(0.65)
    n_archaea = int((genus_info["domain"] == "Archaea").sum())
    ax_heat.axhline(n_archaea - 0.5, color="#7B858C", linewidth=0.9, zorder=4)
    ax_labels.axhline(n_archaea - 0.5, color="#7B858C", linewidth=0.9)

    ax_a_legend.axis("off")
    cax = ax_a_legend.inset_axes([0.015, 0.56, 0.34, 0.18])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation="horizontal", ticks=[-2, -1, 0, 1, 2])
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=5.8, length=2.0, width=0.5)
    cax.set_title("Within-genus z score", fontsize=6.5, fontweight="bold", loc="left", pad=2)

    abundance_values = [0.01, 0.1, 1, 5, 20]
    abundance_areas = abundance_area(np.asarray(abundance_values), max_abundance)
    abundance_handles = [
        ax_a_legend.scatter([], [], s=area, facecolor="white", edgecolor=COL["edge"], linewidth=0.45)
        for area in abundance_areas
    ]
    abundance_legend = ax_a_legend.legend(
        abundance_handles, [f"{value:g}" for value in abundance_values],
        title="Relative abundance (%) · log-scaled area",
        loc="upper left", bbox_to_anchor=(0.40, 0.80), frameon=False,
        ncol=5, handletextpad=0.22, columnspacing=0.55, borderaxespad=0,
        fontsize=5.8,
    )
    abundance_legend.get_title().set_fontsize(6.4)
    abundance_legend.get_title().set_fontweight("bold")
    ax_a_legend.add_artist(abundance_legend)
    annotation_handles = [
        Patch(facecolor=COL["north"], label="North"), Patch(facecolor=COL["south"], label="South"),
        Patch(facecolor=COL["spring"], label="Spring"), Patch(facecolor=COL["summer"], label="Summer"),
        Patch(facecolor=COL["autumn"], label="Autumn"), Patch(facecolor=COL["winter"], label="Winter"),
    ]
    ax_a_legend.legend(annotation_handles, [h.get_label() for h in annotation_handles],
                       loc="lower left", bbox_to_anchor=(0.40, 0.04), frameon=False,
                       ncol=6, handlelength=0.8, columnspacing=0.55, borderaxespad=0,
                       fontsize=5.7)
    ax_a_legend.text(0.365, 0.16, "Basin / season", transform=ax_a_legend.transAxes,
                     ha="right", va="center", fontsize=6.2, fontweight="bold")

    ax_ann.text(-3.65, 2.12, "a", ha="left", va="top", fontsize=13, fontweight="bold", clip_on=False)

    # ------------------------------------------------------------- panels b-e
    right = main_grid[0, 1].subgridspec(
        4, 1, height_ratios=[0.56, 1.08, 0.66, 0.91], hspace=0.42
    )
    b_grid = right[0].subgridspec(2, 1, height_ratios=[0.31, 1.0], hspace=0.03)
    ax_b_head = fig.add_subplot(b_grid[0])
    ax_b = fig.add_subplot(b_grid[1])
    c_grid = right[1].subgridspec(2, 1, height_ratios=[0.22, 1.0], hspace=0.02)
    ax_c_head = fig.add_subplot(c_grid[0])
    ax_c = fig.add_subplot(c_grid[1])
    ax_d = fig.add_subplot(right[2])
    ax_e = fig.add_subplot(right[3])

    # Panel b: primary occurrence decomposition, followed by domain subsets.
    ax_b_head.axis("off")
    ax_b_head.text(0.0, 0.98, "Spring–autumn feature recurrence", ha="left", va="top",
                   fontsize=9.0, fontweight="bold", color=COL["ink"])
    ax_b_head.text(0.0, 0.43, "Occurrence ≥2 across paired libraries; n = 14",
                   ha="left", va="center", fontsize=6.1, color=COL["muted"])
    component_meta = [
        ("spring_only_percent", "Spring only", "#58A976"),
        ("shared_percent", "Shared", "#439793"),
        ("autumn_only_percent", "Autumn only", "#7C63AD"),
    ]
    component_handles = [Patch(facecolor=color, label=label) for _, label, color in component_meta]
    ax_b_head.legend(component_handles, [handle.get_label() for handle in component_handles],
                     loc="center left", bbox_to_anchor=(0.0, -0.05), ncol=3,
                     handlelength=1.0, columnspacing=0.75, fontsize=5.6)
    left_values = np.zeros(len(occurrence), dtype=float)
    y_b = np.arange(len(occurrence))
    for column, label, color in component_meta:
        values = occurrence[column].to_numpy(dtype=float)
        ax_b.barh(y_b, values, left=left_values, height=0.57,
                  color=color, edgecolor="white", linewidth=0.55, zorder=3)
        for index, (left_value, value) in enumerate(zip(left_values, values)):
            if value >= 8.0:
                ax_b.text(left_value + value / 2, index, f"{value:.1f}%",
                          ha="center", va="center", fontsize=5.7,
                          color="white", fontweight="bold", zorder=4)
            else:
                ax_b.text(left_value + value + 0.8, index, f"{value:.1f}%",
                          ha="left", va="center", fontsize=5.4,
                          color=color, fontweight="bold", zorder=4, clip_on=False)
        left_values += values
    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels([
        f"{row.group} (n={int(row.retained_features):,})"
        for row in occurrence.itertuples(index=False)
    ], fontsize=6.0)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 106)
    ax_b.set_xticks([0, 50, 100])
    ax_b.set_xlabel("Feature union (%)", fontsize=6.5, labelpad=1.2)
    style_numeric_axis(ax_b)
    ax_b.grid(axis="y", visible=False)
    ax_b_head.text(-0.17, 1.08, "b", transform=ax_b_head.transAxes,
                   ha="left", va="top", fontsize=13, fontweight="bold")

    # Panel c header and legends.
    ax_c_head.axis("off")
    ax_c_head.text(0.0, 0.98, "Species-level composition", ha="left", va="top",
                   fontsize=9.2, fontweight="bold", color=COL["ink"])
    ax_c_head.text(0.0, 0.58, "Descriptive ordination", ha="left", va="center",
                   fontsize=6.7, color=COL["muted"])
    ax_c_head.text(0.48, 0.58, "Basin", ha="left", va="center",
                   fontsize=6.1, fontweight="bold")
    for x, basin, label in [(0.66, "north", "North"), (0.84, "south", "South")]:
        ax_c_head.scatter(x, 0.58, s=18, marker="o", color=COL[basin], clip_on=False)
        ax_c_head.text(x + 0.035, 0.58, label, ha="left", va="center", fontsize=5.7)
    ax_c_head.text(0.0, 0.12, "Season", ha="left", va="center",
                   fontsize=6.1, fontweight="bold")
    season_positions = [
        (0.20, "o", "Spring"), (0.40, "^", "Summer"),
        (0.61, "s", "Autumn"), (0.82, "D", "Winter"),
    ]
    for x, marker, label in season_positions:
        ax_c_head.scatter(x, 0.12, s=17, marker=marker, color="#1F2224", clip_on=False)
        ax_c_head.text(x + 0.032, 0.12, label, ha="left", va="center", fontsize=5.45)
    ax_c_head.set_xlim(0, 1.12)
    ax_c_head.set_ylim(0, 1)

    marker_map = {"spring": "o", "summer": "^", "autumn": "s", "winter": "D"}
    for (basin, season), group in pcoa.groupby(["basin", "season"], sort=False):
        ax_c.scatter(group["PCoA1"], group["PCoA2"], s=38,
                     marker=marker_map[season], facecolor=COL[basin],
                     edgecolor="white", linewidth=0.55, zorder=3)
    for basin, group in pcoa.groupby("basin"):
        add_covariance_ellipse(ax_c, group["PCoA1"].to_numpy(), group["PCoA2"].to_numpy(), COL[basin])
    p1 = float(pcoa["PCoA1_percent"].iloc[0])
    p2 = float(pcoa["PCoA2_percent"].iloc[0])
    ax_c.set_xlabel(f"PCoA1 ({p1:.1f}%)")
    ax_c.set_ylabel(f"PCoA2 ({p2:.1f}%)")
    ax_c.set_xlim(-0.58, 0.79)
    ax_c.set_ylim(-0.61, 0.46)
    ax_c.text(0.02, 0.03, "n = 21", transform=ax_c.transAxes, color=COL["muted"], fontsize=6.3)
    style_numeric_axis(ax_c)
    ax_c_head.text(-0.17, 1.06, "c", transform=ax_c_head.transAxes,
                   ha="left", va="top", fontsize=13, fontweight="bold")

    # Panel d: paired effect sizes.
    dplot = permanova.assign(
        label=permanova["term"].map({"season": "Season\nwithin station", "basin": "Basin\npaired stations"}),
        y=permanova["term"].map({"season": 1, "basin": 0}),
        color=permanova["term"].map({"season": COL["teal"], "basin": COL["orange"]}),
    )
    for _, row in dplot.iterrows():
        lower = row["R2"] - row["R2_CI_low"]
        upper = row["R2_CI_high"] - row["R2"]
        ax_d.errorbar(row["R2"], row["y"], xerr=np.array([[lower], [upper]]),
                      fmt="o", color=row["color"], ecolor=row["color"],
                      markersize=5.5, elinewidth=1.35, capsize=3.1, capthick=1.2, zorder=3)
        p_label = f"P = {row['P_exact']:.4f}" if row["term"] == "season" else f"P = {row['P_exact']:.1f}"
        ax_d.text(row["R2_CI_high"] + 0.018, row["y"], p_label,
                  ha="left", va="center", fontsize=6.8)
    ax_d.set_yticks([1, 0])
    ax_d.set_yticklabels(["Season\nwithin station", "Basin\npaired stations"])
    ax_d.set_xlim(0, 0.69)
    ax_d.set_ylim(-0.55, 1.42)
    ax_d.set_xlabel("PERMANOVA $R^2$ (95% station-cluster CI)")
    ax_d.text(0.02, 0.96, "Exact station-unit tests; n = 14", transform=ax_d.transAxes,
              ha="left", va="top", fontsize=6.5, color=COL["muted"])
    style_numeric_axis(ax_d)
    ax_d.text(-0.17, 1.10, "d", transform=ax_d.transAxes,
              ha="left", va="top", fontsize=13, fontweight="bold")

    # Panel e: genera retained by both paired models.
    eplot = genera.sort_values("coef", ascending=False).reset_index(drop=True)
    y_positions = np.arange(len(eplot))
    for y, row in eplot.iterrows():
        lower = row["coef"] - row["conf_low"]
        upper = row["conf_high"] - row["coef"]
        ax_e.errorbar(row["coef"], y, xerr=np.array([[lower], [upper]]), fmt="o",
                      color=COL["teal"], ecolor=COL["teal"], markersize=5.2,
                      elinewidth=1.25, capsize=2.6, capthick=1.0, zorder=3)
    ax_e.set_yticks(y_positions)
    ax_e.set_yticklabels(eplot["feature"], fontstyle="italic", fontsize=5.85)
    ax_e.tick_params(axis="y", pad=1.5)
    ax_e.invert_yaxis()
    ax_e.axvline(0, color="#9AA3A9", linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax_e.set_xlim(-0.55, 13.45)
    ax_e.set_xlabel("Spring versus autumn coefficient (95% CI)")
    ax_e.set_title("Spring-positive genera", loc="left", pad=7)
    ax_e.text(0.98, 0.03, "Retained by both paired models; n = 14",
              transform=ax_e.transAxes, ha="right", va="bottom",
              fontsize=6.15, color=COL["muted"])
    style_numeric_axis(ax_e)
    ax_e.text(-0.17, 1.12, "e", transform=ax_e.transAxes,
              ha="left", va="top", fontsize=13, fontweight="bold")

    save_all(fig, output_dir)
    plt.close(fig)

if __name__ == "__main__":
    main()
