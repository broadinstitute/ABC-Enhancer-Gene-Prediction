"""Figures for the H3K27ac counting-mode CRISPR benchmark report.

Reads the CRISPR_comparison performance summaries for the split DHS-only and
ATAC-only comparisons and builds three Nature-style figures comparing
H3K27ac counting geometry (element vs. element+-flank vs. flank-only) and
flank width (150/200/500bp).

Regenerate:
    python analysis/h3k27ac_counting_report/make_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from nature_style import apply_rcparams, figsize, save_fig, add_panel_label

CRISPR_RESULTS = Path(
    "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/workflow/results"
)
PRED_CONFIG_DIR = Path(
    "/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/resources/pred_config"
)
PLOTS_DIR = Path(__file__).parent / "plots"

ASSAYS = {
    "DHS": {"comparison": "2026_0721_h3k27ac_counting_DHS", "label": "DNase", "only_id": "K562_DHS_only"},
    "ATAC": {"comparison": "2026_0721_h3k27ac_counting_ATAC", "label": "ATAC", "only_id": "K562_ATAC_only"},
}

# Plot order: accessibility alone, then element, then element+-flank family, then flank-only family.
VARIANT_ORDER = ["only", "element", "element_flank150", "element_flank200", "element_flank500",
                  "flankonly150", "flankonly200", "flankonly500"]
VARIANT_LABELS = {
    "only": "accessibility only", "element": "element",
    "element_flank150": "elem. +/-150bp", "element_flank200": "elem. +/-200bp", "element_flank500": "elem. +/-500bp",
    "flankonly150": "flank only 150bp", "flankonly200": "flank only 200bp", "flankonly500": "flank only 500bp",
}
DIST_BIN_ORDER = ["[0,10)", "[10,100)", "[100,2.5e+03)"]
DIST_BIN_LABELS = {"[0,10)": "0-10 kb", "[10,100)": "10-100 kb", "[100,2.5e+03)": "100 kb-2.5 Mb"}


def _variant_key(pred_uid: str, only_id: str) -> str | None:
    """Map a pred_uid (e.g. 'K562_DHS_H3K27ac_flankonly500.ABC.Score') to a VARIANT_ORDER key."""
    if pred_uid.startswith(only_id):
        return "only"
    if "H3K27ac_element_flank" in pred_uid:
        return "element_flank" + pred_uid.split("element_flank")[1].split(".")[0]
    if "H3K27ac_flankonly" in pred_uid:
        return "flankonly" + pred_uid.split("flankonly")[1].split(".")[0]
    if "H3K27ac_element" in pred_uid:
        return "element"
    return None


def load_colors(assay: str) -> dict[str, str]:
    pc = pd.read_csv(PRED_CONFIG_DIR / f"pred_config_h3k27ac_counting_{assay}.tsv", sep="\t")
    colors = {}
    for _, row in pc.iterrows():
        key = _variant_key(row["pred_id"] + ".", ASSAYS[assay]["only_id"])
        if key is not None:
            colors[key] = row["color"]
    return colors


def load_overall(assay: str) -> pd.DataFrame:
    df = pd.read_csv(CRISPR_RESULTS / ASSAYS[assay]["comparison"] / "performance_summary.txt", sep="\t")
    df["variant"] = df["pred_uid"].apply(lambda x: _variant_key(x, ASSAYS[assay]["only_id"]))
    df = df[df["variant"].notna()].copy()
    df["variant"] = pd.Categorical(df["variant"], categories=VARIANT_ORDER, ordered=True)
    return df.sort_values("variant")


def load_distance_binned(assay: str) -> pd.DataFrame:
    df = pd.read_csv(CRISPR_RESULTS / ASSAYS[assay]["comparison"] / "performance_summary_distance_binned.txt", sep="\t")
    df = df[(df["metric"] == "auprc") & (df["subset"] != "All")].copy()
    df["variant"] = df["id"].apply(lambda x: _variant_key(x, ASSAYS[assay]["only_id"]))
    df = df[df["variant"].notna()].copy()
    df["variant"] = pd.Categorical(df["variant"], categories=VARIANT_ORDER, ordered=True)
    df["subset"] = pd.Categorical(df["subset"], categories=DIST_BIN_ORDER, ordered=True)
    return df.sort_values(["variant", "subset"])


def bar_panel(ax, df, value_col, lo_col, hi_col, colors, title):
    x = np.arange(len(VARIANT_ORDER))
    vals = df.set_index("variant")[value_col].reindex(VARIANT_ORDER)
    lo = df.set_index("variant")[lo_col].reindex(VARIANT_ORDER)
    hi = df.set_index("variant")[hi_col].reindex(VARIANT_ORDER)
    yerr = np.vstack([vals - lo, hi - vals])
    bar_colors = [colors.get(v, "#999999") for v in VARIANT_ORDER]
    ax.bar(x, vals, yerr=yerr, color=bar_colors, capsize=2,
           error_kw={"elinewidth": 0.75, "capthick": 0.75})
    ax.set_xticks(x)
    ax.set_xticklabels([VARIANT_LABELS[v] for v in VARIANT_ORDER], rotation=40, ha="right", fontsize=6)
    ax.set_title(title)


def fig_auprc(colors_by_assay):
    fig, axes = plt.subplots(1, 2, figsize=figsize(2, aspect=0.55), sharey=True)
    for i, assay in enumerate(ASSAYS):
        df = load_overall(assay)
        bar_panel(axes[i], df, "AUPRC", "AUPRC_lowerCi", "AUPRC_upperCi",
                  colors_by_assay[assay], ASSAYS[assay]["label"])
        add_panel_label(axes[i], "ab"[i])
    axes[0].set_ylabel("AUPRC (all E-G pairs)")
    fig.suptitle("H3K27ac counting geometry vs. CRISPR benchmark AUPRC", y=1.04)
    save_fig(fig, PLOTS_DIR / "fig1_headline_auprc")
    plt.close(fig)


def fig_precision(colors_by_assay):
    fig, axes = plt.subplots(1, 2, figsize=figsize(2, aspect=0.55), sharey=True)
    for i, assay in enumerate(ASSAYS):
        df = load_overall(assay)
        bar_panel(axes[i], df, "PrecMinSens", "PrecMinSens_lowerCi", "PrecMinSens_upperCi",
                  colors_by_assay[assay], ASSAYS[assay]["label"])
        add_panel_label(axes[i], "ab"[i])
    axes[0].set_ylabel("Precision at 70% sensitivity")
    fig.suptitle("H3K27ac counting geometry vs. precision at fixed sensitivity", y=1.04)
    save_fig(fig, PLOTS_DIR / "fig2_precision")
    plt.close(fig)


def fig_distance(colors_by_assay):
    fig, axes = plt.subplots(1, 2, figsize=figsize(2, aspect=0.5), sharey=True)
    for i, assay in enumerate(ASSAYS):
        df = load_distance_binned(assay)
        colors = colors_by_assay[assay]
        for variant in VARIANT_ORDER:
            sub = df[df["variant"] == variant].set_index("subset").reindex(DIST_BIN_ORDER)
            x = np.arange(len(DIST_BIN_ORDER))
            yerr = np.vstack([sub["full"] - sub["lower"], sub["upper"] - sub["full"]])
            axes[i].errorbar(x, sub["full"], yerr=yerr, marker="o", markersize=2.5,
                              color=colors.get(variant, "#999999"), linewidth=0.9,
                              capsize=2, elinewidth=0.6, capthick=0.6,
                              label=VARIANT_LABELS[variant].replace("\n", " "))
        axes[i].set_xticks(np.arange(len(DIST_BIN_ORDER)))
        axes[i].set_xticklabels([DIST_BIN_LABELS[b] for b in DIST_BIN_ORDER], fontsize=6)
        axes[i].set_title(ASSAYS[list(ASSAYS)[i]]["label"])
        add_panel_label(axes[i], "ab"[i])
    axes[0].set_ylabel("AUPRC")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=5)
    fig.suptitle("AUPRC by enhancer-gene distance", y=1.05)
    save_fig(fig, PLOTS_DIR / "fig3_distance_auprc")
    plt.close(fig)


def main():
    apply_rcparams()
    colors_by_assay = {assay: load_colors(assay) for assay in ASSAYS}
    fig_auprc(colors_by_assay)
    fig_precision(colors_by_assay)
    fig_distance(colors_by_assay)
    print(f"Wrote figures to {PLOTS_DIR}")


if __name__ == "__main__":
    main()
