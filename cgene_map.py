#!/usr/bin/env python3
"""
Generate gene locus maps from BLAST hit TSVs of immunoglobulin C gene loci across species.
Reads one or more TSVs of BLAST hits, extracts species and gene labels from qseqid, selects the most-hit scaffold per species,
orienting each species so IGHM (or fallback to IGHD) is on the left if present, normalizes lengths to the longest scaffold, and color‐codes gene families.
Creates one subplot per species (ordered Homo, Canis, Neogale, Mustela), labels scaffold coordinate ends, labels each subplot with species and scaffold name,
places vertical gene labels directly atop each gene, and saves a high‐resolution PNG output.
"""

import os
import re
import argparse
import warnings
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# suppress NumPy longdouble signature warnings
warnings.filterwarnings(
    "ignore",
    message=r".*Signature.*for <class 'numpy.longdouble'>.*",
    module="numpy._core.getlimits"
)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot immunoglobulin C locus maps from one or more BLAST TSVs.")
    parser.add_argument(
        "tsvs",
        nargs='+',
        help="Input TSV file(s) of BLAST hits"
    )
    return parser.parse_args()


def load_data(tsv_path):  # noqa: E501
    """Load a BLAST hit TSV into a DataFrame."""
    return pd.read_csv(tsv_path, sep="\t", dtype={"qseqid": str})


def extract_metadata(df):
    """Extract species and gene labels from qseqid column."""
    df = df.copy()
    df["species"] = df["qseqid"].str.split("|").str[2]
    df["gene"] = df["qseqid"].str.split("|").str[1].str.split("*").str[0]
    df["start_abs"] = df[["sstart", "send"]].min(axis=1)
    df["end_abs"] = df[["sstart", "send"]].max(axis=1)
    return df


def select_scaffold(df):
    """For each species, keep only hits on its most frequently hit scaffold and store scaffold metadata."""
    selected = {}
    for sp, sub in df.groupby("species"):
        top_scaffold = sub["sacc"].value_counts().idxmax()
        sp_df = sub[sub["sacc"] == top_scaffold].copy()
        sp_df.attrs["scaffold_min"] = sp_df["start_abs"].min()
        sp_df.attrs["scaffold_max"] = sp_df["end_abs"].max()
        sp_df.attrs["scaffold_name"] = top_scaffold
        selected[sp] = sp_df
    return selected


def orient_and_normalize(selected):
    """Orient each species so IGHM or IGHD is on the left (if present) and compute relative positions."""
    oriented = {}
    for sp, df in selected.items():
        Smin = df.attrs["scaffold_min"]
        Smax = df.attrs["scaffold_max"]
        L = Smax - Smin
        df["rel_start"] = df["start_abs"] - Smin
        df["rel_end"]   = df["end_abs"]   - Smin
        # choose orientation gene: IGHM preferred, then IGHD
        orient_gene = None
        if (df["gene"] == "IGHM").any():
            orient_gene = "IGHM"
        elif (df["gene"] == "IGHD").any():
            orient_gene = "IGHD"
        # if orientation gene available, flip if needed
        if orient_gene:
            hit = df[df["gene"] == orient_gene].iloc[0]
            center = (hit["rel_start"] + hit["rel_end"]) / 2
            if center > L / 2:
                df["rel_start"], df["rel_end"] = (
                    L - df["rel_end"],
                    L - df["rel_start"]
                )
        # assign for plotting
        df["plot_start"] = df["rel_start"]
        df["plot_end"]   = df["rel_end"]
        df.attrs["length"] = L
        oriented[sp] = df
    return oriented


def assign_colors(oriented):
    """Assign colorblind‐friendly shades to genes sharing family prefixes."""
    all_genes = pd.concat(oriented.values())["gene"].unique()
    families = {}
    for g in sorted(all_genes):
        prefix = re.match(r"([A-Za-z]+)", g).group(1)
        families.setdefault(prefix, []).append(g)
    base_colors = sns.color_palette("colorblind", n_colors=len(families))
    family_to_base = dict(zip(families.keys(), base_colors))
    gene_colors = {}
    for fam, genes in families.items():
        base = family_to_base[fam]
        n = len(genes)
        for i, g in enumerate(sorted(genes, key=lambda x: int(re.sub(r"\D", "", x)) if re.search(r"\d+$", x) else 0)):
            if n > 1:
                factor = 1 - 0.5 * (i / (n - 1))
                shade = tuple(base[j] * factor + (1 - factor) for j in range(3))
            else:
                shade = base
            gene_colors[g] = shade
    for df in oriented.values():
        df["color"] = df["gene"].map(gene_colors)
    return oriented


def plot_loci(oriented, output_path):
    """Plot one horizontal locus map per species with titles, vertical gene labels, and save to PNG."""
    desired = ["Homo", "Canis", "Neogale", "Mustela"]
    species = [sp for sp in desired if sp in oriented] + [sp for sp in oriented if sp not in desired]
    n = len(species)
    global_max = max(df.attrs["length"] for df in oriented.values())
    fig, axes = plt.subplots(n, 1, figsize=(8, 2 * n), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, sp in zip(axes, species):
        df = oriented[sp]
        scaffold = df.attrs.get("scaffold_name", "")
        # draw genes
        for _, row in df.iterrows():
            ax.add_patch(Rectangle(
                (row["plot_start"], 0.1),
                row["plot_end"] - row["plot_start"],
                0.6,
                color=row["color"], ec="black"
            ))
        # vertical gene labels directly on top
        for _, row in df.iterrows():
            x_center = (row["plot_start"] + row["plot_end"]) / 2
            ax.text(
                x_center, 0.75, row["gene"],
                ha="center", va="bottom",
                rotation=90, fontsize=8,
                transform=ax.get_xaxis_transform()
            )
        # title above plot
        ax.set_title(f"{sp} — {scaffold}", loc='left', pad=5, y=1.15)
        ax.set_xlim(0, global_max)
        ax.set_ylim(0, 1.2)
        ax.axis("off")
        # scaffold end labels
        ax.text(0, -0.1, str(df.attrs["scaffold_min"]), ha="center", va="top", transform=ax.get_xaxis_transform())
        ax.text(df.attrs["length"], -0.1, str(df.attrs["scaffold_max"]), ha="center", va="top", transform=ax.get_xaxis_transform())
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    """Orchestrate loading, processing, plotting, and saving of gene locus maps."""
    args = parse_args()
    df_list = [load_data(tsv) for tsv in args.tsvs]
    df = pd.concat(df_list, ignore_index=True)
    df = extract_metadata(df)
    selected = select_scaffold(df)
    oriented = orient_and_normalize(selected)
    colored = assign_colors(oriented)
    if len(args.tsvs) == 1:
        base = os.path.splitext(args.tsvs[0])[0]
    else:
        base = os.path.commonprefix([os.path.basename(tsv) for tsv in args.tsvs]).rstrip("_-") or "combined"
    output_file = f"{base}_locus_map.png"
    plot_loci(colored, output_file)
    print(f"Saved locus map to {output_file}")

if __name__ == "__main__":
    main()