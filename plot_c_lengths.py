#!/usr/bin/env python3
"""
Read one or more FASTA files containing immunoglobulin C gene sequences (one species per file)
and produce a clustered horizontal bar plot of sequence lengths. Each cluster corresponds to an
isotype (IGHM, IGHD, IGHG, IGHE, IGHA, plus any others), colored by species, and bars represent
every individual gene, labelled by gene name and species abbreviation. Uses Biopython to parse FASTA,
pandas to organize the data, and seaborn for theming and palette. Saves a high-resolution PNG (600 DPI),
named by either the user-supplied prefix or the combined input basenames.
"""

import argparse
import os
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def abbreviate_species(species_str):
    """Abbreviate a species string by taking first letters of each word before any underscore."""
    generic = species_str.split('_')[0]
    return ''.join(w[0].upper() for w in generic.split())


def parse_sequences(fasta_file):
    """Parse immunoglobulin C gene records from a FASTA, yielding dicts with isotype, gene, species, length."""
    for rec in SeqIO.parse(fasta_file, "fasta"):
        fields = rec.description.split('|')
        if len(fields) < 3:
            continue
        gene = fields[1]
        isotype = gene[:4]
        species_abbrev = abbreviate_species(fields[2])
        yield {
            'isotype': isotype,
            'gene': gene,
            'species': species_abbrev,
            'length': len(rec.seq)
        }


def plot_sequence_lengths(df, outfile):
    """Manually draw every gene as its own bar within isotype clusters, colored by species."""
    # 1) Determine isotype order
    common = ["IGHM", "IGHD", "IGHG", "IGHE", "IGHA"]
    present = [i for i in common if i in df['isotype'].unique()]
    extra   = sorted(i for i in df['isotype'].unique() if i not in common)
    iso_order = present + extra

    # 2) Species palette
    species_order = sorted(df['species'].unique())
    palette = dict(zip(species_order,
                       sns.color_palette("colorblind", n_colors=len(species_order))))

    # 3) Compute bar geometry
    n_clusters = len(iso_order)
    # max number of genes in any one isotype → use to set bar thickness
    max_bars = df.groupby('isotype').size().reindex(iso_order, fill_value=0).max()
    bar_h = 0.8 / max_bars
    centers = list(range(n_clusters))

    # 4) Setup figure
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, max(6, n_clusters * 0.6)))

    # 5) Draw each gene
    for i, iso in enumerate(iso_order):
        sub = df[df['isotype'] == iso]
        for j, row in sub.reset_index(drop=True).iterrows():
            # center offset so bars straddle the integer y‑tick
            offset = (-(len(sub) - 1) / 2 + j) * bar_h
            y = centers[i] + offset
            ax.barh(y,
                    row['length'],
                    height=bar_h,
                    color=palette[row['species']])
            ax.text(row['length'],
                    y,
                    f"{row['gene']}_{row['species']}",
                    va='center',
                    ha='left',
                    fontsize=6)

    # 6) Final formatting
    ax.set_yticks(centers)
    ax.set_yticklabels(iso_order)
    ax.invert_yaxis()
    ax.set_xlabel("Sequence Length (nt)")
    ax.set_ylabel("Isotype")

    # species legend
    handles = [plt.Line2D([0], [0], color=palette[s], lw=6) for s in species_order]
    ax.legend(handles=handles, labels=species_order,
              title="Species", bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(outfile, dpi=600)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Clustered bar plot of immunoglobulin C gene lengths by isotype, colored by species"
    )
    parser.add_argument('fastas', nargs='+',
                        help='Input FASTA files (one species per file)')
    parser.add_argument('-p', '--prefix', default=None,
                        help='Output file prefix (default: combined input basenames)')
    args = parser.parse_args()

    # collect all records into a DataFrame
    records = []
    for fasta in args.fastas:
        records.extend(parse_sequences(fasta))
    df = pd.DataFrame(records)

    # determine output file name
    bases = [os.path.splitext(os.path.basename(f))[0] for f in args.fastas]
    out_prefix = args.prefix if args.prefix else "_".join(bases)
    out_file = f"{out_prefix}.png"

    # plot and save
    plot_sequence_lengths(df, out_file)
    print(f"Saved plot to {out_file}")


if __name__ == "__main__":
    main()
