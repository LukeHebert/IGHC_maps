#!/usr/bin/env python3
"""
This script filters an IMGT C gene FASTA file to retain only the first non-partial allele for each gene.  
It uses Biopython to parse input sequences, extracts gene names from record descriptions,  
skips any record whose description contains "partial", and writes an output FASTA file  
containing only the first complete allele for each gene. The output file is named based  
on the input file with a '_first_complete_alleles' suffix before the extension.
"""

import argparse
import os
from Bio import SeqIO

def parse_args():
    """Parse command-line arguments to get the input FASTA file path."""
    parser = argparse.ArgumentParser(
        description="Filter IMGT C gene FASTA to include only the first non-partial allele per gene"
    )
    parser.add_argument(
        "input_fasta",
        help="Path to the input IMGT C gene FASTA file"
    )
    return parser.parse_args()

def get_gene_allele(description):
    """Extract the gene name and allele number from a FASTA record description."""
    parts = description.split("|")
    gene_field = parts[1]             # e.g. "IGHA*01"
    gene, allele = gene_field.split("*")
    return gene, allele

def filter_first_complete(records):
    """Yield the first record for each gene whose description does not include 'partial'."""
    seen_genes = set()
    for record in records:
        gene, _ = get_gene_allele(record.description)
        if gene in seen_genes:
            continue
        if "partial" in record.description.lower():
            continue
        seen_genes.add(gene)
        yield record

def main():
    """Main function to parse arguments, filter records, and write output FASTA."""
    args = parse_args()
    # Build output filename by appending '_first_complete_alleles' before the extension
    base, ext = os.path.splitext(args.input_fasta)
    output_path = f"{base}_first_complete_alleles{ext}"
    # Read, filter, and write
    records = SeqIO.parse(args.input_fasta, "fasta")
    complete_first = filter_first_complete(records)
    SeqIO.write(complete_first, output_path, "fasta")
    print(f"Written first-complete-allele records to '{output_path}'")

if __name__ == "__main__":
    main()
