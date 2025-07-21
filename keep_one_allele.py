#!/usr/bin/env python3
"""
This script filters an IMGT C gene FASTA file to retain only the first allele for each gene.  
It uses Biopython to parse input sequences, extracts gene names and allele numbers from record descriptions,  
and writes an output FASTA file containing only allele *01 sequences.  
The output file is named based on the input file with a '_first_alleles' suffix before the extension.  
The script accepts a single input FASTA file via command-line arguments.
"""

import argparse
import os
from Bio import SeqIO

def parse_args():
    """Parse command-line arguments to get the input FASTA file path."""
    parser = argparse.ArgumentParser(
        description="Filter IMGT C gene FASTA to include only first allele (*01) records"
    )
    parser.add_argument(
        "input_fasta",
        help="Path to the input IMGT C gene FASTA file"
    )
    return parser.parse_args()

def get_gene_allele(description):
    """Extract the gene name and allele number from a FASTA record description."""
    # description fields are 'IMGT000001|IGHA*01|...'; we want the second field
    parts = description.split("|")
    gene_field = parts[1]             # e.g. "IGHA*01"
    gene, allele = gene_field.split("*")
    return gene, allele

def filter_first_alleles(records):
    """Return only records corresponding to the first allele (*01) for each gene."""
    seen_genes = set()
    for record in records:
        gene, allele = get_gene_allele(record.description)
        if allele == "01" and gene not in seen_genes:
            seen_genes.add(gene)
            yield record

def main():
    """Main function to parse arguments, filter records, and write output FASTA."""
    args = parse_args()
    # Build output filename by appending '_first_alleles' before the extension
    base, ext = os.path.splitext(args.input_fasta)
    output_path = f"{base}_first_alleles{ext}"
    # Read, filter, and write
    records = SeqIO.parse(args.input_fasta, "fasta")
    first_alleles = filter_first_alleles(records)
    SeqIO.write(first_alleles, output_path, "fasta")
    print(f"Written first-allele records to '{output_path}'")

if __name__ == "__main__":
    main()
