#!/usr/bin/env python3
"""
cgenes_BLASTn_genome.py: A command-line wrapper for NCBI BLAST+ blastn.

This script executes a BLAST nucleotide search (blastn) against a local BLAST database,
adding a header row to the output TSV and logging the command and any stderr
to a timestamped log file. Users can specify the blastn binary location (defaults to
`blastn` in your PATH), query FASTA, and database, with an optional custom output filename.
"""
import argparse
import subprocess
import datetime
import os
import sys

# Fields to include in the BLAST TSV output format 6
OUTFMT_FIELDS = [
    'qseqid', 'sseqid', 'qacc', 'sacc',
    'qstart', 'qend', 'sstart', 'send',
    'sseq', 'length', 'mismatch', 'gaps',
    'evalue', 'bitscore'
]


def get_output_filename(query_path, db_name, custom_name=None):
    """Determine the output filename based on the query and database or use a custom name."""
    if custom_name:
        return custom_name
    query_base = os.path.splitext(os.path.basename(query_path))[0]
    return f"{query_base}_{db_name}.tsv"


def build_blastn_cmd(blastn_path, query_path, db_name, outfmt_fields):
    """Build the BLAST+ blastn command arguments list."""
    fields_str = ' '.join(outfmt_fields)
    return [
        blastn_path,
        '-query', query_path,
        '-db', db_name,
        '-outfmt', f"6 {fields_str}"
    ]


def run_blast(cmd, output_file, log_file, header):
    """Execute BLAST, stream results to output, and log command and stderr."""
    with open(output_file, 'w') as out_f, open(log_file, 'w') as log_f:
        # Log the BLAST command
        log_f.write(f"[{datetime.datetime.now()}] Command: {' '.join(cmd)}\n")
        # Start BLAST process
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Write header to output
        out_f.write(header + '\n')
        # Stream stdout only to output file
        for line in proc.stdout:
            out_f.write(line)
        proc.wait()
        # Capture and log stderr, if any
        stderr = proc.stderr.read()
        if stderr:
            log_f.write("\n[STDERR]\n")
            log_f.write(stderr)
    return proc.returncode


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Run blastn and format output TSV with header and logging.'
    )
    parser.add_argument(
        '--blastn',
        default='blastn',
        help='Path to blastn executable (defaults to blastn in your PATH)'
    )
    parser.add_argument(
        '--query',
        required=True,
        help='Path to query FASTA file'
    )
    parser.add_argument(
        '--db',
        required=True,
        help='Name of local BLAST database'
    )
    parser.add_argument(
        '-o', '--out',
        help='Custom output TSV filename'
    )
    return parser.parse_args()


def main():
    """Main entry point: parse args, build command, and run BLAST with logging."""
    args = parse_args()
    # Determine output filename
    output_file = get_output_filename(args.query, args.db, args.out)
    # Build header and command
    header = '\t'.join(OUTFMT_FIELDS)
    cmd = build_blastn_cmd(args.blastn, args.query, args.db, OUTFMT_FIELDS)
    # Prepare log filename with timestamp
    prefix = os.path.splitext(os.path.basename(output_file))[0]
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = f"{prefix}_{timestamp}.log"
    # Notify user
    print(f"Running BLAST: {' '.join(cmd)}")
    print(f"Writing results to {output_file} and log to {log_file}")
    # Execute BLAST
    rc = run_blast(cmd, output_file, log_file, header)
    if rc != 0:
        print(f"BLAST finished with return code {rc}. Check log file for details.", file=sys.stderr)
        sys.exit(rc)
    print("BLAST search completed successfully.")


if __name__ == '__main__':
    main()