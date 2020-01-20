"""Conting Extender Script

This script allows the ContigExtender software to be run from the command line.
Run `extender_wrapper --help` to see full list of options
"""

import sys

import argparse
import base64
import copy
import os
import os.path
from Bio import SeqIO
import psutil
import extender

MIN_SCORE = 10000
MIN_SCORE2 = 5
MIN_OVERLAP = 15
BRANCH_LIMIT = 1
STOP_LENGTH = 250
COMPLEX_THRESHOLD = 15

PARSER = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

PARSER.add_argument("reference", help="fasta file with contigs to extend")
PARSER.add_argument("reads", help="fastq file with reads to extend with")

PARSER.add_argument("out", nargs="?", help="output directory", default='output')
PARSER.add_argument(
    "--min-overlap-length",
    help="minimum length of overlap between candidate read and contig",
    nargs="?",
    default=MIN_OVERLAP,
    type=int,
)
PARSER.add_argument(
    "--extend-tolerance",
    help="lower numbers require more reads to extend",
    nargs="?",
    default=2.5,
    type=float,
)

PARSER.add_argument(
    "--coverage",
    help="estimate of coverage",
    nargs="?",
    default=10,
    type=float,
)
PARSER.add_argument(
    "--min-branch-score",
    help="minimum score required to create alternative contig",
    nargs="?",
    default=MIN_SCORE2,
    type=int,
)
PARSER.add_argument(
    "--branch-limit",
    help="number of alternative contigs to output",
    nargs="?",
    default=BRANCH_LIMIT,
    type=int,
)
PARSER.add_argument(
    "--stop-length",
    help="terminate extension if substring of this size is repeated within contig",
    nargs="?",
    default=STOP_LENGTH,
    type=int,
)
PARSER.add_argument(
    "--threads",
    help="number of threads to use in computing alignments",
    nargs="?",
    default=psutil.cpu_count(),
)
PARSER.add_argument(
    "--complex-threshold",
    help="[0-100] higher values indicate less complexity. -1 to disable",
    nargs="?",
    default=COMPLEX_THRESHOLD,
    type=int
)
ARGS = PARSER.parse_args(sys.argv[1:])

ALL_OUTPUT = []
with extender.open_file(ARGS.reference) as multi_contigs:
    SEQUENCES = SeqIO.parse(multi_contigs, 'fasta')
    for seq in SEQUENCES:
        fname = base64.urlsafe_b64encode(str.encode(seq.id)).decode()
        with open(fname, 'w') as tempfile:
            tempfile.writelines(['>' + seq.id + '\n', str(seq.seq) + '\n'])
        tempargs = copy.copy(ARGS)
        tempargs.reference = fname
        output_dir = os.path.join(extender.run_script(tempargs), 'contigs')
        branches = []
        num = 1
        for contig in os.listdir(output_dir):
            contig = os.path.join(output_dir, contig)
            seq = list(SeqIO.parse(open(contig), 'fasta'))[0]
            seq.id = seq.id + '_' + str(num)
            branches.append(seq)
            num += 1
        branches = sorted(branches, key=lambda x: len(str(x.seq)), reverse=True)
        ALL_OUTPUT.extend(branches)
if not ALL_OUTPUT:
    print("Empty input")
else:
    SeqIO.write(ALL_OUTPUT, os.path.join(ARGS.out, "contigs.fasta"), 'fasta')
    print("Output is located in " + os.path.join(ARGS.out, "contigs.fasta"))
