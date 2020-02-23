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
import tempfile

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
PARSER.add_argument("reads", nargs="?",
                    help="fastq file with unpaired reads to extend with")
PARSER.add_argument("--m1", nargs="?", help="fastq file with #1 mates")
PARSER.add_argument("--m2", nargs="?", help="fastq file with #2 mates")

PARSER.add_argument("--out", nargs="?",
                    help="output directory", default='output')

PARSER.add_argument("--enable-pair-constraint",
                    dest='pair_constraint', action='store_true', default=False, help="Require paired-end alignments to satisfy orientation and insert size constraints")

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

# No unpaired reads specified, check for paired files
combined_reads, filename = tempfile.mkstemp(dir=os.getcwd())
print(filename)
if not ARGS.reads:
    ARGS.paired = True
    combined_reads_f = os.fdopen(combined_reads, 'w+')

    try:
        m1 = extender.open_file(ARGS.m1)
        m2 = extender.open_file(ARGS.m2)

        while True:
            l1 = []
            l2 = []

            # Read one fastq record from each file
            for line_num in range(4):
                l1.append(next(m1, ''))
                l2.append(next(m2, ''))

            # Check if both files done, and raise exception if only one is
            if l1[0] == '' and l2[0] == '':
                break
            elif l1[0] == '' or l2[0] == '':
                raise RuntimeError(
                    "Different numbers of #1 and #2 mates provided")
            l2[0] = l1[0].strip() + "/2\n"
            l1[0] = l1[0].strip() + "/1\n"
            combined_reads_f.write(''.join(l1))
            combined_reads_f.write(''.join(l2))
        ARGS.reads = filename
    except IOError:
        raise RuntimeError("No read files specified")
else:
    ARGS.paired = False
if not ARGS.pair_constraint:
    ARGS.paired = False
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
        branches = sorted(branches, key=lambda x: len(
            str(x.seq)), reverse=True)
        ALL_OUTPUT.extend(branches)
        os.remove(fname)
if not ALL_OUTPUT:
    print("Empty input")
else:
    SeqIO.write(ALL_OUTPUT, os.path.join(ARGS.out, "contigs.fasta"), 'fasta')
    print("Output is located in " + os.path.join(ARGS.out, "contigs.fasta"))

# Cleanup
os.close(combined_reads)
os.remove(filename)
