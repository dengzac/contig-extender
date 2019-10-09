import sys
import extender
import numpy
import cython
import psutil
import argparse
import multiprocessing
from Bio import SeqIO
MIN_SCORE = 10000
MIN_SCORE2 = 5
MIN_OVERLAP = 15
BRANCH_LIMIT = 1
STOP_LENGTH = 250
NUM_THREADS = multiprocessing.cpu_count()
COMPLEX_THRESHOLD = 15
# extender.run_script(sys.argv[1:])
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("reference", help="fasta file with contigs to extend")
parser.add_argument("reads", help="fastq file with reads to extend with")

parser.add_argument("out", nargs="?", help="output directory", default='output')
parser.add_argument(
    "--min-overlap-length",
    help="minimum length of overlap between candidate read and contig",
    nargs="?",
    default=MIN_OVERLAP,
    type=int,
)
parser.add_argument(
    "--extend-tolerance",
    help="lower numbers require more reads to extend",
    nargs="?",
    default=2.5,
    type=float,
)

parser.add_argument(
    "--coverage",
    help="estimate of coverage",
    nargs="?",
    default=10,
    type=float,
)
parser.add_argument(
    "--min-branch-score",
    help="minimum score required to create alternative contig",
    nargs="?",
    default=MIN_SCORE2,
    type=int,
)
parser.add_argument(
    "--branch-limit",
    help="number of alternative contigs to output",
    nargs="?",
    default=BRANCH_LIMIT,
    type=int,
)
parser.add_argument(
    "--stop-length",
    help="terminate extension if substring of this size is repeated within contig",
    nargs="?",
    default=STOP_LENGTH,
    type=int,
)
parser.add_argument(
    "--threads",
    help="number of threads to use in computing alignments",
    nargs="?",
    default=psutil.cpu_count(),
)
parser.add_argument(
    "--complex-threshold",
    help="[0-100] higher values indicate less complexity. -1 to disable",
    nargs="?",
    default=COMPLEX_THRESHOLD,
    type=int
)
args = parser.parse_args(sys.argv[1:])
import base64
import copy
import os
import os.path
all_output = []
with open(args.reference) as multi_contigs:
    sequences = SeqIO.parse(multi_contigs, 'fasta')
    for seq in sequences:
        fname = base64.urlsafe_b64encode(str.encode(seq.id)).decode()
        with open(fname, 'w') as tempfile:
            tempfile.writelines(['>' + seq.id + '\n', str(seq.seq) + '\n'])
        tempargs = copy.copy(args)
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
        all_output.extend(branches)
if len(all_output) == 0:
    print("Empty input")
else:
    SeqIO.write(all_output, os.path.join(args.out, "contigs.fasta"), 'fasta')
    print("Output is located in " + os.path.join(args.out, "contigs.fasta"))
