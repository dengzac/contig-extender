"""Contig Extender Script

This script allows the ContigExtender software to be run from the command line.
Run `extender_wrapper --help` to see full list of options
"""

import sys
import argparse
import base64
import copy
import os
import os.path
import tempfile
from Bio import SeqIO
import psutil
import extender

# Default parameters
MIN_SCORE = 10000
MIN_SCORE2 = 5
MIN_OVERLAP = 15
BRANCH_LIMIT = 1
STOP_LENGTH = 250
COMPLEX_THRESHOLD = 15
INS_LENGTH = 500
PAIRED = False

class ArgumentParser(argparse.ArgumentParser):
    """Prints help message along with default error message when invalid arguments are given"""
    def error(self, message):
        self.print_help(sys.stderr)
        super().error(message)

def parse_args(argv):
    PARSER = ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    PARSER.add_argument("reference", help="fasta file with contigs to extend")
    PARSER.add_argument(
        "reads",
        nargs="?",
        help="fastq file with unpaired reads to extend with; use the --m1 and --m2 options instead for paired reads",
    )
    PARSER.add_argument("--m1", nargs="?", help="fastq file with #1 mates of paired data")
    PARSER.add_argument("--m2", nargs="?", help="fastq file with #2 mates of paried data")

    PARSER.add_argument("--out", nargs="?", help="output directory", default="output")

    PARSER.add_argument(
        "--enable-pair-constraint",
        dest="pair_constraint",
        action="store_true",
        default=False,
        help="Require paired-end alignments to satisfy orientation and insert size constraints according to Bowtie2. Otherwise, paired data will be treated as two sets of unpaired reads.",
    )

    PARSER.add_argument(
        "--min-overlap-length",
        help="minimum length of overlap between candidate read and the existing contig",
        nargs="?",
        default=MIN_OVERLAP,
        type=int,
    )
    PARSER.add_argument(
        "--extend-tolerance",
        help="This parameter, along with read length and coverage, is used to determine the required score to extend. Lower numbers require better quality alignments to extend",
        nargs="?",
        default=2.5,
        type=float,
    )

    PARSER.add_argument(
        "--coverage",
        help="estimate of sequencing coverage",
        nargs="?",
        default=10,
        type=float,
    )
    PARSER.add_argument(
        "--min-branch-score",
        help="A new branch for an alternative contig will be created if it meets this threshold. Then, both the main and alternative branches will be processed recursively.",
        nargs="?",
        default=MIN_SCORE2,
        type=int,
    )
    PARSER.add_argument(
        "--branch-limit",
        help="This limits the number of alternative contigs that are considered and output. Once the search tree has this many leaf nodes, no new branches will be created",
        nargs="?",
        default=BRANCH_LIMIT,
        type=int,
    )
    PARSER.add_argument(
        "--stop-length",
        help="Terminate extension if any substring of this size is repeated within the contig. This prevents circular genomes from infinite extension.",
        nargs="?",
        default=STOP_LENGTH,
        type=int,
    )
    PARSER.add_argument(
        "--threads",
        help="Number of threads to use in computing alignments",
        nargs="?",
        default=psutil.cpu_count(),
    )
    PARSER.add_argument(
        "--complex-threshold",
        help="[0-100] This parameter is passed to PRINSEQ's DUST complexity filter, is run on the reads. Higher values indicate less complexity. -1 to disable",
        nargs="?",
        default=COMPLEX_THRESHOLD,
        type=int,
    )
    PARSER.add_argument(
        "--maxins",
        help="When using paired-end constraints, this is the maximum fragment length, passed directly to Bowtie2. This length includes both reads ond the gap between them. Refer to the Bowtie2 manual for more information",
        nargs="?",
        default=INS_LENGTH,
        type=int,
    )
    return PARSER.parse_args(argv)

if __name__ == "__main__":
    ARGS = parse_args(sys.argv[1:])
    # No unpaired reads specified, check for paired files
    combined_reads, filename = tempfile.mkstemp(dir=os.getcwd())
    print(filename)
    if not ARGS.reads:
        ARGS.paired = True
        combined_reads_f = os.fdopen(combined_reads, "w+")

        try:
            m1 = extender.open_file(ARGS.m1)
            m2 = extender.open_file(ARGS.m2)

            while True:
                l1 = []
                l2 = []

                # Read one fastq record from each file
                for line_num in range(4):
                    l1.append(next(m1, ""))
                    l2.append(next(m2, ""))

                # Check if both files done, and raise exception if only one is
                if l1[0] == "" and l2[0] == "":
                    break
                if l1[0] == "" or l2[0] == "":
                    raise RuntimeError("Different numbers of #1 and #2 mates provided")
                l2[0] = l1[0].strip() + "/2\n"
                l1[0] = l1[0].strip() + "/1\n"
                combined_reads_f.write("".join(l1))
                combined_reads_f.write("".join(l2))
            ARGS.reads = filename
        except IOError:
            raise RuntimeError("No read files specified")
    else:
        ARGS.paired = False

    # Run in unpaired mode if constraints not enabled
    if not ARGS.pair_constraint:
        ARGS.paired = False
    ALL_OUTPUT = []

    with extender.open_file(ARGS.reference) as multi_contigs:
        SEQUENCES = SeqIO.parse(multi_contigs, "fasta")
        for seq in SEQUENCES:
            # Extract each input contig and pass to extender
            fname = base64.urlsafe_b64encode(str.encode(seq.id)).decode()
            with open(fname, "w") as tempfile:
                tempfile.writelines([">" + seq.id + "\n", str(seq.seq) + "\n"])
            tempargs = copy.copy(ARGS)
            tempargs.reference = fname
            print(tempargs)
            output_dir = os.path.join(extender.run_script(tempargs), "contigs")
            branches = []
            num = 1
            # Output each extended contig, in order of size
            for contig in os.listdir(output_dir):
                contig = os.path.join(output_dir, contig)
                seq = list(SeqIO.parse(open(contig), "fasta"))[0]
                seq.id = seq.id + "_" + str(num)
                branches.append(seq)
                num += 1
            branches = sorted(branches, key=lambda x: len(str(x.seq)), reverse=True)
            ALL_OUTPUT.extend(branches)
            os.remove(fname)

    if not ALL_OUTPUT:
        print("Empty input")
    else:
        SeqIO.write(ALL_OUTPUT, os.path.join(ARGS.out, "contigs.fasta"), "fasta")
        print("Output is located in " + os.path.join(ARGS.out, "contigs.fasta"))

    # Cleanup
    os.close(combined_reads)
    os.remove(filename)
