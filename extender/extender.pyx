import time
import pathlib
import argparse
import tempfile
import shutil
import subprocess
import os
import re
import numpy as np
import heapq
import os.path
import psutil
import gzip
from distutils.version import LooseVersion
from mimetypes import guess_type
import sys
from pathlib import Path
import multiprocessing
import cython
import array
import itertools
from Bio import SeqIO
cimport numpy as np
from libc.stdlib cimport malloc, free

cdef int get_char_pos(char c):
    if c == 'A':
        return 0
    elif c == 'C':
        return 1
    elif c == 'G':
        return 2
    elif c == 'T':
        return 3

def replace_string(string, replace_letter, startPos, endPos):
    return (
        string[0:startPos]
        + (replace_letter * (endPos - startPos + 1))
        + string[endPos + 1:]
    )

cdef int c_max(int a, int b, int c):
    cdef int best = 0
    if a > best:
        best = a
    if b > best:
        best = b
    if c > best:
        best = c
    if a == best:
        return a
    if b == best:
        return b
    if c == best:
        return c
    
def open_file(fname):
    encoding = guess_type(fname)[1]
    return gzip.open(fname, "rt") if encoding == 'gzip' else open(fname)

def check_repeats(contig, length):
    return contig[-length:] in contig[:-length]


cdef char* get_possible_consensus(char[:] read, int pos, int[:] ambigPos, int ambigN):

    cdef char* c_string = <char*>malloc((ambigN+1)*sizeof(char))
    s = ""
    cdef int[:] ambig_arr
    cdef int i, readIndex

    ambig_arr = ambigPos

    cdef int n = len(read)
    cdef int idx
    for idx in range(ambigN):
        i = ambig_arr[idx]
        readIndex = i - pos
        if readIndex < 0 or readIndex >= n:
            c_string[idx] = 'N'
        else:
            c_string[idx] = read[readIndex]
    return c_string

def prepare_directory(name):
    try:
        shutil.rmtree(name)
    except:
        pass
    os.makedirs(name)

cdef int MIN_SCORE = 10000
cdef int MIN_SCORE2 = 5
cdef int MIN_OVERLAP = 15
cdef int BRANCH_LIMIT = 1
cdef int STOP_LENGTH = 250
cdef int NUM_THREADS = multiprocessing.cpu_count()
cdef int COMPLEX_THRESHOLD = 15
cdef int maxlength
paired = False
dirpath = ''
stack = []
maxlength = 0

def check_bowtie_ver(path):
    version = subprocess.run(
            [
                path,
                "--version"
            ],
            stdout=subprocess.PIPE
        )

    try:
        version_string = version.stdout.decode().strip().split('\n')[0].split(' ')[2]
    except:
        pass
    else:
        if LooseVersion(version_string)<LooseVersion("2.3.4"):
            raise RuntimeError("Current bowtie2 version " + version_string + " too old; must be >= 2.3.4")
def get_bowtie_path():
    installed_path = shutil.which("bowtie2")
    if installed_path:
        check_bowtie_ver(installed_path)
        return installed_path
    elif os.name == 'nt' and hasattr(sys, '_MEIPASS'):
        p_path = os.path.join(sys._MEIPASS, 'bowtie2.bat')
        check_bowtie_ver(p_path)
        return p_path
    else:
        raise RuntimeError("bowtie2 not found. Add the executable location to PATH")

def get_bowtie_build_path():
    # Replace last occurence of "bowtie2" with "bowtie2-build"
    li = get_bowtie_path().rsplit("bowtie2", 1)
    return "bowtie2-build".join(li)
def iterate(
    contig,
    inputReads,
    outFile,
    quiet=True,
    exclude=True,
    curPath="",
    allow_alt=True,
    used_reads=set(),
    num_branches=1,
    analyze_mode=False,
    args=[]
):

    cdef int idx
    global dirpath
    global stack
    global maxlength

    with open_file(contig) as reference:
        data = list(SeqIO.parse(reference, 'fasta'))
        referenceData = ['>' + data[0].id, str(data[0].seq)]
        referenceData[1] = referenceData[1].strip("N")

    extendedReference = (
        "A" + ("N" * (maxlength + 1)) + referenceData[1] + ("N" * (maxlength + 1)) + "A"
    )
    origReference = extendedReference
    # Keep only edges of reference sequence in unpaired mode
    if exclude and not PAIRED:
        extendedReference = replace_string(
            extendedReference,
            "N",
            2 + 2 * maxlength + 1*maxlength,
            len(extendedReference) - 2 - 2 * maxlength - 1*maxlength,
        )
    extendedReferenceFile = dirpath + "/extendedRef.fa"

    with open(extendedReferenceFile, "w") as out:
        out.writelines([referenceData[0] + "\n", extendedReference + "\n"])

    with open(outFile+"_raw", "w") as out:
        out.writelines([referenceData[0] + "\n", extendedReference + "\n"])
    if not quiet:
        print("Building index")
    subprocess.run(
        [get_bowtie_build_path(), extendedReferenceFile, dirpath + "/ref"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if not quiet:
        print("Finding matches")

    bamOutput = dirpath + "/out2.sam"
    try: 

        result = subprocess.run(
            [
                get_bowtie_path(),
                "-x",
                dirpath + "/ref",
                "--interleaved" if PAIRED else "-U",
                "-",
                "-L",
                "20",
                "--n-ceil",
                "L,1000,0",
                "--np",
                "0",
                "-p",
                str(NUM_THREADS),
                "-i",
                "S,10,0",
                "--dpad",
                str(maxlength),
                "--gbar",
                "20",
                "--rdg",
                "100000,100000",
                "--rfg",
                "1000000,1000000",
                "--no-unal",
            ],
            stdout=subprocess.PIPE,
            input=inputReads,
            stderr=subprocess.PIPE
        )
    except Exception as e:
        print(e)
        raise e
    if not quiet:
        print("Building consensus")
    frequencies = np.ones((len(extendedReference), 4), dtype=np.int_)
    mult_frequencies = np.zeros((len(extendedReference), 4), dtype=np.int_)
    cdef long[:, :] freq_view = frequencies
    cdef long[:, :] mult_view = mult_frequencies

    mainconsensus = ""
    altconsensus = ""
    alignments = result.stdout.decode().strip().split("\n")
    splitter = re.compile(r"\t+")
    if not quiet:
        print('filter')
    alignments = filter(lambda a: a[0] != "@" and a[3] != "1", alignments)
    if not quiet:
        print('split')
    alignment_fields = [[x for x in re.split(splitter, o)] for o in alignments]

    if not quiet:
        print('pos')

    alignPosList = [int(a[3])-1 for a in alignment_fields]
    cdef long[:] alignPos = np.array(alignPosList)
    cdef int pos, quality, quality2
    cdef char n_char = 'N'
    cdef int c_pos, char_pos
    cdef char[:] r, ext
    cdef int refLength = len(referenceData[1])
    cdef int ii, rLen
    for ii in range(len(alignment_fields)):
        alignment_fields[ii][9] = bytearray(alignment_fields[ii][9].encode())
    ext = bytearray(extendedReference.encode())

    if not quiet:
        print("Computing scores")
    cdef int total_coverage = 0
    with cython.boundscheck(False):
        for ii in range(len(alignment_fields)):
            fields = alignment_fields[ii]

            # Ignore discordant pairs
            if "YT:Z:DP" in fields:
                continue
            elif "YT:Z:CP" in fields:
                pass

            rLen = len(fields[9])

            total_coverage += rLen

            pos = alignPos[ii]

            read = fields[9]
            r = fields[9]
            quality1 = 0
            quality2 = 0

            if pos < (maxlength + 2):
                quality1 = rLen - (maxlength + 2) + pos
            elif pos + rLen > maxlength + 2 + refLength:
                quality2 = maxlength + 2 + refLength - pos
            # Completely within reference

            # else:
            #     for idx in range(rLen):
            #         if r[idx] == ext[idx + pos]:
            #             quality2 += 1
            quality = c_max(quality1, quality2, MIN_OVERLAP)
            if quality == MIN_OVERLAP:
                quality = 0

            for idx in range(rLen):
                if r[idx] != n_char:
                    c_pos = idx + pos

                    char_pos = get_char_pos(r[idx])
                    freq_view[c_pos, char_pos] += quality * quality
                    mult_view[c_pos, char_pos] += 1
    if analyze_mode:
        return total_coverage
    keys = ["A", "C", "G", "T"]
    frequencies = [{"A": x[0], "C": x[1], "G": x[2], "T": x[3]} for x in frequencies]
    mult_frequencies = [{"A": x[0], "C": x[1], "G": x[2], "T": x[3]} for x in mult_frequencies]

    if not quiet:
        print("Building main consensus")
    mainconsensus = list(origReference)

    for i in itertools.chain(range(2 + 2 * maxlength), range(len(extendedReference) - 2 - 2 * maxlength, len(extendedReference))):
        largestVal = sum(frequencies[i].values())
        if largestVal > MIN_SCORE:

            mainconsensus[i] = max(frequencies[i], key=lambda k: frequencies[i][k])

    altconsensus = mainconsensus.copy()
    if not quiet:
        print(len(altconsensus))
    if not quiet:
        print('finding ambiguous')
    ambiguous = []
    for i in range(2*maxlength + 2):
        sortedfreq = heapq.nlargest(
            2, enumerate(mult_frequencies[i].values()), key=lambda x: x[1]
        )
        sortedfreq2 = heapq.nlargest(
            2, enumerate(frequencies[i].values()), key=lambda x: x[1]
        )
        largestVal = sortedfreq[1][1]

        if largestVal > MIN_SCORE2 and sortedfreq2[1][1] > sum(frequencies[i].values())/3:
            ambiguous.append(i)
    ambiguous.sort()
    cdef int[:] ambig_arr = np.array(ambiguous, dtype=np.int32)
    cdef int ambig_size = len(ambiguous)
    scores = {}
    if not quiet:
        print('finding candidates')

    for ii in range(len(alignment_fields)):
        fields = alignment_fields[ii]
        pos = alignPos[ii]

        c_string = get_possible_consensus(fields[9], pos, ambig_arr, ambig_size)
        candidate = unicode(c_string)
        free(c_string)

        if candidate in scores:
            scores[candidate] += 1
        else:
            scores[candidate] = 1

    if not quiet:
        print('altconsensus part 1')
    sorted_candidates = sorted(scores.items(), key=lambda k: k[1])
    sorted_candidates = list(
        filter(lambda x: len(x[0].replace("N", "")) > 0, sorted_candidates)
    )

    if len(sorted_candidates) > 1:
        for indx, pos in enumerate(ambiguous):
            altconsensus[pos] = sorted_candidates[1][0][indx]
    if not quiet:
        print('created alt consensus')
    ambiguous = []

    for i in range(2 + len(referenceData[1]), len(extendedReference)):
        sortedfreq = heapq.nlargest(
            2, enumerate(mult_frequencies[i].values()), key=lambda x: x[1]
        )
        sortedfreq2 = heapq.nlargest(
            2, enumerate(frequencies[i].values()), key=lambda x: x[1]
        )
        largestVal = sortedfreq[1][1]

        if largestVal > MIN_SCORE2 and sortedfreq2[1][1] > sum(frequencies[i].values())/3:
            ambiguous.append(i)
    ambiguous.sort()
    ambig_arr = np.array(ambiguous, dtype=np.int32)
    ambig_size = len(ambiguous)

    scores = {}

    for ii in range(len(alignment_fields)):
        fields = alignment_fields[ii]
        pos = alignPos[ii]

        read = fields[9]
        c_string = get_possible_consensus(fields[9], pos, ambig_arr, ambig_size)
        candidate = unicode(c_string)
        free(c_string)
        if candidate in scores:
            scores[candidate] += 1
        else:
            scores[candidate] = 1
    sorted_candidates = sorted(scores.items(), key=lambda k: k[1], reverse=True)
    sorted_candidates = list(
        filter(lambda x: len(x[0].replace("N", "")) > 0, sorted_candidates)
    )
    if len(sorted_candidates) > 1:

        for indx, pos in enumerate(ambiguous):

            altconsensus[pos] = sorted_candidates[1][0][indx]

    altconsensus = "".join(altconsensus).strip("AN")
    mainconsensus = "".join(mainconsensus).strip("AN")
    if altconsensus != mainconsensus:
        if num_branches == BRANCH_LIMIT:
            pass
        else:
            print("Create branch ")
            with open(os.path.splitext(outFile)[0] + "_alt.fa", "w") as out:
                out.writelines([referenceData[0] + "\n" + altconsensus.strip("N")])
            stack.append(
                (
                    os.path.splitext(outFile)[0] + "_alt",
                    args.reads,
                    curPath + "/" + os.path.splitext(os.path.basename(outFile))[0],
                    num_branches + 1,
                )
            )

    with open(outFile, "w") as out:
        out.writelines([referenceData[0] + "\n" + mainconsensus])
    if check_repeats(mainconsensus, STOP_LENGTH):
        print("Repeat found, stopping")
        return -1
    else:
        return len(mainconsensus)


def regenerate_consensus(contig, inputReads, outFile):
    return iterate(contig, inputReads, outFile, exclude=False, analyze_mode=True)


def _main(args):
    global dirpath
    global stack
    global maxlength
    global MIN_SCORE
    global MIN_SCORE2
    global MIN_OVERLAP
    global BRANCH_LIMIT
    global STOP_LENGTH
    global NUM_THREADS
    global COMPLEX_THRESHOLD
    global PAIRED

    MIN_OVERLAP = args.min_overlap_length
    MIN_SCORE2 = args.min_branch_score
    BRANCH_LIMIT = args.branch_limit
    STOP_LENGTH = args.stop_length
    NUM_THRADS = args.threads
    COMPLEX_THRESHOLD = args.complex_threshold
    PAIRED = args.paired

    print("Using bowtie2 at " + get_bowtie_path())
    if PAIRED:
        print("Running in paired mode")
    output_dir = Path(args.reads).stem + "_" + Path(args.reference).stem
    if args.out is not None:
        output_dir = os.path.join(args.out, output_dir)

    dirpath = tempfile.mkdtemp(dir='')
    print("Temp dir: " + dirpath)
    try:
        shutil.rmtree("backup")
    except:
        pass

    try:
        shutil.copytree(output_dir, "backup")
    except:
        pass

    try:
        shutil.rmtree(output_dir)
    except:
        pass
    os.makedirs(output_dir)

    maxlength = 0
    _, filtered_reads = tempfile.mkstemp()

    if COMPLEX_THRESHOLD != -1:
        try:
            # Look for prinseq library in either pyinstaller package or script location
            p_path = os.path.join(sys._MEIPASS if hasattr(sys, '_MEIPASS') else pathlib.Path().absolute(), 'prinseq-lite.pl')
            if not os.path.isfile(p_path):
                raise RuntimeError("Prinseq library not found")
            filter_input = open_file(args.reads).read().encode()

            res = subprocess.run([shutil.which('perl'), p_path, '-fastq', 'stdin', '-lc_method', 'dust', '-lc_threshold', str(COMPLEX_THRESHOLD), '-out_good', filtered_reads, '-out_bad', 'null'], cwd=os.getcwd(), input=filter_input).args
            # print([p_path, '-fastq', os.path.abspath(args.reads), '-lc_method', 'dust', '-lc_threshold', str(COMPLEX_THRESHOLD), '-out_good', filtered_reads, '-out_bad', 'null'])
            filtered_reads = filtered_reads + '.fastq'
            if not os.path.isfile(filtered_reads):
                print("No prinseq output found. Make sure file EOL sequence is correct")
                raise RuntimeError("No prinseq output found")

            args.reads = filtered_reads
        except Exception as e:
            print(e)
            print("prinseq-lite.pl not found, not filtering")
    with open_file(args.reads) as reads:
        lines = []
        for line in reads:
            lines.append(line.rstrip())
            if len(lines) == 4:
                maxlength = max(maxlength, len(lines[1]))
                lines = []
        reads.seek(0)
        readData = reads.read().encode()

    MIN_SCORE = int(pow(10, -args.extend_tolerance) * maxlength * maxlength * args.coverage)
    print("Using extend threshold " + str(MIN_SCORE))
    matched_reads = set()

    fasta = list(SeqIO.parse(open_file(args.reference), 'fasta'))

    if len(fasta) != 1:
        raise ValueError("this function accepts only one contig per FASTA input")
    shutil.copyfile(args.reference, output_dir + "/0.fa")
    os.makedirs(output_dir + "/contigs")
    prev_len = 0

    stack.append((output_dir + "/0", args.reads, "orig", 1))
    while len(stack) > 0:
        top = stack.pop()
        if top[2] != "":
            os.makedirs(output_dir + "/" + top[2])
        inFile = top[0] + ".fa"
        prev_len = 0
        print("Extending " + inFile)
        for i in range(5000):
            start = time.perf_counter()
            res = iterate(
                inFile,
                readData,
                output_dir + "/" + top[2] + "/" + str(i + 1) + ".fa",
                curPath=top[2],
                allow_alt=(i != 0),
                used_reads=matched_reads,
                num_branches=top[3],
                args=args
            )
            inFile = output_dir + "/" + top[2] + "/" + str(i + 1) + ".fa"
            if res!=-1: print("Iteration " + str(i + 1) + " in " + "{0:.4f}".format(time.perf_counter() - start) + "sec, length " + str(res))
            if prev_len == res or res == -1 or res > 100000:
                if res > 100000:
                    print("Length limit exceeded")
                analyze = regenerate_consensus(inFile, readData, output_dir + '/' +
                                               top[2] + '/consensus_temp.fa')

                shutil.copyfile(inFile, output_dir + '/' +
                                top[2] + '/consensus.fa')
                # print(top[2])
                shutil.copyfile(inFile, output_dir + '/contigs' + '/' + top[2].replace('/', '-') + '_coverage_' + str(float(analyze)/float(max(max(1, res), prev_len))) + '.fa')
                break
            prev_len = res

        # shutil.rmtree(dirpath)
    shutil.rmtree(dirpath)
    return output_dir


if __name__ == "__main__":
    run_script(sys.argv[1:])
    if not cython.compiled:
        print("Warning: not running in Cython. Processing time will be affected.")


def run_script(args):
    return _main(args)
