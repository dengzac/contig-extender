# cython: profile=True
# cython: linetrace=True
# cython: binding=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
import time
import argparse
import tempfile
import shutil
import subprocess
import os
import re
import numpy as np
from scipy.special import softmax
import heapq
import psutil
import sys
from pathlib import Path
import multiprocessing
import pyximport
import cython
import array
cimport numpy as np
import itertools
from libc.stdlib cimport malloc, free
#pyximport.install()
cdef int get_char_pos(char c):
    #print(c)
    if c == 'A':
        return 0
    elif c == 'C':
        return 1
    elif c == 'G':
        return 2
    elif c == 'T':
        return 3

cdef int MIN_SCORE = 10000
cdef int MIN_SCORE2 = 5
cdef int MIN_OVERLAP = 15
cdef int BRANCH_LIMIT = 1
cdef int STOP_LENGTH = 250
cdef int NUM_THREADS = multiprocessing.cpu_count()
cdef int COMPLEX_THRESHOLD = 15
cdef int maxlength
dirpath = ''
def replace_string(string, replace_letter, startPos, endPos):
    return (
        string[0:startPos]
        + (replace_letter * (endPos - startPos + 1))
        + string[endPos + 1 :]
    )

cdef int c_max(int a, int b, int c):
    cdef int best = 0
    if a > best: best = a
    if b > best: best = b
    if c > best: best = c
    if a==best: return a
    if b==best: return b
    if c==best: return c
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
    analyze_mode = False
):

    cdef int idx

    with open(contig) as reference:
        referenceData = reference.readlines()
        referenceData = [x.strip() for x in referenceData]
        referenceData[1] = referenceData[1].strip("N")

    extendedReference = (
        "A" + ("N" * (maxlength + 1)) + referenceData[1] + ("N" * (maxlength + 1)) + "A"
    )
    origReference = extendedReference
    if exclude:
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
        ["bowtie2-build", extendedReferenceFile, dirpath + "/ecoli"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if not quiet:
        print("Finding matches")
    # samOutput = dirpath+"/out.sam"
    bamOutput = dirpath + "/out2.sam"
    #print(NUM_THREADS)
    result = subprocess.run(
        [
            "bowtie2",
            "-x",
            dirpath + "/ecoli",
            "-U",
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
        input=readData.encode(),
        stderr=subprocess.DEVNULL,
    )

    if not quiet:
        print("Building consensus")
    frequencies = np.ones((len(extendedReference), 4), dtype=np.int64)
    mult_frequencies=np.zeros((len(extendedReference), 4), dtype=np.int64)
    cdef long[:,:] freq_view = frequencies
    cdef long[:, :] mult_view = mult_frequencies
    # frequencies = [
    #     {"A": 1, "C": 1, "G": 1, "T": 1} for i in range(len(extendedReference))
    # ]
    # mult_frequencies = [
    #     {"A": 0, "C": 0, "G": 0, "T": 0} for i in range(len(extendedReference))
    # ]

    mainconsensus = ""
    altconsensus = ""
    alignments = result.stdout.decode().strip().split("\n")
    print(len(alignments), 'alignments')
    splitter = re.compile(r"\t+")
    if not quiet:
        print('filter')
    alignments = filter(lambda a: a[0]!="@" and a[3] != "1", alignments)
    if not quiet:
        print('split')
    alignment_fields = np.array([np.array([x for x in re.split(splitter, o, maxsplit=10)]) for o in alignments], dtype=np.object_)
    if not quiet:
        print('pos')
    cdef int[:] alignPos = alignment_fields[:, 3].astype(np.int32)-1
    #np.apply_along_axis(lambda a: bytearray(a.encode()), 9, alignment_fields)
    #print(alignment_fields)
    cdef np.ndarray[dtype=object, ndim=2] align_arr = alignment_fields
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
        for ii in range(len(alignment_fields)):#fields in alignment_fields:
            fields = alignment_fields[ii]
            rLen = len(fields[9])
            total_coverage += rLen
            #fields = re.split(splitter, line)
            pos = alignPos[ii]#int(fields[3]) - 1

            # if not quiet:
            #     print(pos)
            read = fields[9]
            r = fields[9]#read
            quality1 = 0
            quality2 = 0

            if pos < (maxlength + 2):
                quality1 = rLen - (maxlength + 2) + pos
            elif pos + rLen > maxlength + 2 + refLength:
                quality2 = maxlength + 2 + refLength - pos
            else:
                for idx in range(rLen):
                    
                    if r[idx] == ext[idx + pos]:
                        quality2 += 1
            quality = c_max(quality1, quality2, MIN_OVERLAP)
            if quality == MIN_OVERLAP:
                quality = 0

            for idx in range(rLen):
                if r[idx] != n_char: 
                    c_pos = idx + pos


                    char_pos = get_char_pos(r[idx])
                    freq_view[c_pos, char_pos] += quality * quality
                    mult_view[c_pos, char_pos] += 1 #quality
    if analyze_mode: return total_coverage
    keys = ["A", "C", "G", "T"]
    frequencies = [{"A": x[0], "C": x[1], "G": x[2], "T": x[3]}
                for x in frequencies]
    mult_frequencies = [{"A": x[0], "C": x[1], "G": x[2], "T": x[3]}
                for x in mult_frequencies]
    #print (frequencies)
    if not quiet: print("Building main consensus")
    mainconsensus = list(origReference)
    #2 + 2 * maxlength + 2*maxlength,
    
    for i in itertools.chain(range(2 + 2 * maxlength), range(len(extendedReference) - 2 - 2 * maxlength, len(extendedReference))):
        largestVal = sum(frequencies[i].values())
        if largestVal > MIN_SCORE:
            if max(frequencies[i], key=lambda k: frequencies[i][k]) != mainconsensus[i] and mainconsensus[i] != 'N':
                if not quiet: print('replace', i, frequencies[i])
            mainconsensus[i] = max(frequencies[i], key=lambda k: frequencies[i][k])

            #print('skip', i)

            # if i >= maxlength+2 and i < maxlength + 2 + len(referenceData[1]):
            #     mainconsensus[i] = extendedReferenc[i]
            # else:
            #     mainconsensus[i] = "N"
    #print(extendedReference)
    #print(''.join(mainconsensus))
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
        # if largestVal > 1:
        #     print(i, largestVal)
        if largestVal > MIN_SCORE2 and sortedfreq2[1][1] > sum(frequencies[i].values())/3:
            ambiguous.append(i)
    ambiguous.sort()
    cdef int[:] ambig_arr = np.array(ambiguous, dtype=np.int32)
    cdef int ambig_size = len(ambiguous)
    scores = {}
    if not quiet:
        print('finding candidates')
    for ii in range(len(alignment_fields)):#fields in alignment_fields:
        fields = alignment_fields[ii]
        pos = alignPos[ii]#int(fields[3]) - 1

        # if not quiet:
        #     print(pos)
        
        c_string = get_possible_consensus(fields[9], pos, ambig_arr, ambig_size)
        candidate = unicode(c_string)
        free(c_string)
        
        if candidate in scores:
            scores[candidate] += 1
        else:
            scores[candidate] = 1
    # print(scores)
    # print(scores.items())
    if not quiet:
        print('altconsensus part 1')
    sorted_candidates = sorted(scores.items(), key=lambda k: k[1])
    sorted_candidates = list(
        filter(lambda x: len(x[0].replace("N", "")) > 0, sorted_candidates)
    )
    #print(ambiguous)
    if len(sorted_candidates) > 1:
        for indx, pos in enumerate(ambiguous):
            #print(pos, idx)
            altconsensus[pos] = sorted_candidates[1][0][indx]
    if not quiet:
        print('created alt consensus')
    ambiguous = []
    #print(extendedReference)
    for i in range(2 + len(referenceData[1]), len(extendedReference)):
        sortedfreq = heapq.nlargest(
            2, enumerate(mult_frequencies[i].values()), key=lambda x: x[1]
        )
        sortedfreq2 = heapq.nlargest(
            2, enumerate(frequencies[i].values()), key=lambda x: x[1]
        )
        largestVal = sortedfreq[1][1]
        #print(i, sortedfreq)
        # if largestVal >= 1 :
        #     print(i, largestVal)
        if largestVal > MIN_SCORE2 and sortedfreq2[1][1] > sum(frequencies[i].values())/3:
            ambiguous.append(i)
    ambiguous.sort()
    ambig_arr = np.array(ambiguous, dtype=np.int32)
    ambig_size = len(ambiguous)
    #print(ambiguous)
    scores = {}
    #corresponding_reads = {}
    for ii in range(len(alignment_fields)):#fields in alignment_fields:
        fields = alignment_fields[ii]
        pos = alignPos[ii]#int(fields[3]) - 1

        # if not quiet:
        #     print(pos)
        read = fields[9]
        c_string = get_possible_consensus(fields[9], pos, ambig_arr, ambig_size)
        candidate = unicode(c_string)
        free(c_string)
        if candidate in scores:
            scores[candidate] += 1
            #corresponding_reads[candidate].add(fields[0])
        else:
            scores[candidate] = 1
            #corresponding_reads[candidate] = set([fields[0]])
    sorted_candidates = sorted(scores.items(), key=lambda k: k[1], reverse=True)
    sorted_candidates = list(
        filter(lambda x: len(x[0].replace("N", "")) > 0, sorted_candidates)
    )
    if len(sorted_candidates) > 1:
        #print(ambiguous)
        #print(sorted_candidates)
        for indx, pos in enumerate(ambiguous):
            #print(indx, pos)
            altconsensus[pos] = sorted_candidates[1][0][indx]

    altconsensus = "".join(altconsensus).strip("AN")
    mainconsensus = "".join(mainconsensus).strip("AN")
    if altconsensus != mainconsensus:  # branch and allow_alt:
        # print("Alternate detected ")  # + altconsensus)
        if num_branches == BRANCH_LIMIT:
            # print("Reached limit, skipping")
            pass
        else:
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
    if check_repeats(mainconsensus):
        print("Repeat found, stopping")
        return -1
    else:
        return len(mainconsensus)


def check_repeats(contig):
    return contig[-STOP_LENGTH:] in contig[:-STOP_LENGTH]


cdef char* get_possible_consensus(char[:] read, int pos, int[:] ambigPos, int ambigN):

    cdef char* c_string = <char*>malloc((ambigN+1)*sizeof(char))
    s = ""
    cdef int[:] ambig_arr
    cdef int i, readIndex

    ambig_arr = ambigPos#np.array(ambigPos, dtype=np.int32)

    cdef int n = len(read)
    cdef int idx
    for idx in range(ambigN):
        i = ambig_arr[idx]
        readIndex = i - pos
        if readIndex < 0 or readIndex >= n:
            c_string[idx] = 'N'#s += "N"
        else:
            c_string[idx] = read[readIndex]#s += r[readIndex]
    return c_string
    


def regenerate_consensus(contig, inputReads, outFile):
    return iterate(contig, inputReads, outFile, exclude=False, analyze_mode=True)


def prepare_directory(name):
    try:
        shutil.rmtree(name)
    except:
        pass
    os.makedirs(name)

stack = []
import cProfile
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("reference", help="fasta file with contig to extend")
    parser.add_argument("reads", help="fastq file with reads to extend with")

    parser.add_argument("out", nargs="?", help="output directory")
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
        help="[0-100] higher values indicate less complexity",
        nargs="?",
        default=COMPLEX_THRESHOLD,
        type=int
    )
    args = parser.parse_args()

    MIN_OVERLAP = args.min_overlap_length
    #MIN_SCORE = args.extend_tolerance*#args.min_extend_score
    MIN_SCORE2 = args.min_branch_score
    BRANCH_LIMIT = args.branch_limit
    STOP_LENGTH = args.stop_length
    NUM_THRADS = args.threads
    COMPLEX_THRESHOLD = args.complex_threshold

    output_dir = "output"
    if args.out is not None:
        output_dir = args.out
    output_dir = Path(args.reads).stem + "_" + Path(args.reference).stem
    print(output_dir)
    if not shutil.which("bowtie2"):
        print("Error: bowtie2 not found. Add the executable location to PATH")
        sys.exit()
    dirpath =  tempfile.mkdtemp(dir='')
    print("Temp dir", dirpath)
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
    # print(subprocess.run(['./prinseq-lite.pl', '-fastq', args.reads, '-lc_method', 'dust', '-lc_threshold', str(COMPLEX_THRESHOLD), '-out_good', filtered_reads, '-out_bad', 'null']).args)
    filtered_reads = filtered_reads + '.fastq'
    # print(filtered_reads)
    # args.reads = filtered_reads
    with open(args.reads) as reads:
        lines = []
        for line in reads:
            lines.append(line.rstrip())
            if len(lines) == 4:
                maxlength = max(maxlength, len(lines[1]))
                lines = []
        reads.seek(0)
        readData = reads.read()
        # print(readData)
        # print(reads)
    MIN_SCORE = int(pow(10, -args.extend_tolerance) * maxlength * maxlength * args.coverage)
    print("Using extend threshold ", MIN_SCORE)
    matched_reads = set()
    #iterate(args.reference, args.reads, output_dir + "/0.fa")
    shutil.copyfile(args.reference, output_dir + "/0.fa")
    os.makedirs(output_dir + "/contigs")
    prev_len = 0
    stack = []
    stack.append((output_dir + "/0", args.reads, "orig", 1))
    while len(stack) > 0:
        top = stack.pop()
        if top[2] != "":
            os.makedirs(output_dir + "/" + top[2])
        inFile = top[0] + ".fa"
        prev_len = 0
        for i in range(5000):
            start = time.perf_counter()
            print(
                "Read "
                + inFile
                + " output "
                + output_dir + "/"
                + top[2]
                + "/"
                + str(i + 1)
                + ".fa"
            )
            res = iterate(
                inFile,
                args.reads,
                output_dir + "/" + top[2] + "/" + str(i + 1) + ".fa",
                curPath=top[2],
                allow_alt=(i != 0),
                used_reads=matched_reads,
                num_branches=top[3],
            )
            inFile = output_dir + "/" + top[2] + "/" + str(i + 1) + ".fa"
            print(i + 1, (time.perf_counter() - start), res)
            if prev_len == res or res == -1 or res > 10000:
                if res > 10000:
                    print("Length limit exceeded")
                analyze = regenerate_consensus(inFile, args.reads, output_dir + '/' +
                                top[2] + '/consensus_temp.fa')
                print("No progress, finished")
                shutil.copyfile(inFile, output_dir + '/' +
                                top[2] + '/consensus.fa')
                shutil.copyfile(inFile, output_dir + '/contigs' + '/' + top[2].replace('/', '') + '_coverage_' + str(float(analyze)/float(res)) + '.fa')
                break
            prev_len = res

        # shutil.rmtree(dirpath)
    shutil.rmtree(dirpath)
