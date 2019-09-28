# contig-extender
## Description
ContigExtender, was developed to extend contigs, complementing de novo assembly. ContigExtender employs a novel recursive Overlap Layout Candidates (r-OLC) strategy that explores multiple extending paths to achieve longer and highly accurate contigs. ContigExtender is effective for extending contigs significantly in in silico synthesized and real metagenomics datasets.

![extension process](https://i.imgur.com/w4QiDIj.png "extension process")
## Dependencies
###Runtime:
* Bowtie2
* Perl (optional)

###Development: 
* Python 3.2+ (with development headers)
* numpy
* BioPython
* setuptools
* pyinstaller
* Cython
* gcc

## Building
```
git clone https://github.com/dengzac/contig-extender.git
cd contig-extender
./build.sh
```
Executable is ```dist/extender_wrapper```


## Usage
```
usage: extender_wrapper [-h] [--min-overlap-length [MIN_OVERLAP_LENGTH]]
              [--extend-tolerance [EXTEND_TOLERANCE]] [--coverage [COVERAGE]]
              [--min-branch-score [MIN_BRANCH_SCORE]]
              [--branch-limit [BRANCH_LIMIT]] [--stop-length [STOP_LENGTH]]
              [--threads [THREADS]] [--complex-threshold [COMPLEX_THRESHOLD]]
              reference reads [out]

positional arguments:
  reference             fasta file with contigs to extend
  reads                 fastq file with reads to extend with
  out                   output directory (default: None)

optional arguments:
  -h, --help            show this help message and exit
  --min-overlap-length [MIN_OVERLAP_LENGTH]
                        minimum length of overlap between candidate read and
                        contig (default: 15)
  --extend-tolerance [EXTEND_TOLERANCE]
                        lower numbers require more reads to extend (default:
                        2.5)
                        threshold score is proportional to 10^(-tol)
  --coverage [COVERAGE]
                        estimate of coverage (default: 10)
  --min-branch-score [MIN_BRANCH_SCORE]
                        minimum score required to create alternative contig
                        (default: 5)
  --branch-limit [BRANCH_LIMIT]
                        number of alternative contigs to output (default: 1)
  --stop-length [STOP_LENGTH]
                        terminate extension if substring of this size is
                        repeated within contig (default: 250)
  --threads [THREADS]   number of threads to use in computing alignments
                        (default: auto)
  --complex-threshold [COMPLEX_THRESHOLD]
                        [0-100] higher values indicate less complexity. 
                        -1 to disable
                        Uses DUST score from prinseq
                        (default: 15)
```

## Examples
The ```examples``` folder contains a simulated dataset from the BKV genome, with a set of reads and a seed contig. To extend, run the following command:
```
./dist/extender_wrapper --complex-threshold -1 --coverage 50 examples/BKV_seed_1000_867.fa examples/BKV_250_50_0.01_0_.fastq
```
The output contig(s) for each input will be found in the ```output/contigs.fasta```, sorted by length and in the order they appeared in the input. For example, the output contigs for the input sequence ```>reference``` will be named ```>reference_1 reference```, ```>reference_2 reference```, etc. To verify the accuracy of the extended contig, the reference genome is provided in ```BKV.fasta```.
