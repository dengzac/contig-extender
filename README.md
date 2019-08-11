# contig-extender
![extension process](https://i.imgur.com/w4QiDIj.png "extension process")
## Dependencies
* Python 3 (with development headers)
* Cython
* numpy
* setuptools
* gcc

## Building
```
python setup.py build
```
Executable will be named extender/extend

## Usage
```
usage: extend [-h] [--min-overlap-length [MIN_OVERLAP_LENGTH]]
              [--extend-tolerance [EXTEND_TOLERANCE]] [--coverage [COVERAGE]]
              [--min-branch-score [MIN_BRANCH_SCORE]]
              [--branch-limit [BRANCH_LIMIT]] [--stop-length [STOP_LENGTH]]
              [--threads [THREADS]] [--complex-threshold [COMPLEX_THRESHOLD]]
              reference reads [out]

positional arguments:
  reference             fasta file with contig to extend
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
                        [0-100] higher values indicate less complexity. -1 to disable
                        (default: 15)
```
