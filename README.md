# contig-extender
## Description
ContigExtender, was developed to extend contigs, complementing de novo assembly. ContigExtender employs a novel recursive Overlap Layout Candidates (r-OLC) strategy that explores multiple extending paths to achieve longer and highly accurate contigs. ContigExtender is effective for extending contigs significantly in in silico synthesized and real metagenomics datasets.

Binaries for 64-bit Linux and Windows are provided in the ```dist``` folder.
![extension process](https://i.imgur.com/w4QiDIj.png "extension process")
## Dependencies
### Runtime:
* Bowtie2 2.3.4+
* Perl

### Development: 
* Python 3.2+ (with development headers)
* numpy
* BioPython
* setuptools
* pyinstaller
* Cython
* gcc
* pytest (for testing)

## Building
### Linux
```
git clone https://github.com/dengzac/contig-extender.git
cd contig-extender
./build.sh
```
Executable is ```dist/extender_wrapper```

### Windows

Tested  using Windows10 with WLS (Windows Linux Subsystem). FYI 
https://docs.microsoft.com/en-us/windows/wsl/install-win10
. Gzipped input fastq file may not work on some Windows machine. We recommend using unzipped fastq on Windows10/WLS.
Native Windows distribution is no-longer supported.

### Tests
Unit and integration tests can be run at `tests/test.sh`
## Usage
For best results, preprocess fastq by removing adaptors before attemping extension
Filtering based on complexity is an available option (uses DUST method from PRINSEQ)
```
usage: extender_wrapper [-h] [--m1 [M1]] [--m2 [M2]] [--out [OUT]]
                           [--enable-pair-constraint]
                           [--min-overlap-length [MIN_OVERLAP_LENGTH]]
                           [--extend-tolerance [EXTEND_TOLERANCE]]
                           [--coverage [COVERAGE]]
                           [--min-branch-score [MIN_BRANCH_SCORE]]
                           [--branch-limit [BRANCH_LIMIT]]
                           [--stop-length [STOP_LENGTH]] [--threads [THREADS]]
                           [--complex-threshold [COMPLEX_THRESHOLD]]
                           [--maxins [MAXINS]]
                           reference [reads]

positional arguments:
  reference             fasta file with contigs to extend
  reads                 fastq file with unpaired reads to extend with; use the
                        --m1 and --m2 options instead for paired reads
                        (default: None)

optional arguments:
  -h, --help            show this help message and exit
  --m1 [M1]             fastq file with #1 mates of paired data (default:
                        None)
  --m2 [M2]             fastq file with #2 mates of paried data (default:
                        None)
  --out [OUT]           output directory (default: output)
  --enable-pair-constraint
                        Require paired-end alignments to satisfy orientation
                        and insert size constraints according to Bowtie2.
                        Otherwise, paired data will be treated as two sets of
                        unpaired reads. (default: False)
  --min-overlap-length [MIN_OVERLAP_LENGTH]
                        minimum length of overlap between candidate read and
                        the existing contig (default: 15)
  --extend-tolerance [EXTEND_TOLERANCE]
                        This parameter, along with read length and coverage,
                        is used to determine the required score to extend.
                        Lower numbers require better quality alignments to
                        extend (default: 2.5)
  --coverage [COVERAGE]
                        estimate of sequencing coverage (default: 10)
  --min-branch-score [MIN_BRANCH_SCORE]
                        A new branch for an alternative contig will be created
                        if it meets this threshold. Then, both the main and
                        alternative branches will be processed recursively.
                        (default: 5)
  --branch-limit [BRANCH_LIMIT]
                        This limits the number of alternative contigs that are
                        considered and output. Once the search tree has this
                        many leaf nodes, no new branches will be created
                        (default: 1)
  --stop-length [STOP_LENGTH]
                        Terminate extension if any substring of this size is
                        repeated within the contig. This prevents circular
                        genomes from infinite extension. (default: 250)
  --threads [THREADS]   Number of threads to use in computing alignments
                        (default: 8)
  --complex-threshold [COMPLEX_THRESHOLD]
                        [0-100] This parameter is passed to PRINSEQ's DUST
                        complexity filter, is run on the reads. Higher values
                        indicate less complexity. -1 to disable (default: 15)
  --maxins [MAXINS]     When using paired-end constraints, this is the maximum
                        fragment length, passed directly to Bowtie2. This
                        length includes both reads ond the gap between them.
                        Refer to the Bowtie2 manual for more information
                        (default: 500)
```

## Examples
The ```examples``` folder contains a simulated dataset from the BKV genome, with a set of reads and a seed contig. To extend, run the following command:

```
./dist/extender_wrapper --coverage 50 examples/BKV_seed_1000_867.fa examples/BKV_250_50_0.01_0_.fastq
```
The --complex-threshold -1 option ignore prinseq quality checking for faster execution.

```
./dist/extender_wrapper --complex-threshold -1 --coverage 50 examples/BKV_seed_1000_867.fa examples/BKV_250_50_0.01_0_.fastq
```
The output contig(s) for each input will be found in the ```output/contigs.fasta```, sorted by length and in the order they appeared in the input. For example, the output contigs for an input sequence ```>reference``` would be named ```>reference_1 reference```, ```>reference_2 reference```, etc. To verify the accuracy of the extended contig, the reference genome is provided in ```BKV.fasta```.

Paired-end example:
```
./dist/extender_wrapper examples/vir_ref.fa --m1 examples/reads_1.fq --m2 examples/reads_2.fq --enable-pair-constraint
```
Paired alignments with incorrect orientation or fragment lengths longer than 500bp are excluded.

use higher value for --extend-tolerance [2.5, 3.5] (default 2.5) to force extension,  but risks higher chances of false extension.
```
./dist/extender_wrapper  --extend-tolerance 3.0  examples/BKV_seed_1000_867.fa examples/BKV_250_50_0.01_0_.fastq
```


