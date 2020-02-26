import extender
import hashlib
import os
import extender_wrapper
def test_replace_string():
    orig = "AOEURCHEU"
    replaced = extender.replace_string(orig, 'N', 1, 3)
    assert replaced == "ANNNRCHEU"

def test_open_gzip():
    file = extender.open_file("BKV_250_50_0.01_0_.fastq.gz")
    content = file.read()
    assert hashlib.md5(content.encode('utf-8')
                       ).hexdigest() == "b10fe0be40303f7fd93370d9bd5336da"

def test_bowtie_path():
    path = extender.get_bowtie_path()
    build = extender.get_bowtie_build_path()
    assert os.path.dirname(path) == os.path.dirname(build)


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def test_extend():
    extender._main(Namespace(branch_limit=1, complex_threshold=-1, coverage=50.0, extend_tolerance=2.5, m1=None, m2=None, maxins=500, min_branch_score=5,
                                  min_overlap_length=15, out='output', pair_constraint=False, paired=False, reads='BKV_250_50_0.01_0_.fastq', reference='BKV_seed_1000_867.fa', stop_length=250, threads=8))
    outfilename = "output/BKV_250_50_0.01_0__BKV_seed_1000_867/contigs/orig_coverage_49.06945510360706.fa"
    assert os.path.isfile(outfilename)
    file = open(outfilename)
    assert file.readline().strip()==">867"

def test_complex_threshold():
    extender._main(Namespace(branch_limit=1, complex_threshold=1, coverage=50.0, extend_tolerance=2.5, m1=None, m2=None, maxins=500, min_branch_score=5,
                             min_overlap_length=15, out='output', pair_constraint=False, paired=False, reads='BKV_250_50_0.01_0_.fastq', reference='BKV_seed_1000_867.fa', stop_length=250, threads=8))
    outfilename = "output/BKV_250_50_0.01_0__BKV_seed_1000_867/contigs/orig_coverage_24.535123966942148.fa"
    assert os.path.isfile(outfilename)
    file = open(outfilename)
    assert file.readline().strip() == ">867"

def test_score():
    assert extender.calc_min_score(1.5, 30, 500) == 237170


def test_wrapper_unpair():
    args = extender_wrapper.parse_args(['--complex-threshold', '-1', '--coverage',
                                        '50', 'examples/BKV_seed_1000_867.fa', 'examples/BKV_250_50_0.01_0_.fastq'])
    assert args.complex_threshold == -1
    assert args.coverage == 50
    assert args.reference == "examples/BKV_seed_1000_867.fa"
    assert args.reads == "examples/BKV_250_50_0.01_0_.fastq"


def test_wrapper_pair():
    args = extender_wrapper.parse_args(
        ['examples/vir_ref.fa', '--m1', 'examples/reads_1.fq', '--m2', 'examples/reads_2.fq', '--enable-pair-constraint'])

    assert args.reference == "examples/vir_ref.fa"
    assert args.m1 == "examples/reads_1.fq"
    assert args.m2 == "examples/reads_2.fq"
    assert args.pair_constraint
    
