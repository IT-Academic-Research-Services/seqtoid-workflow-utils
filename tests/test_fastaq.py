"""

"""


from src.fastaq import fastq_iterate
from src.file_utils import read_handle

TEST_BLANK = 'tests/data/consensus_genome/blank.fastq.gz'
TEST_NO_HOST = 'tests/data/consensus_genome/no_host_1.fq.gz'
TEST_CT20K = 'tests/data/consensus_genome/Ct20K.fastq.gz'


def test_fastq_iterate():

    fh = read_handle(TEST_NO_HOST)
    lc = 0
    for _ in fastq_iterate(fh):
        lc += 1
    assert lc == 3


    fh = read_handle(TEST_BLANK)
    lc = 0
    for _ in fastq_iterate(fh):
        lc += 1
    assert lc == 0

    fh = read_handle(TEST_CT20K)
    lc = 0
    for header, seq, qual in fastq_iterate(fh):
        if lc == 0:
            assert header == '@d5115431-c39a-46f9-8c6c-94a1853842dc'
        lc += 1
    assert lc == 912

