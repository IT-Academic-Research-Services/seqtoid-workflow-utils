"""

"""


from src.fastaq import fastq_iterate
from src.file_utils import read_handle

TEST_BLANK = 'data/consensus-genome/blank.fastq.gz'
TEST_NO_HOST = 'data/consensus-genome/no_host_1.fq.gz'


def test_fastq_iterate():

    fh = read_handle(TEST_NO_HOST)
    lc = 0
    for _ in fastq_iterate(fh):
        lc += 1
    assert lc == 3
