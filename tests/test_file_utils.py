"""
Tests functions in file_utils.py
"""

from src.io_utils import s3_check
from src.file_utils import read_handle, gzip_check

# -------------------------
# Definitions
# -------------------------


TEST_BLANK = 'tests/data/consensus_genome/blank.fastq.gz'
TEST_NO_HOST = 'tests/data/consensus_genome/no_host_1.fq.gz'
TEST_CT20K = 'tests/data/consensus_genome/Ct20K.fastq.gz'
TEST_LOCAL_TEXT = 'tests/data/io/hello.txt'
TEST_LOCAL_ABSENT = 'tests/data/io/idontexist.txt'
TEST_S3 = 's3://cypherid-public-references/hi.txt'

def test_gzip_check():

    s3, bucket_name, path = s3_check(TEST_NO_HOST)
    assert gzip_check(s3, bucket_name, path)

    s3, bucket_name, path = s3_check(TEST_LOCAL_TEXT)
    assert not gzip_check(s3, bucket_name, path)

    s3, bucket_name, path = s3_check(TEST_S3)
    assert not gzip_check(s3, bucket_name, path)

def test_read_handle():

    handle = read_handle(TEST_NO_HOST)
    assert handle is not None
    handle.close()

    handle = read_handle(TEST_LOCAL_TEXT)
    assert handle is not None
    handle.close()

    handle = read_handle(TEST_NO_HOST)
    assert handle is not None
    handle.close()

    handle = read_handle(TEST_S3)
    assert handle is not None
    handle.close()