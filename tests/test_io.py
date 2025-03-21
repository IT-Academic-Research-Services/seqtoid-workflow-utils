"""
This module contains tests for the io_utils module.
"""

import os
from src.io_utils import s3_check, file_check, get_file

# -------------------------
# Definitions
# -------------------------

TEST_S3 = 's3://cypherid-public-references/hi.txt'
TEST_S3_ABSENT = 's3://cypherid-public-references/idontexist.txt'
TEST_S3_BAD_BUCKET = 's3://idontexist/hi.txt'
TEST_LOCAL = 'tests/data/io/hello.txt'
TEST_LOCAL_ABSENT = 'tests/data/io/idontexist.txt'

def test_s3_check():

    s3, bucket, path = s3_check(TEST_S3)
    assert s3 is not None

    s3, bucket, path = s3_check(TEST_S3_ABSENT)
    assert s3 is not None

    s3, bucket, path = s3_check(TEST_S3_BAD_BUCKET)
    assert s3 is not None

    s3, bucket, path = s3_check(TEST_LOCAL)
    assert s3 is None

    s3, bucket, path = s3_check(TEST_LOCAL_ABSENT)
    assert s3 is None

def test_file_check():

    s3, bucket, path = file_check(TEST_S3)
    assert path is not None

    s3, bucket, path = file_check(TEST_S3_ABSENT)
    assert path is None

    s3, bucket, path = file_check(TEST_S3_BAD_BUCKET)
    assert path is None

    s3, bucket, path = file_check(TEST_LOCAL)
    assert path is not None

    s3, bucket, path = file_check(TEST_LOCAL_ABSENT)
    assert s3 is None

def test_get_file():

    local_path = get_file(TEST_S3)
    assert os.path.isfile(local_path)
    os.remove(local_path)

    local_path = get_file(TEST_S3_ABSENT)
    assert local_path is None

    local_path = get_file(TEST_S3_BAD_BUCKET)
    assert local_path is None

    local_path = get_file(TEST_LOCAL)
    assert os.path.isfile(local_path)

    local_path = get_file(TEST_LOCAL_ABSENT)
    assert local_path is None

