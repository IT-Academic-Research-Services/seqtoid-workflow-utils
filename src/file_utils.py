"""
Functions for manipulating files and file names.
"""

import os
import io
import gzip

from src.defs import GZIP_EXT_SET, GZIP_MAGIC_NUMBER
from src.io_utils import s3_check
from src.logging_utils import get_logger

# -------------------------
# Definitions
# -------------------------


# -------------------------
# Functions
# -------------------------

def extension_remover_gzip(filename):
    """
    Removes extension, works even if file is <file>.ext.gz.
    :param filename: Name of file.
    :return: Base of file.
    """

    if filename.endswith(tuple(GZIP_EXT_SET)):
        first_base, _ = os.path.splitext(filename)
    else:
        first_base = filename

    final_base, _ = os.path.splitext(first_base)
    return final_base

def gzip_check(s3, bucket_name, path):
    """
    Checks if a file is gzipped. Works on S3 or local.
    :param s3: s3 client or None
    :param bucket_name: s3 bucket or None
    :param path: local path or s3 prefix
    :return: file handle
    """

    if path is None:
        get_logger().error("Path cannot be None")
        raise Exception("Path cannot be None")
    if s3 is None:  # a local file
        with open(path, "rb") as f:
            return f.read(2) == GZIP_MAGIC_NUMBER
    else:
        if bucket_name is None:
            get_logger().error("Bucket name cannot be None")
            raise Exception("Bucket name cannot be None")
        obj = s3.get_object(Bucket=bucket_name, Key=path, Range='bytes={}-{}'.format(0, 2 - 1))
        res = obj['Body'].read()
        return res == GZIP_MAGIC_NUMBER

def read_handle(file_path, profile_name=None, force_gz=False):
    """
    Returns a read handle.
    :param file_path: If S3, of form s3://bucket-name/dir/item-name.txt. If local, full path to file.
    :return: Handle for reading file.
    """

    s3, bucket_name, path = s3_check(file_path, profile_name=profile_name)
    gzipped = gzip_check(s3, bucket_name, path)

    if force_gz:
        gzipped = True

    if s3:
        obj = s3.get_object(Bucket=bucket_name, Key=path)
        if gzipped:
            res = gzip.GzipFile(None, 'rb', fileobj=obj['Body'])
        else:
            res = obj['Body']

        return io.TextIOWrapper(res)
    else:
        if gzipped:
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')
