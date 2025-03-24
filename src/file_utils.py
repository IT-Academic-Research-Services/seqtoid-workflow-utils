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
    :return: Boolean indicating if the file is gzipped
    """
    if not isinstance(path, str) or not path.strip():
        get_logger().error("Path must be a non-empty string, got: %s", path)
        return None

    if s3 is None:  # Local file
        try:
            with open(path, "rb") as f:
                return f.read(2) == GZIP_MAGIC_NUMBER
        except FileNotFoundError:
            get_logger().error(f"Local file not found: {path}")
            return None
        except PermissionError:
            get_logger().error(f"Permission denied for file: {path}")
            return None
        except OSError as e:
            get_logger().error(f"Error reading local file: {e}")
            return None
    else:  # S3 file
        if not isinstance(bucket_name, str) or not bucket_name.strip():
            get_logger().error("Bucket name must be a non-empty string, got: %s", bucket_name)
            return None

        try:
            obj = s3.get_object(Bucket=bucket_name, Key=path, Range='bytes=0-1')
            res = obj['Body'].read()
            if len(res) < 2:
                get_logger().warning("S3 object too short to check gzip magic number: %s bytes", len(res))
                return False  # Too short to be gzipped
            return res == GZIP_MAGIC_NUMBER
        except s3.exceptions.NoSuchKey:
            get_logger().error(f"S3 object not found: s3://{bucket_name}/{path}")
            return None
        except Exception as e:
            get_logger().error(f"Error retrieving S3 object: {e}")
            return None


def read_handle(file_path, profile_name=None, force_gz=False, encoding='utf-8'):
    """
    Returns a read handle. Works with S3 or local files.
    :param file_path: If S3, of form s3://bucket-name/dir/item-name.txt. If local, full path to file.
    :param profile_name: AWS profile name for S3 access (optional).
    :param force_gz: If True, treat the file as gzipped regardless of detection.
    :param encoding: Text encoding for reading (default: 'utf-8').
    :return: Handle for reading file. Caller must close the handle (e.g., using `with` statement).
    """
    if not isinstance(file_path, str) or not file_path.strip():
        raise ValueError("file_path must be a non-empty string")

    try:
        s3, bucket_name, path = s3_check(file_path, profile_name=profile_name)
        gzipped = gzip_check(s3, bucket_name, path)
    except ValueError as e:
        get_logger().error(f"Invalid file path or configuration: {e}")
    except Exception as e:
        get_logger().error(f"Error checking file: {e}")

    else:
        if force_gz:
            gzipped = True

        if s3:
            if not path:
                get_logger().error(f"Invalid S3 path: no key specified in {file_path}")
            try:
                obj = s3.get_object(Bucket=bucket_name, Key=path)
            except s3.exceptions.NoSuchKey:
                get_logger().error(f"S3 object not found: s3://{bucket_name}/{path}")
            except Exception as e:
                get_logger().error(f"Error retrieving S3 object: {e}")

            if gzipped:
                try:
                    res = gzip.GzipFile(None, 'rb', fileobj=obj['Body'])
                except OSError as e:
                    get_logger().error(f"Failed to decompress S3 object as gzip: {e}")
            else:
                res = obj['Body']

            try:
                return io.TextIOWrapper(res, encoding=encoding)
            except UnicodeDecodeError as e:
                get_logger().error(f"File cannot be decoded with encoding {encoding}: {e}")
        else:
            try:
                if gzipped:
                    return gzip.open(file_path, 'rt', encoding=encoding)
                else:
                    return open(file_path, 'r', encoding=encoding)
            except FileNotFoundError:
                get_logger().error(f"Local file not found: {file_path}")
            except PermissionError:
                get_logger().error(f"Permission denied for file: {file_path}")
            except OSError as e:
                get_logger().error(f"Failed to open local file as gzip: {e}")