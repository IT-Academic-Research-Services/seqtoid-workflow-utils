"""
Functions for input/output, including S3.
"""

import os
import logging

import boto3
from botocore.exceptions import ClientError

from src.logging_utils import get_logger

# -------------------------
# Definitions
# -------------------------

S3_TAG = 's3:'
S3_BUCKET_POS = 2

# -------------------------
# Setup
# -------------------------


# -------------------------
# Functions
# -------------------------

def s3_client(profile_name=None):

    session = boto3.Session(profile_name=profile_name)
    return session.client('s3')

def s3_check(in_string, profile_name=None):
    """
    Given a string, determine if it is an S3 path or a valid local path.
    If s3, split into an s3 clint, bucket name, and full prefix.
    If local, return None, None, and the absolute path.
    :param in_string: String, expected to be an S3 path or a local path.
    :param profile_name: Optional AWS profile.
    :return: Tuple (s3_client or None, bucket_name or None, prefix or abs_path)
    """

    if not isinstance(in_string, str) or not in_string.strip():
        get_logger().error('Input string is invalid. Returning None.')
        return None, None, None

    in_string = in_string.replace('\\', '/')  # Normalize Windows paths by replacing backslashes

    try:
        string_cols = in_string.split('/')
    except AttributeError:
        get_logger().error('Cannot split input string. Returning None.')
        return None, None, None

    if in_string.startswith(S3_TAG) and len(string_cols) > S3_BUCKET_POS:

        if not in_string.startswith('s3://'):
            get_logger().error(f'Malformed s3 string {in_string} Appears to be missing double forward slash. Returning None.')
            return None, None, None

        bucket_name = string_cols[S3_BUCKET_POS]
        if not bucket_name or not (3 <= len(bucket_name) <= 63) or not bucket_name.islower():
            get_logger().error(f"Improper bucket name {bucket_name}. Returning None.")
            return None, None, None

        try:
            s3 = s3_client(profile_name)
            s3.head_bucket(Bucket=bucket_name)
        except ClientError:
            get_logger().error(f"Bucket name {bucket_name} fails at s3.head_bucket. Returning None.")
            return None, None, None

        prefix_name = '/'.join(string_cols[S3_BUCKET_POS + 1:]).strip('/')

        if len(prefix_name) < 1:
            get_logger().warning(f"S3 valid but not {prefix_name}. Returning None for prefix.")
            return s3, bucket_name, None
        else:
            return s3, bucket_name, prefix_name

    else:
        try:
            abs_path = os.path.abspath(in_string)
            return None, None, abs_path
        except OSError:
            get_logger().error(f"Unable to resolve absolute path of {in_string}. Returning None.")
            return None, None, None

def file_check(in_string, profile_name=None):
    """
    Check if a file exists, either locally or on S3.
    Designed to call s3_check, and then check if the file exists.
    :param in_string: String, expected to be an S3 path or a local path.
    :param profile_name: Optional AWS profile.
    :return: Tuple (s3_client or None, bucket_name or None, prefix or abs_path)
    """

    s3, bucket_name, path = s3_check(in_string)

    if path is None:
        return None, None, None

    if s3 is None:
        if os.path.isfile(path):
            return None, None, path
        else:
            get_logger().error(f"Local file {path} does not exist.")
            return None, None, None
    else:
        try:
            s3.head_object(Bucket=bucket_name, Key=path)
        except ClientError as e:
            get_logger().error(e)
            return None, None, None
        else:
            return s3, bucket_name, path

def get_file(in_string, profile_name=None, overwrite_local=True):
    """
    The function first checks whether the file string is a valid s3 or local path.
    It will receive either a tuple (s3, bucket, prefix) if s3, (None, None, abs_path)
    if local, or None, None, None if there was an error with the check.

    It will then check that the specified file exists, returning None if it does not.

    If the path is local it will then return the path, which will be the absolute path
    as defined in s3_check.

    If the path is S3, it will define a local path in the cwd, with the base name of the
    path as the file name.

    If a local file by that name already exists, it will log a warning and return the existing path
    if overwrite is False.

    If the file is absent locally and/or overwrite is true, it will attempt to download the file
    from S3 to the local location.

    NB: The overwrite process is simple and pays no attention to datestamps.

    :param in_string: String, expected to be an S3 path or a local path.
    :param profile_name: Optional AWS profile.
    :param overwrite_local: If True, overwrite local file if it exists, only when in_string is an S3 location.

    :return: Absolute local path to file.
    """

    s3, bucket_name, path = file_check(in_string, profile_name=profile_name)

    if path is None:  # Represents error in file_check or s3_check
        get_logger().error(f"Invalid path {in_string}.")
        return None

    if s3 is None:  # File is local
        return path
    else:
        local_path = os.path.join(os.getcwd(), os.path.basename(path))

        if os.path.isfile(local_path) and not overwrite_local:
            get_logger().warning(f"Local file {local_path} already exists for input path {in_string}.")
            return local_path

        get_logger().info(f"Attempting to download {in_string} to {local_path}")
        s3.download_file(bucket_name, path, local_path)
        s3.close()

        if os.path.isfile(local_path):
            return local_path
        else:
            get_logger().error(f"Download failed for {in_string}.")
            return None
