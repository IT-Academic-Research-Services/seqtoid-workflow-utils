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

