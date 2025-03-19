"""
Functions for manipulating files and file names.
"""

import os

from src.defs import GZIP_EXT_SET

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
