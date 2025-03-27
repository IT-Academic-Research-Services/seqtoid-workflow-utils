"""
Functions for working with FASTA and FASTQ files.
"""

import os

from src.defs import FASTQ_EXT_SET, FASTA_EXT_SET, R1_TAG_SET, FASTA_TAG, FASTQ_TAG, R1_TAG, R2_TAG, R2_TAG_SET
from src.logging_utils import get_logger


# -------------------------
# Definitions
# -------------------------


# -------------------------
# Setup
# -------------------------



# -------------------------
# Functions
# ------------------------


def fast_a_q_basename(file_name, delims={'_', '-'}):
    """
    Find the base name of a file and determine if it is a FASTA or FASTQ file.
    Will also find the location of the R1 tag in the file name and return it.
    The file name is assumed to be the base of the file itself, stripped of all enclosing directories,
    and the file extension. If an R1 tag is found, the base file name will be truncated
    to just before the R1 tag. Searches through a set of delimiters to find the R1 tag.

    :param file_name: Name of file, which does not need to be absolute.
    :param delims: Set of possible delimiters for the base anme.
    :return: File base name, R1 location or None, file type or None.
    """
    file_basename = os.path.normpath(os.path.basename(file_name)).strip()

    file_type =  None
    file_basename_no_ext = None

    r1_loc = None

    for fq_ext in FASTQ_EXT_SET:
        if file_basename.endswith(fq_ext):
            file_type = FASTQ_TAG
            file_basename_no_ext = file_basename.replace(fq_ext, '')
            break

    if file_type ==  None:
        for fa_ext in FASTA_EXT_SET:
            if file_basename.endswith(fa_ext):
                file_type = FASTA_TAG
                file_basename_no_ext = file_basename.replace(fa_ext, '')
                break

    if file_basename_no_ext is None:
        get_logger().info(f"File type is not FASTA or FASTQ for {file_name}.")
        return None, None, None

    file_delim = None
    for delim in delims:  # First, find an R1 tag in a file starting with file_base if possible
        file_base_cols = file_basename_no_ext.split(delim)

        for i in range(len(file_base_cols)):
            col = file_base_cols[i]

            if col in R1_TAG_SET:
                r1_loc = i
                file_delim = delim
                break

    if file_delim is None:
        get_logger().error("File type not recognized.")
        return None, None, None

    if r1_loc is None:
        final_basename = file_basename_no_ext
    else:
        file_base_cols = file_basename_no_ext.split(file_delim)
        final_basename = file_delim.join(file_base_cols[:r1_loc])

    return final_basename, r1_loc, file_type


def fastq_iterate(file_handle):
    """
    Iterate over a FASTQ file.
    :param file_handle: File handle.
    :return: Generator or None.
    """

    while True:
        header = file_handle.readline().strip()
        if header == '':  # File has ended on a multiple of 4
            break

        seq = file_handle.readline().strip()
        if seq == '':
            get_logger().error(f"header is {header} but no sequence follows")
            return None

        skip = file_handle.readline().strip()
        if skip != '+':
            get_logger().error(f"Third line should be a '+', got: {skip}")
            return None

        qual = file_handle.readline().strip()
        if seq == '':
            get_logger().error(f"header is {header} but no quality line follows")
            return None

        yield header, seq, qual

def fasta_iterate(file_handle):
    """
    Iterate over a FASTA file.
    :param file_handle: File handle.
    :return: Generator or None.
    """
    while True:
        header = file_handle.readline().strip()
        if header == '':
            break
        if not header.startswith('>'):
            get_logger().error(f"header is {header} but does not start with '>'")
            return None

        seq = file_handle.readline().strip()
        if seq == '':
            get_logger().error(f"header is {header} but no sequence follows")
            return None

        yield header, seq

def fast_a_q_iterate(file_handle, fastq=True):
    """
    Iterate over a FASTA or FASTQ file.
    :param file_handle: File handle.
    :param fastq: If True, iterate over FASTQ file. If False, iterate over FASTA file.
    :return: Generator or None.
    """
    if fastq:
        return fastq_iterate(file_handle)
    else:
        return fasta_iterate(file_handle)

def acquire_fast_a_q_files(working_dir, file1=None, file2=None, fastq=True, delims={'_', '-'}):
    """
    Search a given working directory for either FASTQ or FASTA files.
    The file_base is assumed to be the start of the file name itself, stripped
    of all enclosing directories.
    Will locate the R1 file first, and if it does, will try to find a matching R2 file.
    If it finds any file starting with file_base but cannot find an R1 tag, it will
    take the first such example it finds and return it as R1.

    :param working_dir: Directory in which to search.
    :param file_base:
    :param fastq: If True, search for FASTQ files, else search for FASTA files.
    :param delims: Sets of possible file name delimiters
    :return: Dict {file_base: [r1_full_path, r2_full_path], ...}, where {file_base: [r1_full_path, r2_full_path]} means paired-ended.
    """

    working_abspath = os.path.abspath(working_dir)

    if file1 is not None:
        if not os.path.isfile(file1):
            get_logger().error(f"file1 is {file1} but is not a regular file")
            return {}
        file_base, r1_loc, file_type = fast_a_q_basename(file1, delims)
        if file2 is not None:
            if not os.path.isfile(file1):
                get_logger().error(f"file2 is {file2} but is not a regular file")
                return {}

    if fastq:
        ext_set = FASTQ_EXT_SET
    else:
        ext_set = FASTA_EXT_SET

    possible_files = [f for f in os.listdir(working_dir) if f.endswith(tuple(ext_set))]

    if not possible_files:
        get_logger().error(f"Cannot fund any input files in dir {working_abspath}")
        return {}

    output_dict = {}
    added_files = set()
    for file in possible_files:
        if file in added_files:  # Skip, especially for R2 files added below
            continue
        file_basename_no_ext = None
        r1_file = None
        r1_loc = None
        file_delim = None
        ext_used = None
        file_basename = os.path.normpath(os.path.basename(file)).strip()
        for ext in ext_set:
            if file_basename.endswith(ext):
                ext_used = ext
                file_basename_no_ext = file_basename.replace(ext, '')
                break

        if file_basename_no_ext is None:
            continue

        is_r2 = False
        for delim in delims:  # First, find an R1 tag in a file starting with file_base if possible
            file_base_cols = file_basename_no_ext.split(delim)

            for i in range(len(file_base_cols)):
                col = file_base_cols[i]

                if col in R1_TAG_SET:
                    r1_loc = i
                    file_delim = delim
                    break
                elif col in R2_TAG_SET:  # Solo R2 files should be mated with r1 files or skipped
                    is_r2 = True
                    break
            if r1_loc is not None:
                break
            if is_r2:
                break

        if is_r2:
            continue

        added_files.add(file)


        if r1_loc is None:
            output_dict[file_basename_no_ext] = [os.path.join(working_abspath, file), None]

        else:  # search for r2
            r1_file = file
            sample_base_cols = file_basename_no_ext.split(file_delim)[:r1_loc]
            sample_base = file_delim.join(sample_base_cols)

            file_base_cols = file_basename_no_ext.split(file_delim)
            file_base_cols[r1_loc] = R2_TAG
            r2_file_base = file_delim.join(file_base_cols)
            r2_file_name = r2_file_base + ext_used
            if os.path.isfile(os.path.join(working_abspath, r2_file_name)):
                output_dict[sample_base] = [os.path.join(working_abspath, file), os.path.join(working_abspath, r2_file_name)]
                added_files.add(r2_file_name)
            else:
                output_dict[sample_base] = [os.path.join(working_abspath, file), None]

    return output_dict


