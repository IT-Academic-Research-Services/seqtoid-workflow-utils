"""
Functions for working with FASTA and FASTQ files.
"""

import os

from src.defs import FASTQ_EXT_SET, FASTA_EXT_SET, R1_TAG_SET, FASTA_TAG, FASTQ_TAG, R1_TAG, R2_TAG, R2_TAG_SET, R1_R2_TAGS, FAST_A_Q_DELIMS, R1_TAG_PRIORITY
from src.logging_utils import get_logger


# -------------------------
# Definitions
# -------------------------

ONE_TAG = '1'

# -------------------------
# Setup
# -------------------------



# -------------------------
# Functions
# ------------------------


def fast_a_q_basename(file_name):
    """
    Find the base name of a file and determine if it is a FASTA or FASTQ file.
    Will also find the location of the R1 tag in the file name and return it.
    The file name is assumed to be the base of the file itself, stripped of all enclosing directories,
    and the file extension. If an R1 tag is found, the base file name will be truncated
    to just before the R1 tag. Searches through a set of delimiters to find the R1 tag.

    :param file_name: Name of file, which does not need to be absolute.
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
    for delim in FAST_A_Q_DELIMS:  # First, find an R1 tag in a file starting with file_base if possible
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

def acquire_fast_a_q_files(working_dir, file1=None, file2=None, fastq=True):
    """
    Search a given working directory for either FASTQ or FASTA files.
    Applicable to: short read, consensus genome , host genome generation, long read, and phylotree-ng pipelines.

    If given file1, verifies its presence, and then check for file2 if it is also not None. Get the base name
    of the R1 file. Return a dict of just the R1 file, and if file2 is not None, the R2 file. The key will be the
    base name.

    If file1 is None, it will search the directory with the following steps:
    1. List all files in the directory matching all extensions in the set for FASTQ or FASTQ.
    2. For each file, get a base name with extension removed.
    3. For each of the given delimiters, split the base name and look for an R1 tag. If found, mark the location
    of the R1 tag.
    4. If an R1 tag is found, look for file with a matching R2 tag by replacing the R1 tag with the R2 tag.
    5. If a matching R2 file is found, assign it to the R1 file

    Thus, for each file, there will be three outcome and three possible dictionary entry types:
    a: {file_base: [r1_full_path, None]} - R1 file only, R1 tag found, presumed single ended, but possibly an error
    b: file_base: [r1_full_path, None]} - R1 file only, R1 tag not found, presumed single ended
    c: {file_base: [r1_full_path, r2_full_path]} - R1 and R2 file, paired ended

    NB:
    Files identified as R2 files without a matching R1 will be ignored.
    Files are not expected to have more than one valid R1 or R2 tag. e.g. sample_R2_extra_R1.fastq is disallowed.
    Only exactly matching extensions allowed. i.e. sample1_R1.fq does NOT pair with sample1_R2.fq.gz.

    :param working_dir: Directory in which to search.
    :param file1: R1 or single-ended file, if only one file/pair specified.
    :param file2: R2 file if paired ended.
    :param fastq: If True, iterate over FASTQ files. Else, assume FASTA files.
    :return: Dict {file_base: [r1_full_path, r2_full_path], ...}, where {file_base: [r1_full_path, r2_full_path]} means paired-ended.
    """

    working_abspath = os.path.abspath(working_dir)

    if not os.path.isdir(working_abspath):
        get_logger().error(f"Working directory {working_abspath} is not a directory")
        return {}

    if file1 is not None:
        file1_path = os.path.join(working_abspath, file1)
        if not os.path.isfile(file1_path):
            get_logger().error(f"file1 absolute path is {file1_path} but is not a regular file")
            return {}
        file_base, r1_loc, file_type = fast_a_q_basename(file1)
        if file2 is None:
            return {file_base: [file1_path, None]}
        else:
            file2_path = os.path.join(working_abspath, file2)
            if not os.path.isfile(file2_path):
                get_logger().error(f"file2 is {file2} but is not a regular file")
                return {}
            return {file_base: [os.path.join(working_abspath, file1), os.path.join(working_abspath, file2)]}

    if fastq:
        ext_set = FASTQ_EXT_SET
    else:
        ext_set = FASTA_EXT_SET

    possible_files = [f for f in os.listdir(working_dir) if f.endswith(tuple(ext_set))]

    if not possible_files:
        get_logger().error(f"Cannot find any input files in dir {working_abspath}")
        return {}

    output_dict = {}
    added_files = set()
    for file in possible_files:
        if file in added_files:  # Skip, done for R2 files added below
            continue
        file_basename_no_ext = None
        ext_used = None
        file_basename = os.path.normpath(os.path.basename(file)).strip()
        for ext in ext_set:
            if file_basename.endswith(ext):
                ext_used = ext
                file_basename_no_ext = file_basename.replace(ext, '')
                break

        if file_basename_no_ext is None:  # A file without a recognized extension is skipped
            added_files.add(file)
            continue

        best_r1_priority = float('inf')  # Lower is better, start with infinity
        best_r1_loc = None
        best_r1_tag = None
        best_delim = None
        is_r2 = False

        # Scan all delimiters to find the highest-priority R1 tag
        for delim in FAST_A_Q_DELIMS:
            file_base_cols = file_basename_no_ext.split(delim)
            found_r2 = False

            for i in range(len(file_base_cols) - 1, -1, -1):  # Reverse search
                col = file_base_cols[i]
                if col in R1_R2_TAGS.values():  # Check for R2 first to skip early
                    found_r2 = True
                    break
                if col in R1_R2_TAGS.keys():
                    priority = R1_TAG_PRIORITY.get(col, 100)  # Default high for unlisted tags
                    if priority < best_r1_priority:  # Higher priority found
                        best_r1_priority = priority
                        best_r1_loc = i
                        best_r1_tag = col
                        best_delim = delim
                    break  # Move to next delimiter after finding an R1 tag

            if found_r2:  # If R2 tag found, mark as R2-only and stop checking this file
                is_r2 = True
                break

        if is_r2:  # Skip R2-only files
            continue

        added_files.add(file)

        if best_r1_loc is None:  # No R1 tag found, single-ended
            output_dict[file_basename_no_ext] = [os.path.join(working_abspath, file), None]
        else:  # Search for R2 using the highest-priority R1 tag
            r1_file = file
            try:
                r2_tag = R1_R2_TAGS[best_r1_tag]  # Get the mapped R2 tag
            except KeyError:
                get_logger().error(f"R1 tag {best_r1_tag} not found while searching for R2. Function should not have reached this point.")
                continue
            sample_base_cols = file_basename_no_ext.split(best_delim)[:best_r1_loc]
            sample_base = best_delim.join(sample_base_cols)

            file_base_cols = file_basename_no_ext.split(best_delim)
            file_base_cols[best_r1_loc] = r2_tag  # Use the mapped R2 tag, not R2_TAG
            r2_file_base = best_delim.join(file_base_cols)
            r2_file_name = r2_file_base + ext_used
            if os.path.isfile(os.path.join(working_abspath, r2_file_name)):
                output_dict[sample_base] = [os.path.join(working_abspath, file), os.path.join(working_abspath, r2_file_name)]
                added_files.add(r2_file_name)
            else:
                output_dict[sample_base] = [os.path.join(working_abspath, file), None]

    return output_dict
