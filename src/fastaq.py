"""
Functions for working with FASTA and FASTQ files.
"""


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
