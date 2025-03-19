"""
Common definitions.
"""


# ---------------
# Character definitions
# ---------------

GZIP_MAGIC_NUMBER = b'\x1f\x8b'


# ---------------
# General File definitions
# ---------------

GZIP_EXT_SET = {'.gz', '.gzip', '.tgz'}


# ---------------
# Bio File definitions
# ---------------

# Canonical for creating files
FASTA_EXT = '.fasta'
FASTQ_EXT = '.fastq'
R1_TAG = 'R1'
R2_TAG = 'R2'

# Sets for detecting file types
FASTA_EXT_SET = {'.fa', '.fasta', '.fa.gz', '.fasta.gz'}
FASTQ_EXT_SET = {'.fq', '.fastq', '.fastq.gz', '.fq.gz'}
R1_TAG_SET = {'R1', 'r1'}
R2_TAG_SET = {'R2', 'r2'}


# ---------------
# Bio definitions
# ---------------


# ---------------
# Pipeline definitions
# ---------------

PROJECT_ROOT_NAME = "project_root"
