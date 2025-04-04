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
FASTA_TAG = 'FASTA'
FASTQ_TAG = 'FASTQ'
FAST_A_Q_DELIMS = {'_', '-', '.', '@', '#', ':'}
R1_TAG = 'R1'
R2_TAG = 'R2'

# Sets for detecting file types
FASTA_EXT_SET = {'.fa', '.fasta', '.fa.gz', '.fasta.gz'}
FASTQ_EXT_SET = {'.fq', '.fastq', '.fastq.gz', '.fq.gz'}
R1_TAG_SET = {'R1', 'r1', '1'}
R2_TAG_SET = {'R2', 'r2', '2'}
R1_R2_TAGS = {'R1': 'R2', 'r1': 'r2', '1': '2', 'F': 'R','f': 'r', 'FWD':'REV', 'fwd': 'rev', 'PE1': 'PE2', 'pe1': 'pe2', 'READ1': 'READ2', 'read1': 'read2'}
R1_TAG_PRIORITY = {
    'R1': 1, 'r1': 2,  # R1/r1 highest
    'PE1': 3, 'pe1': 4,  # PE1/pe1 next
    'READ1': 5, 'read1': 6,  # READ1/read1
    'FWD': 7, 'fwd': 8,  # FWD/fwd
    '1': 9,  # Numeric '1'
    'F': 10, 'f': 11  # Least specific
}

# ---------------
# Bio definitions
# ---------------


# ---------------
# Pipeline definitions
# ---------------

PROJECT_ROOT_TAG = 'project_root'
LOG_PATH_TAG = 'log_path'
INPUT_DICT_TAG = 'input_dict'

