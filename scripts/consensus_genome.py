"""
Consensus genome run script.
"""

import os
import logging
import json
import argparse
from pathlib import Path
from datetime import datetime

from src.defs import *
from src.config_utils import setup_config
from src.logging_utils import get_logger, set_log_file
from src.pipeline_utils import run_pipeline, common_parser
from src.fastaq import acquire_fast_a_q_files


# -------------------------
# Definitions
# -------------------------

PIPELINE_NAME = "consensus_genome"
INPUT_DIR = "data/input"
TECHNOLOGY_TAG = "technology"
ILLUMINA_TAG = "illumina"
ONT_TAG = "ont"


# -------------------------
# Setup
# -------------------------


# -------------------------
# Functions
# -------------------------

def parse_arguments():
    """
    PArse args specific to this pipeline.
    :return: parser.parse_args()
    """

    parser = argparse.ArgumentParser(parents=[common_parser()])
    parser.add_argument('-i', '-in1', '--input_fastq1', type=str, required=False, help='Path to input FASTQ file 1')
    parser.add_argument('-I', '-in2', '--input_fastq2', type=str, required=False, help='Path to input FASTQ file 2')
    parser.add_argument('--technology', type=str, default=ILLUMINA_TAG, help='Sequencing tech. Options: illumina, ont')
    return parser.parse_args()


# -------------------------
# Pipeline
# -------------------------

args = parse_arguments()

input_dict = acquire_fast_a_q_files(INPUT_DIR, args.input_fastq1, args.input_fastq2, fastq=True)
if input_dict is None:
    get_logger().critical("Error: No input files found.")
    exit(1)

config, config_path = setup_config(args, PIPELINE_NAME)

config[TECHNOLOGY_TAG] = args.technology

input_json = json.dumps(input_dict)
config[INPUT_DICT_TAG] = input_json

get_logger().info("Starting consensus genome pipeline")

project_root = Path(config["project_root"])

snakefile = project_root/ "workflows" / PIPELINE_NAME / "Snakefile"
if not snakefile.exists():
    print(f"Error: Snakefile {snakefile} not found.")
    exit(1)

run_pipeline(config=config, config_path=config_path, pipeline_name=PIPELINE_NAME)

get_logger().info("Finished consensus genome pipeline")

exit(0)
