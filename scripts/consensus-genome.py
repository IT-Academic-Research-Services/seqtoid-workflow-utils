"""
Consensus genome run script.
"""

import sys
import logging
import argparse
from pathlib import Path

from src.config_utils import setup_config
from src.pipeline_utils import run_pipeline, common_parser


# -------------------------
# Definitions
# -------------------------

PIPELINE_NAME = "consensus-genome"


# -------------------------
# Setup
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/


# -------------------------
# Functions
# -------------------------

def parse_arguments():
    """
    PArse args specific to this pipeline.
    :return: parser.parse_args()
    """

    parser = argparse.ArgumentParser(parents=[common_parser()])
    parser.add_argument('-i', '-in1', '--input_fastq1', type=str, required=True, help='Path to input FASTQ file 1')
    return parser.parse_args()


# -------------------------
# Pipeline
# -------------------------

args = parse_arguments()
config = setup_config(PROJECT_ROOT, args.config_file)
if args.log_level is None:
    log_level = getattr(logging, config.get("logging", {}).get("level", "INFO").upper())
else:
    log_level = getattr(logging, args.log_level.upper())


snakefile = PROJECT_ROOT / "workflows" / PIPELINE_NAME / "Snakefile"
if not snakefile.exists():
    print(f"Error: Snakefile {snakefile} not found.")
    sys.exit(1)

run_pipeline(project_root=PROJECT_ROOT, log_level=log_level, pipeline_name=PIPELINE_NAME, dry_run=args.dry_run, extra_args=args.extra_args)

exit(0)
