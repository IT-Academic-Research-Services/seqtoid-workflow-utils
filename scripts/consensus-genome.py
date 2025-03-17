"""
Consensus genome run script.
"""

import os
import logging
import argparse
from pathlib import Path
from datetime import datetime

from src.config_utils import setup_config
from src.logging_utils import setup_logger
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

log_dir = os.path.join(os.getcwd(), 'logs')
log_date = datetime.now().strftime('%Y%m%d_%H%M%S')
log_file = "consensus_genome_" + log_date +".log"
log_path = os.path.join(log_dir, log_file)
logger = setup_logger("consensus_genome", log_path, level=log_level)

logger.info("Starting consensus genome pipeline")

snakefile = PROJECT_ROOT / "workflows" / PIPELINE_NAME / "Snakefile"
if not snakefile.exists():
    print(f"Error: Snakefile {snakefile} not found.")
    exit(1)

run_pipeline(project_root=PROJECT_ROOT, logger=logger, log_path=log_path, pipeline_name=PIPELINE_NAME, dry_run=args.dry_run, extra_args=args.extra_args)

logger.info("Finished consensus genome pipeline")

exit(0)
