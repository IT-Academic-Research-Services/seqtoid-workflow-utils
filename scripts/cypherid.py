"""
CypherID main run script.
"""

import sys
import logging
import argparse
import subprocess
from pathlib import Path
from src.logging_utils import setup_logger
from src.config_utils import setup_config
from src.pipeline_utils import run_pipeline, parse_arguments


# -------------------------
# Definitions
# -------------------------

AVAILABLE_PIPELINES = ["consensus-genome"]
LOG_FILENAME = "logs/cypherid.log"

# -------------------------
# Setup
# -------------------------

# Get the project root directory (where run_workflows.py resides)
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/

args = parse_arguments()

config = setup_config(PROJECT_ROOT, args.config_file)
if args.log_level is None:
    log_level = getattr(logging, config.get("logging", {}).get("level", "INFO").upper())
else:
    log_level = getattr(logging, args.log_level.upper())
logger = setup_logger("cypherid", LOG_FILENAME, level=log_level)

# -------------------------
# Functions
# -------------------------

if __name__ == "__main__":

    print("\n------------\n  CypherID\n------------\n")
    if args.pipeline is None:
        logger.debug(f"Available pipelines: {AVAILABLE_PIPELINES}")

        print("Available pipelines:", ", ".join(AVAILABLE_PIPELINES))
        choice = input("Enter pipeline name to run or Enter for main pipeline: ").strip()
    else:
        choice = args.pipeline

    if len(choice) < 1:
        run_pipeline(logger=logger, project_root=PROJECT_ROOT, pipeline_name=None, dry_run=args.dry_run)
    elif choice in AVAILABLE_PIPELINES:
        run_pipeline(logger=logger, project_root=PROJECT_ROOT, pipeline_name=choice, dry_run=args.dry_run)
    else:
        logger.warning(f"Invalid pipeline name entered: {choice}")
        print("Invalid pipeline name!")

    exit(0)
