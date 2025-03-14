"""
CypherID main run script.
"""

import sys
import logging
import argparse
import subprocess
from pathlib import Path
from src.config_utils import setup_config
from src.pipeline_utils import run_pipeline, common_parser


# -------------------------
# Definitions
# -------------------------

AVAILABLE_PIPELINES = ["consensus-genome"]

# -------------------------
# Setup
# -------------------------

# Get the project root directory (where run_workflows.py resides)
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/

parser = argparse.ArgumentParser(parents=[common_parser()])  # Use common_parser() from pipeline_utils.py only
args = parser.parse_args()

config = setup_config(PROJECT_ROOT, args.config_file)
if args.log_level is None:
    log_level = getattr(logging, config.get("logging", {}).get("level", "INFO").upper())
else:
    log_level = getattr(logging, args.log_level.upper())

# -------------------------
# Functions
# -------------------------

def mock_test(test_num):
    """
    Simple function here as a stub to set up the pytest framework.
    :param test_num: Any integer
    :return: test_num * 2 if an integer is passed.
    """

    return test_num * 2

def main():

    print("\n------------\n  CypherID\n------------\n")
    if args.pipeline is None:

        print("Available pipelines:", ", ".join(AVAILABLE_PIPELINES))
        choice = input("Enter pipeline name to run or Enter for main pipeline: ").strip()
    else:
        choice = args.pipeline

    if len(choice) < 1:
        run_pipeline(project_root=PROJECT_ROOT, log_level=log_level, pipeline_name=None, dry_run=args.dry_run,
                     extra_args=args.extra_args)
    elif choice in AVAILABLE_PIPELINES:
        run_pipeline(project_root=PROJECT_ROOT, log_level=log_level, pipeline_name=choice, dry_run=args.dry_run, extra_args=args.extra_args)
    else:
        print("Invalid pipeline name!")
        exit(1)


if __name__ == "__main__":

    main()

    exit(0)
