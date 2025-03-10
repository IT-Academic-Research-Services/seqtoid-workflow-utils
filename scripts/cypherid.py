"""
CypherID main run script.
"""

import os
import sys
import logging
import argparse
import subprocess
import yaml
from pathlib import Path
from src.logging_utils import setup_logger


# -------------------------
# Definitions
# -------------------------

AVAILABLE_PIPELINES = ["consensus-genome"]


# -------------------------
# Setup
# -------------------------

# Get the project root directory (where run_workflows.py resides)
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/

parser = argparse.ArgumentParser()
parser.add_argument("--log-level", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
parser.add_argument("-p", "--pipeline", help="Workflow pipeline file to run (e.g., 'consensus-genome') Defaults to the main pipeline.", default=None)
parser.add_argument("--dry-run", action="store_true", help="Perform a dry run")
# parser.add_argument("--cores", type=int, default=1, help="Number of cores to use")

args = parser.parse_args()
level = getattr(logging, args.log_level.upper())

# Use a relative log file path; setup_logger will place it in cwd/logs/
log_file = "logs/cypherid.log"

# Optionally load config for level, but log dir is handled by setup_logger
config_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config", "config.yaml")
if os.path.exists(config_file):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
else:
    config = {}

logger = setup_logger("cypherid", log_file, level=config.get("logging", {}).get("level", "INFO"))


# -------------------------
# Functions
# -------------------------

def run_pipeline(pipeline_name=None, dry_run=False, **kwargs):
    """
        Run a CypherID workflow.
        :param wpipeline_name: Name of the workflow file (e.g., 'workflow1.smk') or None for main Snakefile
        :param dry_run: If True, perform a dry run (-n flag)
        :param kwargs: Additional Snakemake CLI arguments
        """

    logger.info(f"Attempting to run pipeline: {pipeline_name}")

    snakefile = PROJECT_ROOT / "Snakefile" if pipeline_name is None else PROJECT_ROOT / "workflows" / pipeline_name / "Snakefile"
    if not snakefile.exists():
        print(f"Error: Snakefile {snakefile} not found.")
        sys.exit(1)

    cmd = ["snakemake", "--snakefile", str(snakefile)]
    if dry_run:
        cmd.append("-n")  # Dry run
    for key, value in kwargs.items():
        cmd.extend([f"--{key}", str(value)])

    try:
        subprocess.run(cmd, shell=False, check=True)
        logger.info(f"Successfully ran {pipeline_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run {pipeline_name}: {str(e)}")

if __name__ == "__main__":

    print("\n------------\n  CypherID\n------------\n")
    if args.pipeline is None:
        logger.debug(f"Available pipelines: {AVAILABLE_PIPELINES}")

        print("Available pipelines:", ", ".join(AVAILABLE_PIPELINES))
        choice = input("Enter pipeline name to run or Enter for main pipeline: ").strip()
    else:
        choice = args.pipeline

    if len(choice) < 1:
        run_pipeline(pipeline_name=None, dry_run=args.dry_run)
    elif choice in AVAILABLE_PIPELINES:
        run_pipeline(pipeline_name=choice, dry_run=args.dry_run)
    else:
        logger.warning(f"Invalid pipeline name entered: {choice}")
        print("Invalid pipeline name!")

    exit(0)
