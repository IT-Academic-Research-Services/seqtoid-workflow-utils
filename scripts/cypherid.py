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

parser = argparse.ArgumentParser()
parser.add_argument("--log-level", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
parser.add_argument("-p", "--pipeline", help="Workflow pipeline file to run (e.g., 'consensus-genome') Defaults to the main pipeline.", default=None)
parser.add_argument("-c", "--config_file", default=None, choices=["local", "cluster", "cluster_submit"])
parser.add_argument("--dry-run", action="store_true", help="Perform a dry run")

args = parser.parse_args()
level = getattr(logging, args.log_level.upper())

config = setup_config(PROJECT_ROOT, args.config_file)
logger = setup_logger("cypherid", LOG_FILENAME, level=config.get("logging", {}).get("level", "INFO"))

# -------------------------
# Functions
# -------------------------

def run_pipeline(pipeline_name=None, dry_run=False, **kwargs):
    """
        Run a CypherID workflow.
        :param pipeline_name: Name of the workflow file (e.g., 'workflow1.smk') or None for main Snakefile
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
