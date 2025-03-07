"""
CypherID main run script.
"""

import os
import logging
import argparse
import subprocess
from src.logging_utils import setup_logger

parser = argparse.ArgumentParser()
parser.add_argument("--log-level", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
args = parser.parse_args()
level = getattr(logging, args.log_level.upper())
logger = setup_logger("cypherid", "logs/cypherid.log")

def run_pipeline(pipeline_name):
    logger.info(f"Attempting to run pipeline: {pipeline_name}")
    cmd = f"snakemake -s workflows/{pipeline_name}/Snakefile --cores 8"
    try:
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f"Successfully ran {pipeline_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run {pipeline_name}: {str(e)}")

if __name__ == "__main__":
    pipelines = ["consensus-genome"]
    logger.debug(f"Available pipelines: {pipelines}")
    print("\n-----------\n  CypherID\n  UCSF\n-----------\n")
    print("Available pipelines:", ", ".join(pipelines))
    choice = input("Enter pipeline name to run: ").strip()
    if choice in pipelines:
        run_pipeline(choice)
    else:
        logger.warning(f"Invalid pipeline name entered: {choice}")
        print("Invalid pipeline name!")
