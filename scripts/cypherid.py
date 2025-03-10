"""
CypherID main run script.
"""

import os
import logging
import argparse
import subprocess
import yaml
from src.logging_utils import setup_logger

parser = argparse.ArgumentParser()
parser.add_argument("--log-level", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
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
