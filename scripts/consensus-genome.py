"""
Consensus genome run script.
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

LOG_FILENAME = "logs/consensus-genome.log"
PIPELINE_NAME = "consensus-genome"

# -------------------------
# Setup
# -------------------------


PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/

parser = argparse.ArgumentParser()
parser.add_argument("--log-level", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
parser.add_argument("-p", "--pipeline", help="Workflow pipeline file to run (e.g., 'consensus-genome') Defaults to the main pipeline.", default=None)
parser.add_argument("-c", "--config_file", default=None, choices=["local", "cluster", "cluster_submit"])
parser.add_argument("--dry-run", action="store_true", help="Perform a dry run")

args = parser.parse_args()
level = getattr(logging, args.log_level.upper())

config = setup_config(PROJECT_ROOT, args.config_file)
logger = setup_logger("consensus-genome", LOG_FILENAME, level=snakemake.config.get("logging", {}).get("level", "INFO"))
logger.info(f"Starting consensus genomes with input: {snakemake.input[0]}")

snakefile =  PROJECT_ROOT / "workflows" / PIPELINE_NAME / "Snakefile"
if not snakefile.exists():
    print(f"Error: Snakefile {snakefile} not found.")
    sys.exit(1)



with open(snakemake.output[0], "w") as f:
    f.write("Done")

logger.info(f"Consensus genome run completed.")
