"""
Consensus genome run script.
"""

import sys
import logging
from pathlib import Path
from src.logging_utils import setup_logger
from src.config_utils import setup_config
from src.pipeline_utils import run_pipeline, parse_arguments


# -------------------------
# Definitions
# -------------------------

LOG_FILENAME = "logs/consensus-genome.log"
PIPELINE_NAME = "consensus-genome"


# -------------------------
# Setup
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from scripts/

args = parse_arguments()

config = setup_config(PROJECT_ROOT, args.config_file)
if args.log_level is None:
    log_level = getattr(logging, config.get("logging", {}).get("level", "INFO").upper())
else:
    log_level = getattr(logging, args.log_level.upper())
logger = setup_logger("cypherid", LOG_FILENAME, level=log_level)

logger.info(f"Starting consensus-genome pipeline.")


# -------------------------
# Pipeline
# -------------------------

snakefile =  PROJECT_ROOT / "workflows" / PIPELINE_NAME / "Snakefile"
if not snakefile.exists():
    print(f"Error: Snakefile {snakefile} not found.")
    sys.exit(1)

run_pipeline(logger=logger, project_root=PROJECT_ROOT, pipeline_name=PIPELINE_NAME, dry_run=args.dry_run)

# with open(snakemake.output[0], "w") as f:
#     f.write("Done")

logger.info(f"Consensus-genome run completed.")


# -------------------------
# Functions
# -------------------------

