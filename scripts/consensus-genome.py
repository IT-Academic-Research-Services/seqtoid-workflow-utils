"""
Consensus genome run script.
"""

from src.logging_utils import setup_logger

# Use relative log file path; setup_logger will place it in cwd/logs/
log_file = "logs/consensus-genome.log"

logger = setup_logger("consensus-genome", log_file, level=snakemake.config.get("logging", {}).get("level", "INFO"))
logger.info(f"Starting consensus genomes with input: {snakemake.input[0]}")

logger.info(f"Consensus genome run completed.")

with open(snakemake.output[0], "w") as f:
    f.write("Done")