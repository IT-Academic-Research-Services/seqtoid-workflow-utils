"""
Log handler for snakemake logs.
"""
import sys
import json
import os
import logging

# Get the log file path from the environment variable
log_file = os.environ.get("SNAKEMAKE_LOG_FILE")
if not log_file:
    raise ValueError("SNAKEMAKE_LOG_FILE environment variable not set")

# Configure the logger
logger = logging.getLogger("snakemake")
handler = logging.FileHandler(log_file)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.handlers = []
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# Process log messages from Snakemake (sent as JSON lines via stdin)
for line in sys.stdin:
    msg = json.loads(line)
    level = msg["level"].lower()
    if level == "debug":
        logger.debug(msg["msg"])
    elif level == "info":
        logger.info(msg["msg"])
    elif level == "warning":
        logger.warning(msg["msg"])
    elif level == "error":
        logger.error(msg["msg"])
    else:
        logger.info(f"Unknown log level '{level}': {msg['msg']}")