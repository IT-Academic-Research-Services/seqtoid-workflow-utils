"""
Logging functions.
"""

import logging
import os
from datetime import datetime


# -------------------------
# Functions
# -------------------------

def snakemake_log_level(log_level):
    """
    Determine the appropriate Snakemake log level flag based on the logging level.
    :param log_level:
    :return: log flag suitable for a snakemake cli invocation
    """

    log_flag = ""
    if log_level <= logging.DEBUG:  # 10
        log_flag = "--debug"  # Very verbose, includes debug info
    elif log_level <= logging.INFO:  # 20
        log_flag = "--verbose"  # Show detailed execution info
    elif log_level <= logging.WARNING:  # 30
        log_flag = ""  # Default behavior (moderate output)
    else:  # ERROR (40) or higher
        log_flag = "--quiet"  # Minimal output
    return log_flag
