"""
Logging functions.
"""

import logging
import os


# -------------------------
# Setup
# -------------------------

_logger = None


# -------------------------
# Setup
# -------------------------

LOGGER_NAME = 'cypherid'
LOGGER_DEFAULT_FILE = 'logs/cypherid.log'


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

def get_logger(log_file=LOGGER_DEFAULT_FILE, level=logging.DEBUG):
    """
    Set up a logger with file and console handlers.


    :param  name: Name of the logger
    :param  log_file: Path to the log file
    :param level: Logging level

    :return  logging.Logger: Configured logger instance
    """

    global _logger

    if _logger is None:

        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        _logger = logging.getLogger(LOGGER_NAME)
        _logger.setLevel(level)

        if _logger.handlers:  # Avoid duplicate handlers if logger is reused
            _logger.handlers.clear()


        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)

        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        _logger.addHandler(file_handler)
        _logger.addHandler(console_handler)

    return _logger

def set_log_file(log_file, level=logging.DEBUG):
    get_logger(log_file, level)
