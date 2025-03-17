"""
Logging functions.
"""

import logging
import os


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

def setup_logger(name, log_file, level=logging.DEBUG):
    """
    Set up a logger with file and console handlers.


    :param  name: Name of the logger
    :param  log_file: Path to the log file
    :param level: Logging level

    :return  logging.Logger: Configured logger instance
    """

    # Ensure logs directory exists
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid duplicate handlers if logger is reused
    if logger.handlers:
        logger.handlers.clear()

    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(level)

    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    # Define log format
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

# Example usage within this module (optional, for testing)
if __name__ == "__main__":
    logger = setup_logger("test", "logs/test.log")
    logger.debug("This is a debug message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
