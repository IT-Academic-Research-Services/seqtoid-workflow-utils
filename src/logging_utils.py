"""
Logging functions.
"""

import logging
import os
from datetime import datetime


def setup_logger(name, log_file, level=logging.DEBUG):
    """
    Set up a logger with file and console handlers.

    Args:
        name (str): Name of the logger (e.g., 'pipeline1', 'run_all')
        log_file (str): Path to the log file (e.g., 'logs/pipeline1.log')
        level (int): Logging level (e.g., logging.DEBUG, logging.WARNING)

    Returns:
        logging.Logger: Configured logger instance
    """

    # Use current working directory
    cwd = os.getcwd()
    # If log_file is absolute, use it as-is; otherwise, place in cwd/logs/
    if not os.path.isabs(log_file):
        log_file = os.path.join(cwd, "logs", os.path.basename(log_file))

    # Ensure logs directory exists in the current working directory
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