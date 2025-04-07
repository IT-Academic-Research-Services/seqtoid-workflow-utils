"""
Configuration functions shared across pipelines.
"""

import os
import logging
import yaml
from pathlib import Path
from datetime import datetime

from src.defs import *
from src.logging_utils import set_log_file

# -------------------------
# Definitions
# -------------------------

BASE_CONFIG_FILE = 'config.yaml'

# -------------------------
# Setup
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # Two levels up from src/

# -------------------------
# Functions
# -------------------------

def setup_config(args, pipeline_name):
    """
    Sets up the config file by reading from a YAML in the config dir.
    :param project_root: Absolute path to the project root folder.
    :param config_file: Name without the extension of YAML config file
    :return: Config dict.
    """

    if args.config_file is None:
        config_name = BASE_CONFIG_FILE
    else:
        config_name = args.config_file + '.yaml'

    config_path = PROJECT_ROOT / "config" / config_name

    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    else:
        print(f"Config file {args.config_file} not found")
        config = {}

    # Set up other params
    if args.log_level is None:
        log_level = getattr(logging, config.get("logging", {}).get("level", "INFO").upper())
    else:
        log_level = getattr(logging, args.log_level.upper())
    log_dir = os.path.join(os.getcwd(), 'logs')
    log_date = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = pipeline_name + "_" + log_date + ".log"
    log_path = os.path.join(log_dir, log_file)
    set_log_file(log_path)
    config[LOG_PATH_TAG] = str(log_path)

    config[PROJECT_ROOT_TAG] = str(PROJECT_ROOT)

    return config, config_path
