"""
Configuration functions shared across pipelines.
"""

import os
import yaml

# -------------------------
# Definitions
# -------------------------

BASE_CONFIG_FILE = 'config.yaml'

# -------------------------
# Functions
# -------------------------

def setup_config(project_root, config_file):
    """
    Sets up the config file by reading from a YAML in the config dir.
    :param project_root: Absolute path to the project root folder.
    :param config_file: Name without the extension of YAML config file
    :return: Config dict.
    """

    if config_file is None:
        config_name = BASE_CONFIG_FILE
    else:
        config_name = config_file + '.yaml'

    config_path = project_root / "config" / config_name

    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    else:
        print(f"Config file {config_file} not found")
        config = {}

    return config, config_path
