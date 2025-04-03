"""
Common functions for running all pipelines.
"""

import sys
import subprocess
import argparse
import json

from src.logging_utils import get_logger


# -------------------------
# Definitions
# -------------------------

STANDARD_CONFIG_FILE = "config/config.yaml"


# -------------------------
# Functions
# -------------------------

def common_parser():
    """
    Parse command line arguments for all pipelines.
    :return: parser.parse_args() output
    """

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--log-level", default=None, choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("-p", "--pipeline",
                        help="Workflow pipeline file to run (e.g., 'consensus_genome') Defaults to the main pipeline.",
                        default=None)
    parser.add_argument("-c", "--config_file", default=None, choices=["local", "cluster", "cluster_submit"])
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run")
    parser.add_argument("extra_args", nargs=argparse.REMAINDER, help="Additional arguments to pass to Snakemake")

    return parser

def run_pipeline(project_root, log_path, config_dict, input_dict=None, config_path=None, pipeline_name=None, dry_run=False, extra_args=None, **kwargs):
    """
        Run a CypherID workflow.
        :param project_root: project root directory.
        :param log_path: Absolute path to the log file.
        :param pipeline_name: Name of the workflow file (e.g., 'workflow1.smk') or None for main Snakefile
        :param dry_run: If True, perform a dry run (-n flag)
        :param extra_args: Additional arguments to pass to Snakemake
        :param kwargs: Additional Snakemake CLI arguments
        """

    get_logger().info("Starting run_pipeline")

    snakefile = project_root / "Snakefile" if pipeline_name is None else project_root/ "workflows" / pipeline_name / "Snakefile"
    if not snakefile.exists():
        print(f"Error: Snakefile {snakefile} not found.")
        sys.exit(1)

    # Get necessary run parameters from the config file
    execution = config_dict.get("execution", {})
    mode = execution.get("mode", "local").lower()
    cores = execution.get("cores", 1)  # NB: this determines number of simultaneous processes if local
    jobs = execution.get("jobs", 1)
    latency_wait = execution.get("latency_wait", 30)
    dry_run = execution.get("dry_run", False)

    config_args = [f"project_root={project_root}", f"log_path={log_path}"]
    if input_dict is not None:
        input_json = json.dumps(input_dict)
        config_args.append(f"input_dict={input_json}")

    cmd = [
        "snakemake",
        "--snakefile", snakefile,
        "--config", *config_args,
    ]

    config_files = []
    if config_path:  # e.g cluster_submit.yaml, local.yaml
        config_files.append(str(config_path))
    if pipeline_name:  # Add pipeline-specific config
        pipeline_config = project_root / "config" / f"{pipeline_name}.yaml"
        if pipeline_config.exists():
            config_files.append(str(pipeline_config))
        else:
            get_logger().warning(f"Warning: Pipeline config {pipeline_config} not found, using only {config_path}")
    if config_files:
        cmd.extend(["--configfile"] + config_files)

    if mode == "slurm":
        cmd.extend(["--executor", "slurm"])
        cmd.extend(["--jobs", str(jobs)])
        cmd.extend(["--cores", str(cores)])  # Cores for the main process
        cmd.extend(["--latency-wait", str(latency_wait)])
        if "cluster_config" in execution:
            cmd.extend(["--cluster-config", str(project_root / execution["cluster_config"])])

    elif mode == "local":
        cmd.extend(["--cores", str(cores)])
    else:
        print(f"Error: Unknown execution mode '{mode}' in config. Use 'local' or 'slurm'.")
        sys.exit(1)

    if dry_run:
        cmd.append("-n")  # Dry run
    for key, value in kwargs.items():
        cmd.extend([f"--{key}", str(value)])
    if extra_args is not None:
        cmd.extend(extra_args)

    try:
        subprocess.run(cmd, shell=False, check=True)
    except subprocess.CalledProcessError as e:
        get_logger().critical(f"Failed to run {pipeline_name}: {str(e)}")
        sys.exit(1)
    else:
        get_logger().info(f"Finished run_pipeline {pipeline_name}")
