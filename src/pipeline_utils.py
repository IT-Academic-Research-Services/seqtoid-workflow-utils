"""
Common functions for running all pipelines.
"""

import sys
import subprocess
import argparse

import json
from pathlib import Path
from typing import Sequence

from snakemake.api import SnakemakeApi, ResourceSettings, OutputSettings, ConfigSettings, ExecutionSettings
from snakemake.settings.types import DefaultResources

from src.defs import *
from src.logging_utils import get_logger


# -------------------------
# Definitions
# -------------------------

STANDARD_CONFIG_FILE = "config/config.yaml"

# -------------------------
# Setup
# -------------------------

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

def run_pipeline(config, config_path=None, pipeline_name=None, dry_run=False, extra_args=None, **kwargs):
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

    try:
        project_root = config[PROJECT_ROOT_TAG]
    except KeyError:
        get_logger().critical("Error: project_root not found in config.")
        sys.exit(1)

    project_root = Path(project_root)

    snakefile = project_root / "Snakefile" if pipeline_name is None else project_root/ "workflows" / pipeline_name / "Snakefile"
    if not snakefile.exists():
        print(f"Error: Snakefile {snakefile} not found.")
        sys.exit(1)

    # Get necessary run parameters from the config file
    execution = config.get("execution", {})
    mode = execution.get("mode", "local").lower()
    cores = execution.get("cores", 1)  # NB: this determines number of simultaneous processes if local
    jobs = execution.get("jobs", 1)
    latency_wait = execution.get("latency_wait", 30)
    dry_run = execution.get("dry_run", False)


    config_files = []
    if config_path:  # e.g cluster_submit.yaml, local.yaml
        config_files.append(str(config_path))
    if pipeline_name:  # Add pipeline-specific config
        pipeline_config = project_root / "config" / f"{pipeline_name}.yaml"
        if pipeline_config.exists():
            config_files.append(str(pipeline_config))
        else:
            get_logger().warning(f"Warning: Pipeline config {pipeline_config} not found, using only {config_path}")

    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(mode="w", suffix=".json", delete=False) as temp_config:
        json.dump(config, temp_config)
        temp_config_path = temp_config.name
        config_files.append(temp_config_path)
    config_path_sequence: Sequence[Path] = [Path(s) for s in config_files]
    config_settings = ConfigSettings(configfiles=config_path_sequence, replace_workflow_config=True)


    default_resources = DefaultResources()
    default_resources.cores = cores
    default_resources.memory = execution.get("memory", 0)  # Default memory in MB
    default_resources.cpu = execution.get("cpu", 1)
    default_resources.threads = execution.get("threads", 1)
    resource_settings = ResourceSettings(
        nodes=jobs if mode == "slurm" else 1,  # Max nodes (jobs) for cluster
        default_resources=default_resources,
    )

    execution_settings = ExecutionSettings(latency_wait=latency_wait)

    executor = "local" if mode == "local" else "cluster-generic" if mode == "slurm" else None
    if not executor:
        print(f"Error: Unknown execution mode '{mode}'.")
        sys.exit(1)


    try:
        with SnakemakeApi(
                OutputSettings(
                    verbose=False,
                    show_failed_logs=True,
                    dryrun=dry_run,
                ),


        ) as snakemake_api:

            workflow_api = snakemake_api.workflow(
                snakefile=snakefile,
                workdir=project_root,
                resource_settings=resource_settings,
                config_settings=config_settings

            )

    except Exception as e:
        print(f"Error during execution: {e} type: {type(e)}")
        sys.exit(1)
    else:
        print(workflow_api)

        try:
            dag = workflow_api.dag()
        except Exception as e:
            print(f"Error during execution: {e} type: {type(e)}")
            sys.exit(1)

        #
        # success = dag.execute_workflow(
        #     workflow=workflow_api,
        #     cores=cores,
        #     nodes=jobs if executor == "cluster-generic" else None,
        #     executor=executor,
        #     dryrun=dry_run,
        #     configfiles=config_files,
        #     **kwargs
        # )
        # if not success:
        #     print(f"Failed to run {pipeline_name or 'main workflow'}.")
        #     sys.exit(1)
        # print(f"Finished running {pipeline_name or 'main workflow'}.")

