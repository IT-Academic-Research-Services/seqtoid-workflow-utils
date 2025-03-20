"""
Common functions for running all pipelines.
"""
import logging
import os
import sys
import subprocess
import argparse
import snakemake
from src.logging_utils import get_logger
from src.defs import FASTQ_EXT_SET, FASTA_EXT_SET, R1_TAG_SET, R2_TAG_SET, R1_TAG, R2_TAG

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
                        help="Workflow pipeline file to run (e.g., 'consensus-genome') Defaults to the main pipeline.",
                        default=None)
    parser.add_argument("-c", "--config_file", default=None, choices=["local", "cluster", "cluster_submit"])
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run")
    parser.add_argument("extra_args", nargs=argparse.REMAINDER, help="Additional arguments to pass to Snakemake")

    return parser

def run_pipeline(project_root, log_path, config_dict, config_path=None, pipeline_name=None, dry_run=False, extra_args=None, **kwargs):
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
    cores = execution.get("cores", 1)
    jobs = execution.get("jobs", 1)
    latency_wait = execution.get("latency_wait", 30)
    dry_run = execution.get("dry_run", False)
    cluster_config = execution.get("cluster_config", None)

    cmd = [
        "snakemake",
        "--snakefile", snakefile,
        "--config", f"project_root={project_root}", f"log_path={log_path}",
        "--configfile", config_path
    ]

    if mode == "slurm":
        cmd.extend(["--executor", "slurm"])
        cmd.extend(["--jobs", str(jobs)])
        cmd.extend(["--cores", str(cores)])  # Cores for the main process
        cmd.extend(["--latency-wait", str(latency_wait)])
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


def acquire_fast_a_q_files(working_dir, file_base, fastq=True, delims={'_', '-'}, without_r1=True):
    """
    Search a given working directory for either FASTQ or FASTA files.
    The file_base is assumed to be the start of the file name itself, stripped
    of all enclosing directories.
    Will locate the R1 file first, and if it does, will try to find a matching R2 file.
    If it finds any file starting with file_base but cannot find an R1 tag, it will
    take the first such example it finds and return it as R1.

    :param working_dir: Directory in which to search.
    :param file_base:
    :param fastq: If True, search for FASTQ files, else search for FASTA files.
    :param delims: Sets of possible file name delimiters
    :return: Dict {R1: full_file_path, R2: full_file_path}, where {R1: full_file_path, R2: None} means single-ended.
    """

    working_abspath = os.path.abspath(working_dir)

    if fastq:
        ext_set = FASTQ_EXT_SET
    else:
        ext_set = FASTA_EXT_SET

    possible_files = [f for f in os.listdir(working_dir) if f.endswith(tuple(ext_set))]

    if not possible_files:
        return {R1_TAG, None, R2_TAG, None}

    file_delim = None
    r1_loc = None
    ext_used = None
    file_basename_no_ext = None
    r1_file = None
    for file in possible_files:

        file_basename = os.path.normpath(os.path.basename(file)).strip()

        for ext in ext_set:
            if file_basename.endswith(ext):
                ext_used = ext
                file_basename_no_ext = file_basename.replace(ext, '')
                break

        if file_basename_no_ext is None:
            continue

        if file_basename_no_ext.startswith(file_base):
            r1_file = file  # Assign this the r1_file whether or not it has an R1 tag

            for delim in delims:  # First, find an R1 tag in a file starting with file_base if possible
                file_base_cols = file_basename_no_ext.split(delim)

                for i in range(len(file_base_cols)):
                    col = file_base_cols[i]

                    if col in R1_TAG_SET:
                        r1_loc = i
                        file_delim = delim
                        break

                if r1_loc is not None:
                    break

        if r1_loc is not None:
            break

    if r1_file is None:  # No file found with the given base
        return {R1_TAG: None, R2_TAG: None}
    else:  # An R1 file has been found
        r1_full_path = os.path.join(working_abspath, r1_file)
        if r1_loc is None:  # R1 file lacks an R1 tag, assume its single ended
            if without_r1:  # Only return the r1 file here is without_r1 is True
                return {R1_TAG: r1_full_path, R2_TAG: None}
            else:
                return {R1_TAG: None , R2_TAG: None}
        else:  # Search for matching R2 file

            if file_delim is None or ext_used is None:  # Unlikely, but just in case
                return {R1_TAG: None, R2_TAG: None}
            file_base_cols = file_basename_no_ext.split(file_delim)
            file_base_cols[r1_loc] = R2_TAG

            r2_file_base = file_delim.join(file_base_cols)
            r2_file_basename = r2_file_base + ext_used
            r2_full_path = os.path.join(working_abspath, r2_file_basename)

            if os.path.isfile(r2_full_path):
                return {R1_TAG: r1_full_path, R2_TAG: r2_full_path}
            else:
                return {R1_TAG: r1_full_path, R2_TAG: None}
