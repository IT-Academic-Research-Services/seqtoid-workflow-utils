import sys
import subprocess


def run_pipeline(logger, project_root, pipeline_name=None, dry_run=False, **kwargs):
    """
        Run a CypherID workflow.
        :param pipeline_name: Name of the workflow file (e.g., 'workflow1.smk') or None for main Snakefile
        :param dry_run: If True, perform a dry run (-n flag)
        :param kwargs: Additional Snakemake CLI arguments
        """

    logger.info(f"Attempting to run pipeline: {pipeline_name}")

    snakefile = project_root / "Snakefile" if pipeline_name is None else project_root/ "workflows" / pipeline_name / "Snakefile"
    if not snakefile.exists():
        print(f"Error: Snakefile {snakefile} not found.")
        sys.exit(1)

    cmd = ["snakemake", "--snakefile", str(snakefile)]
    if dry_run:
        cmd.append("-n")  # Dry run
    for key, value in kwargs.items():
        cmd.extend([f"--{key}", str(value)])

    try:
        subprocess.run(cmd, shell=False, check=True)
        logger.info(f"Successfully ran {pipeline_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run {pipeline_name}: {str(e)}")