"""
Snakefile for CypherID project.
"""

import os

project_root = os.path.dirname(os.path.abspath(workflow.snakefile))

fail_message = "Placeholder error message."

# include: "workflows/consensus_genome/Snakefile"

rule all:
    input:
        "data/output/cypherid.done"

rule create_done_file:
    output:
        "data/output/cypherid.done"
    shell:
        """
        mkdir -p data/output
        echo "Done!" > {output}
        """