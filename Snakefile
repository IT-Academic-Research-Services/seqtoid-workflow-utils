"""
Snakefile for CypherID project.
"""

import os

project_root = os.path.dirname(os.path.abspath(workflow.snakefile))

include: "workflows/consensus-genome/Snakefile"


rule all:
    input:
        # Default target (optional, can be overridden by interactive selection)
        "data/output/consensus-genome.done",
