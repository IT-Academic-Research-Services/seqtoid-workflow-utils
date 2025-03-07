"""
Snakefile for CypherID project.
"""


include: "workflows/consensus-genome/Snakefile"


rule all:
    input:
        # Default target (optional, can be overridden by interactive selection)
        "data/output/consensus-genome.done",
