#!/usr/bin/env python3
"""
CZID Phylo Pipeline Comparison
"""

import os
import glob
import pandas as pd
import hashlib
from typing import List, Tuple
import csv
from collections import defaultdict
from ete3 import Tree

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES = [
    'sample',
]

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'no diffs':     'identical',
    'identical':    'T',
    'differ':       'F'
}




# ────────────────────────────────────────────────────────────────
# Step functions
# ────────────────────────────────────────────────────────────────

czid_tree_path = os.path.join(CZID_DIR, "phylo-dbl-2_phylotree.nwk")
seqtoid_tree_path = os.path.join(SEQTOID_DIR, "validation-01-07-2026-phylo-2_phylotree.nwk")


t1 = Tree(czid_tree_path)   # or Tree("((A,B),C);")
t2 = Tree(seqtoid_tree_path)

# Robinson-Foulds distance
rf, max_rf, common_leaves, parts_t1, parts_t2, discarded_t1, discarded_t2 = t1.robinson_foulds(
    t2,
    unrooted_trees=True
)

print(f"RF distance: {rf} / {max_rf}")
print(f"Normalized RF: {rf / max_rf if max_rf > 0 else 0:.3f}")
print(f"Common leaves: {len(common_leaves)}")
print("Splits only in tree1:", len(parts_t1 - parts_t2))
print("Splits only in tree2:", len(parts_t2 - parts_t1))

