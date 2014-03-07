# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Get results using the JS method.

# <codecell>

import sys, os, time
from collections import deque, Counter
import argparse
import imp
import numpy as np

# <codecell>

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# <markdowncell>

# Utility functions.

# <codecell>

def import_anywhere(module_name, paths):
    """import methods from any folder"""
    try:
        f, filename, desc = imp.find_module(module_name, paths)
        return imp.load_module(module_name, f, filename, desc)
    finally:
        # Since we may exit via an exception, close fp explicitly.
        if f:
            f.close()

# <markdowncell>

# Parse command line arguments.

# <codecell>

parser = argparse.ArgumentParser(description="Get results based on the JS method.")
parser.add_argument("k", metavar="NUM_SNP", help="number of SNPs to output", type=int)
parser.add_argument("e", metavar="EPSILON", help="privacy budget epsilon", type=float)
parser.add_argument("infile", help="input file of JS distances")
parser.add_argument("-s", help="sensitivity of the scoring function", default=1, type=int)
args = parser.parse_args()

if not os.path.isfile(args.infile):
    sys.exit("The follwoing file does not exist: {}".format(args.infile))

# <markdowncell>

# Setup data.

# <codecell>

js_dist_tuples = []
with open(args.infile, 'r') as infile:
    # skip header line
    garbage = infile.readline()
    for line in infile:
        name, distance = line.strip().split()
        js_dist_tuples.append((name, int(distance)))

indexed_snp_name_dict = dict(enumerate([name for name, dd in js_dist_tuples]))
snp_scores =  np.array([-dd if dd >= 0 else -dd - 1 for name, dd in js_dist_tuples])

# <markdowncell>

# Perform JS algorithm.

# <codecell>

loc_sig = import_anywhere('loc_sig', [SCRIPT_DIR])
from loc_sig import loc_sig

# <codecell>

results_indices = loc_sig(args.e, args.k, args.s, snp_scores)
results_names = map(indexed_snp_name_dict.get, results_indices)

# <codecell>

for nn in results_names:
    print nn

