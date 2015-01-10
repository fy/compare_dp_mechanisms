# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Get results using the exponential mechanism.

# <codecell>

import sys, os, time
from collections import deque
import argparse
import numpy as np
import imp

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

# <codecell>

utility_functions = import_anywhere('utility_functions', [SCRIPT_DIR])
from utility_functions import get_chisq_sensitivity

# <markdowncell>

# Parse command line arguments.

# <codecell>

parser = argparse.ArgumentParser(description="Get results based on the JS method.")
parser.add_argument("k", metavar="NUM_SNP", help="number of SNPs to output", type=int)
parser.add_argument("e", metavar="EPSILON", help="privacy budget epsilon", type=float)
parser.add_argument("n_case", help="number of cases", type=int)
parser.add_argument("n_control", help="number of controls", type=int)
parser.add_argument("infile", help="input file of chisquare statistics")
parser.add_argument("-s", help="sensitivity of the scoring function", default=1, type=int)
args = parser.parse_args()

if not os.path.isfile(args.infile):
    sys.exit("The follwoing file does not exist: {}".format(args.infile))

# <markdowncell>

# Setup data.

# <codecell>

name_score_tuples = []
with open(args.infile, 'r') as infile:
    # skip header line
    garbage = infile.readline()
    for line in infile:
        name, score = line.strip().split()
        name_score_tuples.append((name, float(score)))

indexed_snp_name_dict = dict(enumerate([name for name, ss in name_score_tuples]))
snp_scores =  np.array([ss for name, ss in name_score_tuples])

# <markdowncell>

# Perform Laplace mechanism.

# <codecell>

sensitivity = get_chisq_sensitivity(args.n_case, args.n_control)
scale = sensitivity * 2.0 * args.k / args.e
scores_perturbed = snp_scores + np.random.laplace(scale=scale, 
                                                  size=len(snp_scores))
results_indices = np.argsort(scores_perturbed)[::-1][:args.k]
results_names = map(indexed_snp_name_dict.get, results_indices)

# <codecell>

for nn in results_names:
    print nn

