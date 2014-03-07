# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Write the JS distance to a file. Each line will have the following fields:
# 
# * name: SNP name
# * distance: JS distance

# <codecell>

import sys, os, time
from collections import deque, Counter
import argparse
import numpy as np
import scipy as sp
import imp

# <codecell>

parser = argparse.ArgumentParser(description="write JS distance")
parser.add_argument("infile", help="input genotype table file")
parser.add_argument("outfile", help="JS distance")
parser.add_argument("-p", help="threshold p-value", default=0.9, type=float)
args = parser.parse_args()

if not os.path.isfile(args.infile):
    sys.exit("The follwoing file does not exist: {}".format(args.infile))

# <markdowncell>

# Utility functions.

# <codecell>

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

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

# Import some functions.

# <codecell>

get_distance_to_significance = import_anywhere('get_distance_to_significance', [SCRIPT_DIR])
from get_distance_to_significance import greedy_distance_to_significance_flip
get_distance_to_significance.DEBUG = False

# <markdowncell>

# Set p-value.

# <codecell>

pval = args.p

# <markdowncell>

# Write JS distance.

# <codecell>

utility_functions = import_anywhere('utility_functions', [SCRIPT_DIR])
from utility_functions import check_table_valid

# <codecell>

start_time = time.time()

with open(args.infile, 'r') as infile, open(args.outfile, 'w') as outfile:
    outfile.write("{}\t{}\n".format(*['name', 'distance']))
    headers = infile.readline().split()
    for line in infile:
        dd = dict(zip(headers, line.split()))
        input_table = np.array([[int(dd['case_0']), 
                                 int(dd['case_1']),
                                 int(dd['case_2'])], 
                                [int(dd['ctrl_0']), 
                                 int(dd['ctrl_1']),
                                 int(dd['ctrl_2'])],])
        if not check_table_valid(input_table):
            continue
        print dd['name']
        dist = greedy_distance_to_significance_flip(input_table, pval)
        outfile.write('{}\t{}\n'.format(*[dd['name'], dist]))

print('Time spent: {} minutes.\n'.format(round((time.time() - start_time) / 60, 2)))

