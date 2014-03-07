# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Write the $\chi^2$-statistics to a file. Each line will have the following fields:
# 
# * name: SNP name
# * score: $\chi^2$-statistics

# <codecell>

import sys, os, time
from collections import deque, Counter
import argparse
import numpy as np
import scipy as sp
import imp

# <codecell>

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# <codecell>

parser = argparse.ArgumentParser(description="write chisquare distance")
parser.add_argument("infile", help="input genotype table file")
parser.add_argument("outfile", help="chisquare statistics")
args = parser.parse_args()

if not os.path.isfile(args.infile):
    sys.exit("The follwoing file does not exist: {}".format(args.infile))

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

class Dummy(dict):
    pass

# <markdowncell>

# Write $\chi^2$-statistics

# <codecell>

utility_functions = import_anywhere('utility_functions', [SCRIPT_DIR])
from utility_functions import check_table_valid, chisq_stat

# <codecell>

start_time = time.time()

with open(args.infile, 'r') as infile, open(args.outfile, 'w') as outfile:
    outfile.write("{}\t{}\n".format(*['name', 'chisquare']))
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
        outfile.write('{}\t{}\n'.format(*[dd['name'], chisq_stat(input_table)]))

print('Time spent: {} minutes.\n'.format(round((time.time() - start_time) / 60, 2)))

