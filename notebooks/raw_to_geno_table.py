# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# This script will convert raw genotype data to genotype tables.

# <codecell>

import sys, os, time
import argparse
from collections import Counter
import imp
import numpy as np

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

# Parse arguments.

# <codecell>

parser = argparse.ArgumentParser()
parser.add_argument("case_file", help="raw case genotype data file")
parser.add_argument("control_file", help="raw control genotype data file")
parser.add_argument("outfile", help="output file")
args = parser.parse_args()

# <markdowncell>

# Process the raw data.

# <codecell>

utility_functions = import_anywhere('utility_functions', [SCRIPT_DIR])
from utility_functions import check_table_valid

# <codecell>

## read case data
case_data = {}
with open(args.case_file, 'r') as infile:
    for line in infile:
        fields = line.split()
        snp_name = fields[1]
        case_data[snp_name] = Counter(fields[4:])

## read control data
control_data = {}
with open(args.control_file, 'r') as infile:
    for line in infile:
        fields = line.split()
        snp_name = fields[1]
        control_data[snp_name] = Counter(fields[4:])

## combine case and control data and write as table   
with open(args.outfile, 'w') as outfile:
    headers = "name, case_0, case_1, case_2, ctrl_0, ctrl_1, ctrl_2".split(', ')
    line_template = '\t'.join(['{}'] * len(headers)) + '\n'
    outfile.write(line_template.format(*headers))
    for snp_name in case_data.keys():
        ## check whether the input table is a 2x3 contingency table with positive margins first
        input_table = np.array([[int(case_data[snp_name]['0']), 
                                 int(case_data[snp_name]['1']), 
                                 int(case_data[snp_name]['2'])], 
                                [int(control_data[snp_name]['0'],), 
                                 int(control_data[snp_name]['1'],),
                                 int(control_data[snp_name]['2'],)],])
        if not check_table_valid(input_table):
            continue
        ## the input table is valid. proceed.
        outfile.write(line_template.format(*[snp_name,
                                             case_data[snp_name]['0'],
                                             case_data[snp_name]['1'],
                                             case_data[snp_name]['2'],
                                             control_data[snp_name]['0'],
                                             control_data[snp_name]['1'],
                                             control_data[snp_name]['2'],]))

