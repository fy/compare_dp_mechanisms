# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# This function will calculate the distance to flipping the significance of a SNP.
# 
# Reference: Johnson & Shmatikov (2013)

# <markdowncell>

# Assumptions:
# 
# * The MAF of the controls are fixed. 

# <codecell>

import numpy as np
import scipy.stats as stats
import os
import imp

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

# <headingcell level=3>

# Some constants and utility functions

# <codecell>

DEBUG = True

# <codecell>

utility_functions = import_anywhere('utility_functions', [SCRIPT_DIR])
from utility_functions import chisq_stat

# <codecell>

## direction vectors for the case row and the genotype total of a 2x3 genotype table: (r_0, r_1, n_0, n_1)

def augment_direction_vector(uu):
    """Augment the direction vectors so that they are compatible with 2x3 input 
    tables.
    Args:
        uu: direction vector for the cases
    """
    uu_for_r = np.array([uu[0], uu[1], -(uu[0] + uu[1])])
    ## check if uu will change the first row
    if np.any(uu_for_r != 0):
        uu_for_s = np.zeros(3)
    else:
        uu_for_s = np.array([uu[2], uu[3], -(uu[2] + uu[3])])
    return np.array([uu_for_r, uu_for_s])

# <codecell>

def move_is_legal(input_table, uu):
    """Chcek whether moving the input table in the direction of uu is legal.
    Args:
        input_table: A 2x3 numpy matrix.
        uu: A direction vector compatible with the input table.

    Returns:
        True/False.
    """
    new_table = input_table + uu
    ## check negative cells
    if np.any(new_table < 0):
        return False
    ## check zero margins
    colsum = np.array(map(np.sum, new_table.T))
    if np.any(colsum == 0):
        return False
    return True

# <codecell>

def greedy_distance_to_significance_flip(input_table, threshold_pval):
    """Calculate the distance to flip the significance.

    Distance is defined as the Hamming distance in the space of all databases. 

    Args:
        input_table: A 2x3 numpy matrix.
        threshold_pval: Threshold p-value.

    Returns:
        The Hamming distance. Return "None" if no such distance can be found.
    """
    threshold_chisq = stats.chi2.isf(threshold_pval, 1) # !!! allelic chisquare
    ## If sig_direction > 0, make table significant.
    ## If sig_direction < 0, make table insignificant.
    sig_direction = 1 if chisq_stat(input_table) < threshold_chisq else -1
    ## if table is insignificant, then we want to make it significant
    if sig_direction >= 1:
        (H1, H2) = 0, 0
        # find H1: decrease r_0 then decrease r_1
        curr_table = input_table.copy()     ## make a copy of the input table
        while 1:
            # if r_0 can be decreased 
            if move_is_legal(curr_table, augment_direction_vector(np.array((-1, 0, -1, 0)))):
                next_move = np.array((-1, 0, -1, 0))
            # if r_0 cannot be decreased anymore
            else:
                # if r_1 can be decreased
                if move_is_legal(curr_table, augment_direction_vector(np.array((0, -1, 0, -1)))):
                    next_move = np.array((0, -1, 0, -1))
                else:
                    # if all things fail
                    break
            H1 += 1
            curr_table += augment_direction_vector(next_move)
            if chisq_stat(curr_table) >= threshold_chisq:
            ## Significance has flipped!
                 break
        # find H2: increase r_0 then increase r_0 and decrease r_1 at the same time
        curr_table = input_table.copy()    ## make a copy of the input table
        while 1:
            # if r_0 can be increased 
            if move_is_legal(curr_table, augment_direction_vector(np.array((1, 0, 1, 0)))):
                next_move = np.array((1, 0, 1, 0))
            # if r_0 cannot be increased anymore
            else:
                # if r_0 can be increased with r_1 being decreased
                if move_is_legal(curr_table, augment_direction_vector(np.array((1, -1, 1, -1)))):
                    next_move = np.array((1, -1, 1, -1))
                else:
                    # if all things fail
                    break
            H2 += 1
            curr_table += augment_direction_vector(next_move)
            if chisq_stat(curr_table) >= threshold_chisq:
            ## Significance has flipped!
                 break       
        return sig_direction * min([H1, H2])
    ## if table is significant, then we want to make it insignificant
    if sig_direction <= -1: 
        dist = 0
        # if table is to the right of red line (decreases as r_0 decreases)
        if (2 * input_table[0,0] + input_table[0,1]) * np.sum(input_table[1]) > \
            (2 * input_table[1,0] + input_table[1,1]) * np.sum(input_table[0]):
            curr_table = input_table.copy()    ## make a copy of the input table
            while 1:
                next_move = np.array((-1, 0, -1, 0))
                # if r_0 can be increased 
                if not move_is_legal(curr_table, augment_direction_vector(next_move)):
                    break
                curr_table += augment_direction_vector(next_move)
                dist += 1
                if chisq_stat(curr_table) < threshold_chisq:
                ## Significance has flipped!
                     break              
        # if table is to the right of red line (decreases as r_0 increases)
        else:
            curr_table = input_table.copy()    ## make a copy of the input table
            while 1:
                next_move = next_move = np.array((1, 0, 1, 0))
                # if r_0 can be increased 
                if not move_is_legal(curr_table, augment_direction_vector(next_move)):
                    break
                curr_table += augment_direction_vector(next_move)
                dist += 1
                if chisq_stat(curr_table) < threshold_chisq:
                ## Significance has flipped!
                     break      
        return dist * sig_direction
    return None

