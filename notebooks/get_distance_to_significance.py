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
from utility_functions import chisq_stat, chisq_gradient

# <codecell>

## direction vectors for the case row and the genotype total of a 2x3 genotype table: (r_0, r_1, n_0, n_1)
DIRECTION_VECTORS = np.array([
    (0, 0, -1, 1), (0, 0, -1, 0), (0, 0, 1, -1), (0, 0, 1, 0),
    (0, 0, 0, -1), (0, 0, 0, 1),
    (-1, 1, -1, 1), (-1, 0, -1, 0), (1, -1, 1, -1), (1, 0, 1, 0),
    (0, -1, 0, -1), (0, 1, 0, 1),
])

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

    At any point, say v, if the absolute value of the directional
    derivative is maximized in the direction u, then move v in the direction
    of u.
    """
    def find_best_legal_move(input_table, sig_direction):
        """Find the best legal move conditioning on the direction of the change
            in significance.
        """
        ## determine which moves are legal
        legal_moves = np.array([move_is_legal(input_table, augment_direction_vector(ee))
                                for ee in DIRECTION_VECTORS])
        if not np.any(legal_moves):
            ## return if not legal move is possible
            return None
        legal_move_indices = np.arange(len(DIRECTION_VECTORS))[legal_moves]
        ## get directional_derivatives
        gradient = chisq_gradient(input_table)
        directional_derivatives = np.array([
            np.dot(gradient, ee) for ee in np.array(DIRECTION_VECTORS)[legal_moves]])
        assert np.max(directional_derivatives * sig_direction) >= 0,\
        "The direction of the significance flip is %d " % (sig_direction) +\
        "but all of the directional derivatives are of the opposite sign." +\
        "The directional derivatives are: {}".format(str(directional_derivatives))
        if sig_direction > 0:
            best_move_idx = legal_move_indices[np.argmax(directional_derivatives)]
        else:
            best_move_idx = legal_move_indices[np.argmin(directional_derivatives)]
        if gradient == 0:
            print input_table
        return DIRECTION_VECTORS[best_move_idx]
    ## make a copy of the input table
    input_table = input_table.copy()
    dist = 0
    threshold_chisq = stats.chi2.isf(threshold_pval, 2)
    ## If sig_direction > 0, make table significant.
    ## If sig_direction < 0, make table insignificant.
    sig_direction = 1 if chisq_stat(input_table) < threshold_chisq else -1
    curr_table = input_table.copy()
    while 1:
        ## find best direction vector
        uu = find_best_legal_move(curr_table, sig_direction)
        if uu is None:
            break
        uu = augment_direction_vector(uu)
        curr_table += uu
        if DEBUG:
            print("chi-sqaure stat = {}, curr_table={}, uu={}".format(
                  chisq_stat(curr_table), 
                  str(curr_table.tolist()),
                  str(uu.tolist())))
        dist += 1
        if sig_direction * (chisq_stat(curr_table) - threshold_chisq) >= 0:
            ## Significance has flipped!
            return dist * sig_direction
    return None

# <codecell>

def main():
    print("Begin to debug.")
    ## check chisq_stat()
    tt1 = np.array([[16, 17, 18],
                    [10, 20, 5]])    
    assert round(chisq_stat(tt1), 3) == 6.214, "Error with chisq_stat()."
    ## check chisq_gradient()
    grad1 = chisq_gradient(tt1)
    assert map(round, grad1, [3] * len(grad1)) == [-0.334, -0.646, 0.234, 0.401],\
        "Error with chisq_gradient()."
    tt2 = np.array([[0, 17, 18],
                    [10, 20, 5]]) 
    grad2 = chisq_gradient(tt2)
    assert map(round, grad2, [3] * len(grad2)) == [-1.565, -0.646, 0.612, 0.401],\
        "Error with chisq_gradient()."
    ## check move_is_legal()
    assert move_is_legal(tt1, augment_direction_vector((1, -1, 1, -1))),\
        "Error with move_is_legal()."
    assert not move_is_legal(tt2, augment_direction_vector((-1, 0, -1, 0))),\
        "Error with move_is_legal()."
    print("Finish debugging.")

if __name__ == "__main__":
    main()

# <codecell>


