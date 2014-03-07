# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Utility functions.

# <codecell>

import numpy as np

# <codecell>

def get_chisq_sensitivity(NN_case, NN_control):
    """sensitivity for the chi-square statistic based on 2x3 genotype tables"""
    NN = NN_case + NN_control  # total number of subjects
    CC_max = max(NN_case, NN_control)
    CC_min = min(NN_case, NN_control)
    sensitivity = 1. * NN**2 / (CC_min * (CC_max + 1))  # sensitivity of chisq
    return sensitivity


# <codecell>

def get_allelic_test_sensitivity(NN_case, NN_control):
    """sensitivity for the chi-square statistic based on 2x2 allelic tables derived from 2x3 genotype tables"""
    def sensitivity_type_1(SS, RR):
        NN = SS + RR
        return 1.0 * 8 * NN**2 * SS / \
                (RR * (2 * SS + 3) * (2 * SS + 1))
    
    def sensitivity_type_2(SS, RR):
        NN = SS + RR
        return 1.0 * 4 * NN**2 * ((2 * RR**2 - 1) * (2 * SS - 1) - 1) / \
                (SS * RR * (2 * RR + 1) * (2 * RR - 1) * (2 * SS + 1))

    return np.max([sensitivity_type_1(NN_case, NN_control),
                   sensitivity_type_1(NN_control, NN_case),
                   sensitivity_type_2(NN_case, NN_control),
                   sensitivity_type_2(NN_control, NN_case)])

# <codecell>

def check_table_valid(input_table):
    """Make sure that the margins (row sums and column sums ) are all positive.
    Args:
        input_table: A 2x3 numpy matrix.
    """
    ## check zero margins
    rowsum = np.array(map(np.sum, input_table))
    colsum = np.array(map(np.sum, input_table.T))
    if np.any(rowsum == 0) or np.any(colsum == 0):
        return False
    else:
        return True

# <codecell>

def chisq_stat(input_table):
    """Calculate the Pearson's chi-square staitsitc.
    Args:
        input_table: A 2x3 numpy matrix.

    Returns:
        A tuple (chisquare_statistics, degree_of_freedom).
    """
    input_table = input_table.astype(float)
    rowsum = np.array(map(np.sum, input_table))
    colsum = np.array(map(np.sum, input_table.T))
    expected = np.outer(rowsum, colsum) / np.sum(rowsum)
    # df = (len([1 for rr in rowsum if rr > 0]) - 1) * \
    #     (len([1 for cc in colsum if cc > 0]) - 1)
    chisq = np.sum(np.array(input_table[expected > 0] -
                            expected[expected > 0]) ** 2 /
                   expected[expected > 0])
    # return (chisq, df)
    return chisq

# <codecell>

def chisq_gradient(input_table):
    """Return the changable part of the gradient of the chi-square staitsitc.

    Args:
        input_table: A 2x3 numpy matrix.

    Returns:
        A four-element tuple consisting of the partial derivatives based on the
        parametrization the chi-square statistic by (r0, r1, n0, n1). The
        full parametrization would be
        (r0, r1, r2, s0, s1, s2, n0, n1, n2), where ri + si = ni. The returned 
        value will be scaled down by N^2 / (R * S).
    """
    input_table = input_table.astype(float)
    colsum = np.array(map(np.sum, input_table.T))
    ## divide each cell by colsum
    fraction_table = input_table / colsum
    dy_dr0, dy_dr1 = [2 * fraction_table[0, ii] - 2 * fraction_table[0, 2] for
                      ii in [0, 1]]
    dy_dn0, dy_dn1 = [-fraction_table[0, ii] ** 2 + fraction_table[0, 2] ** 2 for
                      ii in [0, 1]]
    return (dy_dr0, dy_dr1, dy_dn0, dy_dn1)

