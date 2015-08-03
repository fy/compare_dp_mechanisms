# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# This function implements Johnson & Shmatikov (2013)'s LocSig function.

# <codecell>

import numpy as np

# <codecell>

def loc_sig(epsilon, kk, sensitivity, all_scores):
    """Return the index of the top kk SNPs.
    
    Algorithm:
    (1) Weight the sampling probability of each indexed element by e^score. 
    (2) Sample an index using the sampling probabilities.
    (3) Store the index and set its sampling probability to be 0.
    (4) Repeat (1) to (3) until kk unique indices have been obtained.

    Args:
        epsilon: Privacy budget.
        kk: The top k SNPs to output. 
        sensitivity: sensitivity of the scoring function.
        all_scores: The scores to significance (+) or insignificance (-). 

    Returns:
        A list indices of the top kk SNPs.
    """
    def get_sampling_weights(exponent_vec):
        max_exponent = np.max(exponent_vec)
        sampling_weights = np.array([0 if ss is None else np.exp(ss - max_exponent + 50) 
                                     for ss in exponent_vec])
        sampling_weights = sampling_weights / np.sum(sampling_weights)
        return sampling_weights
    
    ## get the exponents used to calculate the sampling weights
    exponent_vec = [None if ss is None else 1. * ss * (1. * epsilon / kk) / (2 * sensitivity)
                    for ss in all_scores]
    sampling_weights = get_sampling_weights(exponent_vec)
    loc_vec = []
    for ii in xrange(kk):
        loc = np.random.choice(np.arange(len(sampling_weights)), 
                               p=sampling_weights)
        loc_vec.append(loc)
        exponent_vec[loc] = None
        sampling_weights = get_sampling_weights(exponent_vec)
    return loc_vec

