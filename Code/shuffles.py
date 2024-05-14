import numpy as np
import random
from sklearn.preprocessing import binarize
from scipy.stats import norm

# ----------------------------------------------------------------------------------------
# SHUFFLES
# ----------------------------------------------------------------------------------------

def fixed_pn_fixed_kc(W):
    """ Fix pn out-degree and kc in-degree

    Args:
        W (_type_): connectivity matrix (PN x KC)

    Returns:
        _type_: shuffled binary matrix 
    """    
    # Binarize connectivity
    mat_W = binarize(np.array(W))

    n_PNs = mat_W.shape[0]
    n_KCs = mat_W.shape[1]

    # In-degree to the KCs
    indeg = np.sum(mat_W>0,0)

    # Connection probability of PNs
    cprobs = np.mean(mat_W>0,1) # conn probabilities for each PN
    cprobs = cprobs / np.sum(cprobs) # so that all connection probabilities add up to 1

    Wshuf = np.zeros([n_PNs,n_KCs])

    for kc in range(n_KCs):
        num_inputs = np.random.choice(indeg) # use observed in-degrees to kcs
        inds = np.random.choice(n_PNs,num_inputs,p=cprobs, replace=False) # pick incident PNs without replacement, using observed connection probabilities 
        Wshuf[inds, kc] = 1

    return Wshuf


def fixed_pn(W):
    """ Fix just conneciton probabilities

    Args:
        W (_type_): connectivity matrix (PN x KC)

    Returns:
        _type_: shuffled binary matrix 
    """    
    # Binarize connectivity
    mat_W = binarize(np.array(W))

    n_PNs = mat_W.shape[0]
    n_KCs = mat_W.shape[1]

    # In-degree to the KCs
    indeg = np.sum(mat_W>0,0)
    mindeg = norm.ppf(0.05, loc=np.mean(indeg), scale=np.std(indeg)) # 5% percentile assuming normal distribution
    maxdeg = norm.ppf(0.95, loc=np.mean(indeg), scale=np.std(indeg)) # 95% percentile assuming normal distribution

    # Connection probability of PNs
    cprobs = np.mean(mat_W>0,1) # conn probabilities for each PN
    cprobs = cprobs / np.sum(cprobs) # so that all connection probabilities add up to 1

    Wshuf = np.zeros([n_PNs,n_KCs])

    for kc in range(n_KCs):
        num_inputs = np.random.randint(mindeg, maxdeg) # determine number of inputs to this KC, random int in range of observed kc in-degrees
        inds = np.random.choice(n_PNs,num_inputs,p=cprobs, replace=False) # pick incident PNs without replacement, using observed connection probabilities 
        Wshuf[inds, kc] = 1

    return Wshuf

def fixed_kc(W):
    """ Fix just kc in-degree distribution

    Args:
        W (_type_): connectivity matrix (PN x KC)

    Returns:
        _type_: shuffled binary matrix
    """    
    # Binarize connectivity
    mat_W = binarize(np.array(W))

    n_PNs = mat_W.shape[0]
    n_KCs = mat_W.shape[1]

    # In-degree to the KCs
    indeg = np.sum(mat_W>0,0)

    Wshuf = np.zeros([n_PNs,n_KCs])

    for kc in range(n_KCs):
        num_inputs = np.random.choice(indeg) # determine number of inputs to this KC, random choice (with replacement) from observed in-degrees to KCs
        inds = np.random.choice(n_PNs,num_inputs,replace=False) # pick incident PNs without replacement, uniform probability across pns
        Wshuf[inds, kc] = 1

    return Wshuf


def shufmat_old(mat_W):
    M = mat_W.shape[0]
    N = mat_W.shape[1]

    indeg = np.sum(mat_W>0,1)
    
    cprobs = np.mean(mat_W>0,0)
    cprobs = cprobs / np.sum(cprobs)

    Wshuf = np.zeros([M,N])

    for mi in range(M):
        num_inputs = np.random.choice(indeg)
        inds = np.random.choice(N,num_inputs,p=cprobs, replace=False)
        Wshuf[mi,inds] = 1

    return Wshuf