#
# Code adapted from https://github.com/alitwinkumar/connectome_analyses
#
import numpy as np
import matplotlib.pyplot as plt

def shuf(J):
    """
    shuffles adjacency matrix J preserving row-wise in-degrees and column-wise connection probabilities
    J is an adjacency matrix with 0/1 entries. 
    returns: Jshuf, a binary shuffled weight matrix
    """
    assert np.in1d(J,[0,1]).all() #J must be a binary adjacency matrix

    M,N = J.shape
    cprobs = np.sum(J,0)/np.sum(J)

    indeg = np.sum(J,1) #in-degree sequence

    Jshuf = np.zeros([M,N])
    for ii in range(M):
        inds = np.random.choice(N,int(indeg[ii]),p=cprobs,replace=False)
        Jshuf[ii,inds] = 1

    return Jshuf

def compare_spectrum_shuf(J,shuf_func, label, Nshuf=200):
    """
    compares the squared singular value spectrum of J to that of shuffled matrices
    J is an adjacency matrix with 0/1 entries. shuf_func is a function that returns a random shuffle of J. Nshuf is the number of shuffles.
    """
    assert np.in1d(J,[0,1]).all() #J must be a binary adjacency matrix
    _,s,_ = np.linalg.svd(J) #singular values
    s = s**2/np.sum(s**2)

    R = len(s)
    sshuf = np.zeros([R,Nshuf]) #shuffled singular values
    for si in range(Nshuf):
        Jshuf = shuf_func(J)
        _,s_rand,_ = np.linalg.svd(Jshuf)
        s_rand = s_rand**2/np.sum(s_rand**2)
        sshuf[:,si] = s_rand

    m = np.mean(sshuf,1)
    #95% confidence intervals
    qmin = np.quantile(sshuf,0.05,axis=1)
    qmax = np.quantile(sshuf,0.95,axis=1) 

    # ax1.plot(range(n_components), lam_rand_mean, label=shuffle_key, marker='o', color='red')
    # ax1.errorbar(range(n_components), lam_rand_mean, errors, fmt='none', color='red')
    # ax1.plot(range(n_components), lam, label=dataset_title, marker='o', color='k')
    # ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.figure()
    plt.plot(1+np.arange(R),s,color="black",marker='o')
    plt.plot(1+np.arange(R),m,color="red",marker='o')
    plt.errorbar(1+np.arange(R),m,yerr=np.vstack([m-qmin,qmax-m]),color="red",ecolor="red")

    plt.ylabel('Fraction Variance Explained')
    plt.xlabel('Component')
    plt.legend(["data","shuffle"])
    plt.xlim(-1,15)


    plt.tight_layout()
    plt.savefig(label+'.svg')

    return s, m, qmin, qmax, sshuf

def compare_spectrum_shuf_normalized(J,shuf_func, label, Nshuf=200):
    """
    compares the squared singular value spectrum of J to that of shuffled matrices
    J is an adjacency matrix with 0/1 entries. shuf_func is a function that returns a random shuffle of J. Nshuf is the number of shuffles.
    """
    assert np.in1d(J,[0,1]).all() #J must be a binary adjacency matrix
    Jnorm = J - np.mean(J)
    _,s,_ = np.linalg.svd(Jnorm) #singular values
    s = s**2/np.sum(s**2)

    R = len(s)
    sshuf = np.zeros([R,Nshuf]) #shuffled singular values
    for si in range(Nshuf):
        Jshuf = shuf_func(J)
        Jshuf = Jshuf - np.mean(Jshuf)
        _,s_rand,_ = np.linalg.svd(Jshuf)
        s_rand = s_rand**2/np.sum(s_rand**2)
        sshuf[:,si] = s_rand

    m = np.mean(sshuf,1)
    #95% confidence intervals
    qmin = np.quantile(sshuf,0.05,axis=1)
    qmax = np.quantile(sshuf,0.95,axis=1) 

    # ax1.plot(range(n_components), lam_rand_mean, label=shuffle_key, marker='o', color='red')
    # ax1.errorbar(range(n_components), lam_rand_mean, errors, fmt='none', color='red')
    # ax1.plot(range(n_components), lam, label=dataset_title, marker='o', color='k')
    # ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.figure()
    plt.plot(1+np.arange(R),s,color="black",marker='o')
    plt.plot(1+np.arange(R),m,color="red",marker='o')
    plt.errorbar(1+np.arange(R),m,yerr=np.vstack([m-qmin,qmax-m]),color="red",ecolor="red")

    plt.ylabel('Fraction Variance Explained')
    plt.xlabel('Component')
    plt.legend(["data","shuffle"])
    plt.xlim(-1,15)
    plt.ylim(0, 1.0)

    plt.tight_layout()
    plt.savefig(label+'.svg')

    return s, m, qmin, qmax, sshuf
