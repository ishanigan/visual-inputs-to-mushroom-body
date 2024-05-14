import numpy as np
import random
from sklearn.preprocessing import binarize
from sklearn.decomposition import PCA
import utils_general
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import shuffles
import scipy.linalg



# ----------------------------------------------------------------------------------------
# STRUCTURE ANALYSIS METHODS
# ----------------------------------------------------------------------------------------
def PCA_struct_analysis(dataset, shuffle_key, dataset_title='', n_components=20, n_shuffles=1000):

    # run PCA on original data
    W_ = np.copy(dataset)
    W_ = W_.T # need data to be in (KC, PN) format for PCA 
    W_ = binarize(np.array(W_), threshold=5)
    pca = PCA(n_components = n_components) 
    covar_matrix = pca.fit(W_)
    variances = covar_matrix.explained_variance_ratio_ # Var ratios

    # run PCA on random shuffles
    shuffle = getattr(shuffles, shuffle_key)
    variances_rand_list = []
    for i in range(n_shuffles):
        W_rand = np.copy(dataset)
        W_rand = shuffle(W_rand)
        W_rand = W_rand.T # need data to be in (KC, PN) format for PCA 
        pca_rand = PCA(n_components = n_components)
        covar_matrix_rand = pca_rand.fit(W_rand)
        variances_rand = covar_matrix_rand.explained_variance_ratio_ # Var ratios
        variances_rand_list.append(variances_rand)

    variances_rand_list = np.array(variances_rand_list) # (n, PC_components)
    variances_mean = np.mean(variances_rand_list, axis=0)

    # confidence intervals
    lowers = []
    highers = []
    for i in range(n_components):
        lower, higher = utils_general.confidence_interval(variances_rand_list[:, i])
        lowers.append(lower)
        highers.append(higher)
    errors = np.array([variances_mean-lowers, highers-variances_mean])

    # plotting
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_subplot(111)

    
    ax1.plot(range(n_components), variances_mean, label=shuffle_key, marker='o', color='red')
    ax1.errorbar(range(n_components), variances_mean, errors, fmt='none', color='red')
    ax1.plot(range(n_components), variances, label=dataset_title, marker='o', color='k')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    title = dataset_title + '_' + shuffle_key
    plt.ylabel('Fraction Variance Explained')
    plt.xlabel('Component')
    plt.title(title)


    plt.legend()
    plt.savefig(title+".pdf")
    plt.show()

    return variances, variances_rand_list

def PCA_struct_analysis_weighted(dataset, shuffle_key, dataset_title='', n_components=20, n_shuffles=1000):

    # run PCA on original data
    W_ = np.copy(dataset)
    W_ = W_.T # need data to be in (KC, PN) format for PCA 
    W_ = binarize(np.array(W_), threshold=5)
    pca = PCA(n_components = n_components) 
    covar_matrix = pca.fit(W_)
    variances = covar_matrix.explained_variance_ratio_ # Var ratios

    # run PCA on random shuffles
    shuffle = getattr(shuffles, shuffle_key)
    variances_rand_list = []
    for i in range(n_shuffles):
        W_rand = np.copy(dataset)
        W_rand = shuffle(W_rand)
        W_rand = W_rand.T # need data to be in (KC, PN) format for PCA 
        pca_rand = PCA(n_components = n_components)
        covar_matrix_rand = pca_rand.fit(W_rand)
        variances_rand = covar_matrix_rand.explained_variance_ratio_ # Var ratios
        variances_rand_list.append(variances_rand)

    variances_rand_list = np.array(variances_rand_list) # (n, PC_components)
    variances_mean = np.mean(variances_rand_list, axis=0)

    # confidence intervals
    lowers = []
    highers = []
    for i in range(n_components):
        lower, higher = utils_general.confidence_interval(variances_rand_list[:, i])
        lowers.append(lower)
        highers.append(higher)
    errors = np.array([variances_mean-lowers, highers-variances_mean])

    # plotting
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_subplot(111)

    ax1.plot(range(n_components), variances, label=dataset_title, marker='o')
    ax1.plot(range(n_components), variances_mean, label=shuffle_key, marker='o')
    ax1.errorbar(range(n_components), variances_mean, errors, fmt='none')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    title = dataset_title + '_' + shuffle_key
    plt.ylabel('Fraction Variance Explained')
    plt.xlabel('Component')
    plt.title(title)


    plt.legend()
    plt.savefig(title+".svg")
    plt.show()

    return variances, variances_rand_list

def PCA_struct_analysis_w_random(dataset, shuffle_key, dataset_title='', n_components=20, n_shuffles=1000):
    # input data is (KC, PN)

    # PCA on data
    U,S,V = scipy.linalg.svd(dataset)
    lam = S*S
    lam = lam/np.sum(lam)

    # PCA on shuffles 
    # run PCA on random shuffles
    shuffle = getattr(shuffles, shuffle_key)
    lam_rand_list = []
    for i in range(n_shuffles):
        W_rand = np.copy(dataset)
        W_rand = shuffle(W_rand.T)

        U,S_rand,V = scipy.linalg.svd(W_rand)
        lam_rand = S_rand*S_rand

        lam_rand = lam_rand/np.sum(lam_rand)
    
        lam_rand_list.append(lam_rand)

    lam_rand_list = np.array(lam_rand_list) # (n, PC_components)
    lam_rand_mean = np.mean(lam_rand_list, axis=0)

    # confidence intervals
    lowers = []
    highers = []
    for i in range(n_components):
        lower, higher = utils_general.confidence_interval(lam_rand_list[:, i])
        lowers.append(lower)
        highers.append(higher)
    errors = np.array([lam_rand_mean-lowers, highers-lam_rand_mean])

    # plotting
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_subplot(111)
    
    ax1.plot(range(n_components), lam_rand_mean, label=shuffle_key, marker='o', color='red')
    ax1.errorbar(range(n_components), lam_rand_mean, errors, fmt='none', color='red')
    ax1.plot(range(n_components), lam, label=dataset_title, marker='o', color='k')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    title = dataset_title + '_' + shuffle_key
    plt.xlim(-1,15)
    plt.ylabel('Fraction Variance Explained')
    plt.xlabel('Component')
    plt.title(title)

    plt.legend()
    plt.savefig(title+".svg")
    plt.show()

    return lam

