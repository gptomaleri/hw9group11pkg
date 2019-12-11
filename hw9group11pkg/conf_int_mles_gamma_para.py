import numpy as np
import bebi103
import pandas as pd


def gen_gamma(params, size, rg):
"""    
Generate samples from the gamma distribution.

    Parameters
    ----------
    params: parameterized alpha and beta
    size:s number of samples that are drawn.
    rg: numpy.random.Generator instance
    
    Returns
    -------
    output : numpy array
        Randon samples from the parameterized gamma distribution.
    """

    alpha, beta = params
    return rg.gamma(alpha, 1 / beta, size=size)


def gamma_bs_MLE(conc_data):
    
"""
Generate parametric bootstrap replicates of maximum likelihood estimator
from gamma distribution using bebi103.draw_bs_reps_mle.

    Parameters
    ----------
    conc_data: MLEs data used for generating parametric bootstrap replicates
    
    Returns
    -------
    output : numpy array
    Bootstrap replicates of MLEs.

"""
       
    bs_reps = bebi103.draw_bs_reps_mle(
    mle_fun=gamma_mle, 
    gen_fun=gen_gamma, 
    data=conc_data, 
    size=5000, 
    n_jobs=3,
    progress_bar=True,
)
    alpha_vals = bs_reps[:,0]
    beta_vals = bs_reps[:,1]

    return [np.mean(alpha_vals), np.mean(beta_vals)], np.percentile(
        bs_reps, [2.5, 97.5], axis=0
    )