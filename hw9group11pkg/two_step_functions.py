import numpy as np
import scipy.optimize
import warnings
import scipy.stats
import scipy

def two_step_log_like(params, t):
    
    """Log likelihood for a two successive Poisson processes with
    input being beta1 and beta2 parameters as well as counts.

    Parameters
    ----------
    params : array of the parameters beta1 and beta2.
    t : array array of counts.

    Returns
    -------
    output : float
        Log-likelihood.
    """
    
    
    beta_1, delta_beta = params

    n = len(t)

    #return negative infinity if delta_beta < 0, or if beta 1 is lower than beta1
    if delta_beta < 0 or beta_1 <= 0:
        return -np.inf
    
    #Use limiting gamma if difference between beta values very small
    if delta_beta < 1e-9:
        return 2 * n * np.log(beta_1) + np.sum(np.log(t)) - beta_1 * np.sum(t)

    
    out = n * (np.log(beta_1) + np.log(beta_1 + delta_beta) - np.log(delta_beta))
    out -= beta_1 * np.sum(t)
    
    return out + np.sum(np.log(1 - np.exp(-delta_beta * t)))
    


def two_step_mle(t):
    
    """
    Perform maximum likelihood estimates for parameters for
    a two successive Poisson processes, parametrized by beta1 and beta2.
    
    Parameters
    ----------
    t : array array of counts.

    Returns
    -------
    output : float
         maximum likelihood estimates for beta1 and beta2.
    """
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            lambda params, t: -two_step_log_like(params, t),
            np.array([0.005, 0.001]),
            args=(t,),
            method="powell"
        )
    
    if res.success:
        return res.x
    else:
        raise RuntimeError('Convergence failed with message', res.message)
        
        
