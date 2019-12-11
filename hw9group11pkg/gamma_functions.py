import numpy as np
import scipy.stats as st
import scipy
import scipy.optimize
import warnings


def log_like_gamma(params, t):
    
    """Log likelihood for a gamma distribution with
    input being alpha and beta parameters as well as counts.

    Parameters
    ----------
    params : array of the parameters alpha and beta.
    t : array array of counts.

    Returns
    -------
    output : float
        Log-likelihood.
    """
    
    alpha, beta = params
    
    if alpha <= 0 or beta <= 0:
        return -np.inf
    
    return st.gamma.logpdf(t, alpha, loc=0, scale=1/beta).sum()


def gamma_mle(t):
    
    """
    Perform maximum likelihood estimates for parameters for
    gamma measurements, parametrized by alpha and beta.
    
    Parameters
    ----------
    t : array array of counts.

    Returns
    -------
    output : float
         maximum likelihood estimates for alpha and beta.
    """
    
    
    # Initial guess
    t_bar = np.mean(t)
    beta_guess = t_bar / np.var(t)
    alpha_guess = t_bar * beta_guess

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            lambda params, t: -log_like_gamma(params, t),
            (alpha_guess, beta_guess),
            args=(t,),
            method="powell",
        )

    if res.success:
        return res.x
    else:
        raise RuntimeError('Convergence failed with message', res.message)