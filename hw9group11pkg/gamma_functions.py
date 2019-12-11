import numpy as np
import scipy.stats as st
import scipy
import scipy.optimize
import warnings
import bebi103
import pandas as pd


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
       
	bs_reps = bebi103.draw_bs_reps_mle(mle_fun=gamma_mle, gen_fun=gen_gamma, data=conc_data, size=5000, n_jobs=3)
    
	alpha_vals = bs_reps[:,0]
	beta_vals = bs_reps[:,1]
	return [np.mean(alpha_vals), np.mean(beta_vals)], np.percentile(bs_reps, [2.5, 97.5], axis=0)

def get_gamma_conf_intervals(concentrations, dataframe, conc_column, time_column):    
	
	"""
	generate bootstrap replicates for gamma MLE estimation of catastrophe
	data and return a dataframe with all relevant information, along 
	with arrays containing confidence intervals and MLE averages.

		Parameters
		----------
		concentrations: list conaining all tubulin concentrations of
		interest
		dataframe: dataframe with catastrophe times under colunns labeled
		with tubulin concentratiosn
		conc_column: the label used for the column with tubulin conc.
		in the dataframe
		time_column: the label used for the column with catastrophe
		times in the dataframe
		
		Returns
		-------
		output : dataframe, array, array
		dataframe with all relevant info
		array with confidence interavls
		array with MLE estimates

	"""
	
	#Make lists for both the MLE values and the confidence intervals
	MLE_vals = []
	conf_ints = []

	#Loop through the different tubulin concentrations
	for conc in concentrations:

		#Grab the relevant data for each different concentration
		data = dataframe[dataframe[conc_column] == conc][time_column]

		#Convert this data into an array for ease of analysis
		array = np.array(data)

		#Use our Gamma MLE funtion to pull out both MLE values and confidence intervals for our data
		mle_val, conf_int = gamma_bs_MLE(array)

		#Add the MLE values and confidence intervals to respective lists
		MLE_vals.append(mle_val)
		conf_ints.append(conf_int)

	#Convert our MLE values to an array for ease of use
	MLE_array = np.array(MLE_vals)

	#Make a dataframe containing these values and the relevant concentrations of tubulin that correspond to the measurements
	MLE_df = pd.DataFrame(MLE_array, columns=['alpha MLE', 'beta MLE'])
	MLE_df['Tubulin Concentration (uM)'] = concentrations
    
	alpha_low, alpha_high, beta_low, beta_high = [], [], [], []

	#Examine each confidence interval separately
	for piece in conf_ints:
		
		#Grab the alphas
		alphas = piece[:,0]
		
		#grab the beta values
		betas = piece[:,1]
		
		#append the high and lower bounds into separate lists
		alpha_low.append(alphas[0])
		alpha_high.append(alphas[1])
		beta_low.append(betas[0])
		beta_high.append(betas[1])
    
	MLE_df['low alpha'] = alpha_low
	MLE_df['high alpha'] = alpha_high
	MLE_df['low beta'] = beta_low
	MLE_df['high beta'] = beta_high
    
	return MLE_df, conf_ints, MLE_array
