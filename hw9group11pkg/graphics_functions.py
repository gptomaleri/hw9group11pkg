import numpy as np
import bokeh
import bokeh_catplot
bokeh.io.output_notebook()
import scipy.stats as st
import pandas as pd
import holoviews as hv
hv.extension('bokeh')
import bebi103

def two_step_CDF_comparison(beta_1, beta_2, data, conc):
	
	"""Graphical comparison of an empirical ECDF compared to a 
	theoretical two-step poisson with beta values for the two poisson
	processes as inputs along with an array of experimental waiting 
	times, as well as a string stating the concentration.

	Parameters
	----------
	beta_1 : rate of the first poisson process.
	beta_2 : rate of the seccond poisson process.
	data : array of experimental arrival times for a given tubulin 
	concentration.
	conc : string that states the chosen concentration for analysis

	Returns
	-------
	output : interactive bokeh graph of two-step CDF compared to 
	empirical data
	"""
	
	t_theor = np.linspace(0, 2000, 200)

	cdf = (
		beta_1 * beta_2 / (beta_2 - beta_1)
		* (
			(1 - np.exp(-beta_1 * t_theor)) / beta_1
			- (1 - np.exp(-beta_2 * t_theor)) / beta_2
		)
	)

	p = bokeh_catplot.ecdf(
		data=data, 
		val=conc, 
		conf_int=True
	)

	p.line(t_theor, 
		   cdf, 
		   line_width=2, 
		   color="orange"
		  )
	bokeh.io.show(p)
	
	
def gamma_function_CDF_comparison(alpha, beta, data, conc):
	
	"""Graphical comparison of an empirical ECDF compared to a 
	theoretical gamma dist CDF with beta values for the two poisson
	processes as inputs along with an array of experimental waiting 
	times, as well as a string stating the concentration.

	Parameters
    ----------
	alpha : alpha value for the gamma distribution.
	beta : rate of the poisson process.
	data : array of experimental arrival times for a given tubulin 
	concentration.
	conc : string that states the chosen concentration for analysis

	Returns
	-------
	output : interactive bokeh graph of gamma CDF compared to 
	experimental data.
	"""
	
	t_theor = np.linspace(0, 2000, 200)

	p = bokeh_catplot.ecdf(data, val=conc, conf_int=True)

	cdf = st.gamma.cdf(t_theor, alpha, loc=0, scale=1/beta)
	p.line(t_theor, cdf, line_width=2, color='orange')

	bokeh.io.show(p)
	
def plot_alpha_beta(dataframe):
    
	"""Graphical representation of alpha and beta values from gamma
	MLE optimization. The input is a dataframe containing a columns of
	'alpha MLE', 'beta MLE', and 'Tubulin Concentration (uM)'.

	Parameters
    ----------
	dataframe : dataframe containing information about the alpha and 
	beta MLE vals, along with a column of the tubulin conc.

	Returns
	-------
	output : Holoviews plot of alpha vs. beta 
	"""
    
    plot = hv.Points(
        data=dataframe,
        kdims=['alpha MLE', 'beta MLE'],
        vdims=['Tubulin Concentration (uM)'],
    ).groupby(
        'Tubulin Concentration (uM)'
    ).opts(
        title='beta MLE vs. alpha MLE',
        height=300,
        width=500,
        xlabel='α MLE',
        ylabel='β MLE',
    ).overlay()
    return plot


def show_alpha_beta_with_conf(conf_ints, MLE_centers):
	
	"""Graphical representation of alpha and beta values from gamma
	MLE optimization with condfidence intervals. The input is arrays
	containing the MLE confidence intervals and MLE centers.

	Parameters
    ----------
	conf_ints : array containing confidence intervals
	MLE_centers : array containign MLE centers for alpha and beta

	Returns
	-------
	output : two plots comparing alpha and beta values for all 
	concentrations 
	"""
	
	#Turn the list of confidence intervals into an array
	conf_array = np.array(conf_ints)

	#Make a plot of the alpha values with their respective confidence intervals
	p = bebi103.viz.plot_with_error_bars(
		centers=MLE_centers[:,0],
		confs=conf_array[:,:,0],
		names=["α_7", "α_9", "α_10", "α_12", "α_14"],
		x_axis_label="Events (α)",
	)

	bokeh.io.show(p)

	#Make a plot of the beta values with their respective confidence intervals
	q = bebi103.viz.plot_with_error_bars(
		centers=MLE_centers[:,1],
		confs=conf_array[:,:,1],
		names=["β_7", "β_9", "β_10", "β_12", "β_14"],
		x_axis_label="Arrival rate s-1 (β)",
	)

	bokeh.io.show(q)
