"""
Plot number of zero components erroneously estimated, versus adaptation rate,
at beginning and end of whiff encounters.

Created by Nirag Kadakia at 22:30 01-28-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
sys.path.append('../src')
import matplotlib.pyplot as plt
from load_specs import read_specs_file
from save_load_data import load_temporal_errors, \
							save_nonzero_pct_err_vs_adapt_rate
from figure_plot_formats import perf_vs_adapt_rate_subfigures


def plot_perf_vs_adapt_rate(data_flag, iter_var_idxs_to_plot=None, 
								whf_thresh=0.1, ylims=[0, 100]):
	"""
	"""
	
	# Load data and get iterated variables and their dimensions; plot vars
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))
	
	assert len(iter_vars_dims) == 2, "Need 2 iterated variables only"
	assert 'temporal_adaptation_rate' in list_dict['iter_vars'], \
		"'temporal_adaptation_rate' must be an iterated variable"
	
	if iter_var_idxs_to_plot is None:
		iter_var_idxs_to_plot = range(iter_vars_dims[0])
	else:
		pass
	avg_var_idxs = range(iter_vars_dims[1])
	
	# Relevant data needed
	adapt_rates =  list_dict['iter_vars']['temporal_adaptation_rate']\
					[iter_var_idxs_to_plot]
	signal = data['signal']
	Nn = list_dict['params']['Nn']
	Kk = list_dict['params']['Kk']
	
	# Get indices of beginning and end of whiff
	binary_hits = 1.*(signal > whf_thresh)
	enc_on_idxs = sp.where(sp.diff(binary_hits) > 0)[0] + 1
	enc_off_idxs = sp.where(sp.diff(binary_hits) < 0)[0] 
	rates_enc_beg = []
	rates_enc_end = []		
	
	# Figure properties
	color_lo = plt.cm.Purples(0.85)
	color_hi = plt.cm.Greens(0.5)
	
	# Plot zero errors, averaging over encounters and average index (2nd index)
	fig = perf_vs_adapt_rate_subfigures(ylims=ylims)
	for iter_var_idx in iter_var_idxs_to_plot:

		errors_beg = 0
		errors_end = 0
		for avg_var_idx in avg_var_idxs:
			
			# Calculate error value: use actual signal at on and off of whiff
			# Error values are mean absolute percentage error of the nonzero
			# components
			nonzero_idxs = sp.where(data['dSs'][0, iter_var_idx, 
									avg_var_idx, :] != 0)[0]
			dSs = data['dSs'][:, iter_var_idx, avg_var_idx, 
								nonzero_idxs][enc_on_idxs, :]
			dSs_est = data['dSs_est'][:, iter_var_idx, avg_var_idx, \
										nonzero_idxs][enc_on_idxs, :]
			errors_beg += sp.sum(sp.absolute((dSs - dSs_est)/dSs))/\
									len(enc_on_idxs)/Kk*100
			
			dSs = data['dSs'][:, iter_var_idx, avg_var_idx, 
								nonzero_idxs][enc_off_idxs]
			dSs_est = data['dSs_est'][:, iter_var_idx, avg_var_idx, \
										nonzero_idxs][enc_off_idxs, :]
			errors_end += sp.sum(sp.absolute((dSs - dSs_est)/dSs))/\
									len(enc_off_idxs)/Kk*100
			
		rates_enc_beg.append(errors_beg/len(avg_var_idxs))
		rates_enc_end.append(errors_end/len(avg_var_idxs))
	
	plt.ylim(0, 25)
	plt.plot(adapt_rates, rates_enc_beg, color=color_lo, lw=3)
	plt.plot(adapt_rates, rates_enc_end, color=color_hi, lw=3)
	
	save_nonzero_pct_err_vs_adapt_rate(fig, data_flag, whf_thresh)

	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_perf_vs_adapt_rate(data_flag, iter_var_idxs_to_plot=sp.arange(4, 12),
							ylims=[0, 25])