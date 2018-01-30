"""
Plot absolute percentage error of nonzero components of estiamted
signal as a function of adaptation rate, for discrimination tasks,
3rd axis of 3D plot is over complexity of background.

Created by Nirag Kadakia at 22:30 01-29-2018
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
from mpl_toolkits.mplot3d import Axes3D
from load_specs import read_specs_file
from save_load_data import load_temporal_errors, \
							save_adapt_rate_discr_num_wrong_zero_3d_fig, \
							save_adapt_rate_discr_pct_error_nonzero_3d_fig
from figure_plot_formats import adapt_rate_discrimination_3d_fig


def plot_errors_vs_adapt_rate_discrimination(data_flag, 
								iter_var_idxs_to_plot=None, 
								whf_threshs_1=[0.1, 100], 
								whf_threshs_2=[0.2, 100],
								zlims_zero=[0, 100], 
								zlims_nonzero=[0, 100]):
	
	
	# Load data and get iterated variables and their dimensions; plot vars
	data = load_temporal_errors(data_flags[0])
	list_dict = read_specs_file(data_flags[0])
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))
	
	assert len(iter_vars_dims) == 2, "Need 2 iterated variables only"
	assert 'temporal_adaptation_rate' in list_dict['iter_vars'], \
		"'temporal_adaptation_rate' must be an iterated variable"
	assert 'Kk_split' in list_dict['params'], "Need Kk_split in params"
	
	if iter_var_idxs_to_plot is None:
		iter_var_idxs_to_plot = range(iter_vars_dims[0])
	else:
		pass
	avg_var_idxs = range(iter_vars_dims[1])
	
	# Data structures and some variables needed for plotting
	adapt_rates =  list_dict['iter_vars']['temporal_adaptation_rate']\
						[iter_var_idxs_to_plot]
	Nn = list_dict['params']['Nn']
	Kk = list_dict['params']['Kk']
	data_rates = sp.zeros((len(iter_var_idxs_to_plot), len(data_flags), 2))
	log_adapt_rates = sp.log(adapt_rates)/sp.log(10)
	
	# Only need to load encounters and signal once
	signal_1 = data['signal']
	signal_2 = data['signal_2']
	binary_hits = (signal_2 < whf_threshs_2[1])*(signal_2 > whf_threshs_2[0])*\
					(signal_1 < whf_threshs_1[1])*(signal_1 > whf_threshs_1[0])
	
	cmap = plt.cm.viridis
	
	for data_idx, data_flag in enumerate(data_flags):
		
		# Load data for each flag
		data = load_temporal_errors(data_flag)
		zero_errors = []
		nonzero_errors = []
		list_dict = read_specs_file(data_flag)
		Kk_split = list_dict['params']['Kk_split']
	
		for iter_var_idx in iter_var_idxs_to_plot:
		
			# Use fluctuation in 2nd signal component 
			zero_errors.append(sp.average(100 - data['zero_errors_2']\
								[binary_hits, iter_var_idx, :])/100.*(Nn - Kk))
			
			# For nonzero components, only use nonzero components in 2nd signal
			nonzero_errors_avg_idxs = 0
			for avg_var_idx in avg_var_idxs:
				nonzero_idxs = data['idxs_2'][0, iter_var_idx, avg_var_idx, :]\
								.astype(int)
				dSs = data['dSs'][:, iter_var_idx, avg_var_idx, 
									nonzero_idxs][binary_hits, :]
				dSs_est = data['dSs_est'][:, iter_var_idx, avg_var_idx, \
											nonzero_idxs][binary_hits, :]
				nonzero_errors_avg_idxs += sp.sum(sp.absolute((dSs - \
										dSs_est)/dSs))/sp.sum(binary_hits)\
										/Kk_split*100
				
			nonzero_errors.append(nonzero_errors_avg_idxs/len(avg_var_idxs))
			
		data_rates[:, data_idx, 0] = zero_errors
		data_rates[:, data_idx, 1] = nonzero_errors
	
	fig, ax = adapt_rate_discrimination_3d_fig(data_rates[:, :, 0],
										adapt_rates, zlims=zlims_zero)
	X, Y = sp.meshgrid(log_adapt_rates, range(len(data_flags)))
	ax.plot_surface(X.T, Y.T, data_rates[:, :, 0], cmap=cmap, 
						lw=0, alpha=0.8, antialiased=False)
	save_adapt_rate_discr_num_wrong_zero_3d_fig(fig, data_flag, 
											whf_threshs_1, whf_threshs_2)
	
	fig, ax = adapt_rate_discrimination_3d_fig(data_rates[:, :, 1], 
										adapt_rates, zlims=zlims_nonzero)
	X, Y = sp.meshgrid(log_adapt_rates, range(len(data_flags)))
	ax.plot_surface(X.T, Y.T, data_rates[:, :, 1], cmap=cmap, 
						lw=0, alpha=0.5, antialiased=False)
	save_adapt_rate_discr_pct_error_nonzero_3d_fig(fig, data_flag, 
											whf_threshs_1, whf_threshs_2)
	
	
if __name__ == '__main__':
	data_flags = sys.argv[1:]
	plot_errors_vs_adapt_rate_discrimination(data_flags, 
					iter_var_idxs_to_plot=sp.arange(4, 12),
					zlims_zero=[0, 30], zlims_nonzero=[0, 25])