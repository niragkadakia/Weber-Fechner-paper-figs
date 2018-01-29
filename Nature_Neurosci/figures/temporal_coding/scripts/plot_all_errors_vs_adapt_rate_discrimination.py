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
							save_nonzero_pct_err_vs_adapt_rate
from figure_plot_formats import perf_vs_adapt_rate_subfigures


def plot_errors_vs_adapt_rate_discrimination(data_flag, 
								iter_var_idxs_to_plot=None, 
								whf_thresh=0.1, ylims=[0, 100]):
	
	
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
	
	adapt_rates =  list_dict['iter_vars']['temporal_adaptation_rate']\
						[iter_var_idxs_to_plot]
	Nn = list_dict['params']['Nn']
	Kk = list_dict['params']['Kk']
	data_beg = sp.zeros((len(iter_var_idxs_to_plot), len(data_flags)))
	data_end = sp.zeros((len(iter_var_idxs_to_plot), len(data_flags)))
	
	# Only need to load encounters and signal once
	signal = data['signal_2']
	binary_hits = 1.*(signal > whf_thresh)
	enc_on_idxs = sp.where(sp.diff(binary_hits) > 0)[0] + 1
	enc_off_idxs = sp.where(sp.diff(binary_hits) < 0)[0] 
	
	print data['Tt'][enc_on_idxs] - data['Tt'][0]
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	cmap_beg = plt.cm.Purples
	cmap_end = plt.cm.Greens
		
	for data_idx, data_flag in enumerate(data_flags):
		
		# Load data for each flag
		data = load_temporal_errors(data_flag)
		zero_errors = data['zero_errors']
		
		rates_enc_beg = []
		rates_enc_end = []	
		
		# Plot zero errors, averaging over encounters and average index (2nd index)
		for iter_var_idx in iter_var_idxs_to_plot:
			rates_enc_beg.append(sp.average((Nn - Kk) - zero_errors\
								[enc_on_idxs, iter_var_idx, :]/100.*(Nn - Kk)))
			rates_enc_end.append(sp.average((Nn - Kk) - zero_errors\
								[enc_off_idxs, iter_var_idx, :]/100.*(Nn - Kk)))
		
		data_beg[:, data_idx] = rates_enc_beg
		data_end[:, data_idx] = rates_enc_end
		
		
	X, Y = sp.meshgrid(sp.log(adapt_rates)/sp.log(10), range(len(data_flags)))
	ax.plot_wireframe(X.T, Y.T, data_beg, color=cmap_beg(0.7))
	ax.plot_wireframe(X.T, Y.T, data_end, color=cmap_end(0.7))
	plt.show()	
	#save_nonzero_pct_err_vs_adapt_rate(fig, data_flag, whf_thresh)

	
if __name__ == '__main__':
	data_flags = sys.argv[1:]
	plot_errors_vs_adapt_rate_discrimination(data_flags, 
					iter_var_idxs_to_plot=sp.arange(4, 12),
					ylims=[0, 25])