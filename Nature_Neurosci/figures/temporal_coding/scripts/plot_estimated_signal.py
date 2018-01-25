"""
Plot estimated signals at a particular point on the time trace.

Created by Nirag Kadakia at 22:30 01-23-2018
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
from save_load_data import load_temporal_errors, save_est_signal_zeros_fig, \
							save_est_signal_nonzeros_fig
from figure_plot_formats import est_signal_zeros_subfigures, \
								est_signal_zeros_subfigures


def plot_estimated_signal(data_flag, iter_vars_idxs_to_plot=[0], 
						avg_var_idx_to_plot=0, dts_to_plot=[1.0, 1.5], 
						zero_ylims=[-1, 1], nonzero_ylims=[-1, 1]):
	"""
	Plot estimated signal zero and nonzero components at two distinct 
	times in a temporal coding trace.
	
	Args:
		data_flag: string; data identifier.
		iter_vars_idxs_to_plot: list; indices of the iterated variables to be
			plotted in successive plots. The iterated variables must be the 
			first variable of the iter_vars in specs.
		avg_var_idx_to_plot: int; index of the average variable  to be 
			plotted in successive plots. This variable must be the 
			second one of the iter_vars in specs.
		dts_to_plot: 2-entry list; times in the signal trace for which to plot.
			The nearest index of the time vector will be matched to this time.
		zero_ylims: 2-entry list; lower and upper limit of the signal 
			estimation plot for zero components. First value should be 
			negative, and second positive.
		nonzero_ylims: 2-entry list; lower and upper limit of the signal 
			estimation plot for nonzero components. First value should be 
			negative, and second positive.
		
	"""
	
	assert len(dts_to_plot) == 2, "dts_to_plot must be length-2 list"
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	
	# Get time axes and time indices at which to plot
	Tt = data['Tt'] - data['Tt'][0]
	iT_lo = (sp.absolute(Tt - dts_to_plot[0])).argmin()
	iT_hi = (sp.absolute(Tt - dts_to_plot[1])).argmin()
	
	for iter_vars_idx_to_plot in iter_vars_idxs_to_plot:
		
		fig = est_signal_zeros_subfigures(zero_ylims)
		
		Nn = list_dict['params']['Nn']
		Kk = list_dict['params']['Kk']
		range_lo = sp.arange(0, 2*list_dict['params']['Nn'], 2)
		range_hi = sp.arange(1, 2*Nn + 1, 2) - 0.2
		color_lo = plt.cm.Purples(0.85)
		color_hi = plt.cm.Greens(0.5)
		
		# Get indices of zero components of dSs; sort
		zero_idxs = (sp.where(data['dSs'][iT_lo, iter_vars_idx_to_plot, 
											avg_var_idx_to_plot, :] == 0))[0]
		est_signal_data_lo = data['dSs_est'][iT_lo, iter_vars_idx_to_plot, 
												avg_var_idx_to_plot, zero_idxs]
		sorted_idxs_zeros = sp.argsort(est_signal_data_lo)[::-1]
		
		# Get the estimated signals at the two times for the zero components
		est_signal_data_lo = est_signal_data_lo[sorted_idxs_zeros]
		est_signal_data_hi = data['dSs_est'][iT_hi, iter_vars_idx_to_plot, 
									avg_var_idx_to_plot, zero_idxs]
		est_signal_data_hi = est_signal_data_hi[sorted_idxs_zeros]
				
		# Plot as filled plot for easier visualization
		plt.fill_between(range_lo[Kk:], 0, est_signal_data_lo, 
							color=color_lo)
		plt.fill_between(range_hi[Kk:], 0, est_signal_data_hi, 
							color=color_hi, alpha=0.75)
		plt.xlim(2*Kk - 1, 2*Nn)
		
		save_est_signal_zeros_fig(fig, data_flag, dts_to_plot, 
									avg_var_idx_to_plot, 
									iter_vars_idx_to_plot)
		
		# Indicate the minimum value of the true signal
		true_signal_lo = sp.sort(data['dSs'][iT_lo, iter_vars_idx_to_plot, 
									avg_var_idx_to_plot, :])[-Kk:][::-1]
		true_signal_hi = sp.sort(data['dSs'][iT_hi, iter_vars_idx_to_plot, 
									avg_var_idx_to_plot, :])[-Kk:][::-1]
		print 'true signal nonzero component mins at two times: %.4f, %.4f'\
				% (true_signal_lo[-1], true_signal_hi[-1])
	
		# Plot nonzero components
		fig = est_signal_zeros_subfigures(nonzero_ylims)
		
		nonzero_idxs = (sp.where(data['dSs'][iT_lo, iter_vars_idx_to_plot, 
											avg_var_idx_to_plot, :] != 0))[0]
		sorted_idxs_nonzeros = sp.argsort(data['dSs'][iT_lo, iter_vars_idx_to_plot, 
									avg_var_idx_to_plot, :])[-Kk:][::-1]
		est_signal_data_lo = data['dSs_est'][iT_lo, iter_vars_idx_to_plot, 
												avg_var_idx_to_plot, 
												sorted_idxs_nonzeros]
		est_signal_data_hi = data['dSs_est'][iT_hi, iter_vars_idx_to_plot, 
												avg_var_idx_to_plot, 
												sorted_idxs_nonzeros]
		
		# Make bar graphs for nonzero components -- looks better for small #
		plt.bar(range_lo[-Kk:], est_signal_data_lo, color=color_lo)
		plt.bar(range_hi[-Kk:], est_signal_data_hi, color=color_hi, alpha=0.75)
		plt.bar(range_lo[-Kk:], true_signal_lo, color='black', 
					fill=False, alpha=0.75)
		plt.bar(range_hi[-Kk:], true_signal_hi, color='black', 
					fill=False, alpha=0.75)
		
		plt.xlim(2*Nn - 2*Kk - 1, 2*Nn)
		save_est_signal_nonzeros_fig(fig, data_flag, dts_to_plot, 
									avg_var_idx_to_plot, 
									iter_vars_idx_to_plot)
		
		
	
if __name__ == '__main__':
	data_flags = sys.argv[1]
	plot_estimated_signal(data_flags, iter_vars_idxs_to_plot=[7, 8, 9], 
							avg_var_idx_to_plot=7, dts_to_plot=[23.86, 25.17], 
							zero_ylims=[-0.01, 0.01], nonzero_ylims=[0, 0.12])