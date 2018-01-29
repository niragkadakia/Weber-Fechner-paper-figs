"""
Plot temporal trace errors and success ratios.

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
from save_load_data import load_temporal_errors, save_signal_trace_fig, \
							save_epsilon_trace_fig, save_zero_errors_trace_fig, \
							save_nonzero_errors_trace_fig
from figure_plot_formats import signal_trace_subfigures, \
								epsilon_trace_subfigures, \
								errors_trace_subfigures


def plot_temporal_data(data_flag, iter_var_axis=0, avg_var_axis=1, 
						iter_var_idxs_to_plot=None, avg_var_idxs_to_plot=None,
						xlims=[0, 1]):
	"""
	"""
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	# Set variables in iter_var_axis dimension to plot; verify in bounds
	if iter_var_idxs_to_plot is None:
		iter_var_idxs_to_plot = range(iter_vars_dims[iter_var_axis])
	assert max(iter_var_idxs_to_plot) <= iter_vars_dims[iter_var_axis] - 1, \
		"Plotting index %s is beyond iterated variable range" \
			% max(iter_var_idxs_to_plot)
	assert min(iter_var_idxs_to_plot) >= 0, "Plotting index must be >= 0"
	
	# Set variables in avg_var_axis dimension to average over; verify in bounds
	if avg_var_idxs_to_plot is None:
		avg_var_idxs_to_plot = range(iter_vars_dims[avg_var_axis])
	assert max(avg_var_idxs_to_plot) <= iter_vars_dims[avg_var_axis] - 1, \
		"Averaging index %s is beyond iterated variable range" \
			% max(avg_var_idxs_to_plot)
	assert min(avg_var_idxs_to_plot) >= 0, "Averaging index must be >= 0"
	
	# Get only the arrays in the avg_var_axis in the index array avg_var_idxs_to_plot
	avg_idxs_proj = sp.zeros(iter_vars_dims[avg_var_axis])
	avg_idxs_proj[avg_var_idxs_to_plot] = 1
	
	# Average over the avg_var_axis for given averaging indices
	full_array_shape = data['nonzero_errors'].shape
	for key in data.keys():
		if data[key].shape == full_array_shape:
			data[key] = sp.ndarray.compress(data[key], avg_idxs_proj, 
											axis=avg_var_axis + 1)
			data[key] = sp.average(data[key], axis=avg_var_axis + 1)
			assert len(data[key].shape) == 2, 'Cannot handle >2 iter_vars yet'
	
	# Get time axes, select colors from blue (non-adaptive) to red (adapted)
	Tt = data['Tt'] - data['Tt'][0]
	cmap = plt.cm.Reds
	colors = cmap(sp.linspace(0.2, 1, len(iter_var_idxs_to_plot)))
	lws = sp.linspace(3, 1.5, len(iter_var_idxs_to_plot))
	x_range_to_plot = range(int(xlims[0]*len(Tt)), int(xlims[1]*len(Tt)))
	
	# Plot signal
	fig = signal_trace_subfigures(xlims)
	plt.plot(Tt[x_range_to_plot], data['signal'][x_range_to_plot], 
				color='slategray', lw=2.5)
	plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
	save_signal_trace_fig(fig, data_flag, xlims)
	
	if 'Kk_split' in list_dict['params'].keys():
		if list_dict['params']['Kk_split'] > 0:
			fig = signal_trace_subfigures(xlims)
			plt.plot(Tt[x_range_to_plot], data['signal_2'][x_range_to_plot],
						color='slategray', lw=2.5)
			plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
			save_signal_trace_fig(fig, data_flag, xlims, dual=True)
		
	# Plot epsilons
	fig = epsilon_trace_subfigures(xlims)
	for iVar, iter_var_idx in enumerate(iter_var_idxs_to_plot):
		plt.plot(Tt[x_range_to_plot], data['avg_eps'][:, iter_var_idx]
				[x_range_to_plot], color=colors[iVar], lw=lws[iVar], alpha=0.8)
	plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
	save_epsilon_trace_fig(fig, data_flag, xlims)
	
	# Plot non-present sparse signal component errors
	fig = errors_trace_subfigures(xlims)
	for iVar, iter_var_idx in enumerate(iter_var_idxs_to_plot):
		plt.plot(Tt[x_range_to_plot], data['zero_errors'][:, iter_var_idx]
				[x_range_to_plot], color=colors[iVar], lw=lws[iVar], alpha=0.8)
	plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
	save_zero_errors_trace_fig(fig, data_flag, xlims)
	
	# Plot sparse component odorant errors
	fig = errors_trace_subfigures(xlims)
	for iVar, iter_var_idx in enumerate(iter_var_idxs_to_plot):
		plt.plot(Tt[x_range_to_plot], data['nonzero_errors'][:, iter_var_idx]
				[x_range_to_plot], color=colors[iVar], lw=lws[iVar], alpha=0.8)
	plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
	save_nonzero_errors_trace_fig(fig, data_flag, xlims)
	
	if 'Kk_split' in list_dict['params'].keys():
		if list_dict['params']['Kk_split'] > 0:
			fig = errors_trace_subfigures(xlims)
			for iVar, iter_var_idx in enumerate(iter_var_idxs_to_plot):
				plt.plot(Tt[x_range_to_plot], data['nonzero_errors_2']
						[:, iter_var_idx][x_range_to_plot], 
						color=colors[iVar], lw=lws[iVar], alpha=0.8)
			plt.xlim(Tt[0] + Tt[int(xlims[0]*len(Tt))], xlims[-1]*Tt[-1])
			save_nonzero_errors_trace_fig(fig, data_flag, xlims, dual=True)
				
	
if __name__ == '__main__':
	data_flags = sys.argv[1]
	plot_temporal_data(data_flags, iter_var_idxs_to_plot=[3, 5, 7, 8, 9, 10], 
						xlims=[.1, 0.14])