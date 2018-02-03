"""
Plot sample signals to illustrate estimations.

Created by Nirag Kadakia at 17:55 01-10-2017
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
from matplotlib import cm
from utils import get_flags, project_tensor, polygon_under_graph
from load_specs import read_specs_file
from save_load_data import save_signal_estimation_zeros_fig, \
							save_signal_estimation_nonzeros_fig, \
							load_signal_decoding_weber_law, \
							load_aggregated_object_list
from figure_plot_formats import signal_estimation_subfigures


def plot_signal_estimation_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict()):
	
	# Define the plot indices
	assert len(data_flags) == 2, \
		"Need command line arguments to be two, alternating " \
		"non-Weber law and Weber law."
			
	cmaps = [cm.Blues, cm.Reds]
	alphas = [1.0, 0.8]
	mu_dSs_to_plot = [99]#sp.arange(0, 100, 5)
	seed_Kk2_to_plot = sp.arange(0, 100, 5)
	true_signal_lw = 1.0
	true_signal_color = 'black'
	CS_object_array = []
	
	# Load both object arrays only once
	for Weber_idx, data_flag in enumerate(data_flags):
	
		list_dict = read_specs_file(data_flag)
		iter_vars = list_dict['iter_vars']
		Nn = list_dict['params']['Nn']
		Kk = list_dict['params']['Kk']
		
		data = load_signal_decoding_weber_law(data_flag)
		
		# Load CS objects for single stimuli plotting
		iter_vars_dims = []
		for iter_var in iter_vars:
			iter_vars_dims.append(len(iter_vars[iter_var]))		
		print ('Loading object list for single stimulus plot...'),
		CS_object_array.append(load_aggregated_object_list(iter_vars_dims, 
								data_flag))
		print ('...loaded.')
	
	# Nonzero components
	for dSs_idx in mu_dSs_to_plot:
		for Kk2_idx in seed_Kk2_to_plot:
			
			# Ready the plotting window; colormaps; colors; signals to plot
			fig = signal_estimation_subfigures(nonzero=True)
			for Weber_idx, data_flag in enumerate(data_flags):

				# Blue for non-adapted; red for adapted
				color = cmaps[Weber_idx](0.6)

				# Plot the bar graphs for true signal and estimates, top and 
				#  bottom; order by signal strength
				true_signal = CS_object_array[Weber_idx][dSs_idx, Kk2_idx].dSs
				est_signal = CS_object_array[Weber_idx][dSs_idx, Kk2_idx].dSs_est
				sorted_idxs = sp.argsort(true_signal)[::-1]
				plt.bar(sp.arange(Kk), true_signal[sorted_idxs][:Kk]*(-1)**Weber_idx, 
					lw=true_signal_lw, edgecolor=true_signal_color,
					zorder=100, width=1.0, fill=False)
				plt.bar(sp.arange(Kk), est_signal[sorted_idxs][:Kk]*(-1)**Weber_idx, 
					color=color, zorder=2, width=1.0)
				
			for data_flag in data_flags:
				save_signal_estimation_nonzeros_fig(fig, data_flag, dSs_idx, Kk2_idx)
			plt.close()
		
	# Zero components
	for dSs_idx in mu_dSs_to_plot:
		for Kk2_idx in seed_Kk2_to_plot:
			
			# Ready the plotting window; colormaps; colors; signals to plot
			fig = signal_estimation_subfigures(nonzero=False)
			for Weber_idx, data_flag in enumerate(data_flags):
			
				# Blue for non-adapted; red for adapted
				color = cmaps[Weber_idx](0.6)
				alpha = alphas[Weber_idx]
				
				# Ordering here is actually not necessary
				true_signal = CS_object_array[Weber_idx][dSs_idx, Kk2_idx].dSs
				est_signal = CS_object_array[Weber_idx][dSs_idx, Kk2_idx].dSs_est
				sorted_idxs = sp.argsort(true_signal)[::-1]
				est_signal_nonzeros = est_signal[sorted_idxs][Kk:]
				est_signal_sorted_idxs = sp.argsort(est_signal_nonzeros)
				plt.fill_between(sp.arange(Nn - Kk), 0, est_signal_nonzeros\
					[est_signal_sorted_idxs], color=color, alpha=alpha)
			
			for data_flag in data_flags:
				save_signal_estimation_zeros_fig(fig, data_flag, dSs_idx, Kk2_idx)
			plt.close()
				
				
if __name__ == '__main__':
	data_flags = get_flags()
	plot_signal_estimation_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict(normal_eps_tuning_width=15))