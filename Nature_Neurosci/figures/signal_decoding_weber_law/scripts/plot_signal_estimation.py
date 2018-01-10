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
from save_load_data import save_signal_estimation_fig, \
							load_signal_decoding_weber_law, \
							load_aggregated_object_list
from figure_plot_formats import signal_estimation_subfigures


def plot_signal_estimation_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict()):
	
	# Define the plot indices
	assert len(data_flags) == 2, \
		"Need command line arguments to be two, alternating " \
		"non-Weber law and Weber law."
		
	# Ready the plotting window; colormaps; colors; signals to plot
	fig = signal_estimation_subfigures()
	cmaps = [cm.Blues, cm.Reds]
	
	mu_dSs_to_plot = 60
	seed_Kk2_to_plot = 10
	true_signal_lw = 1.0
	true_signal_color = 'black'
	
	# Save same plot in both Weber Law and non-Weber Law folders
	
	for Weber_idx, data_flag in enumerate(data_flags):
	
		list_dict = read_specs_file(data_flag)
		iter_vars = list_dict['iter_vars']
		Nn = list_dict['params']['Nn']
		iter_plot_var = iter_vars.keys()[axes_to_plot[0]]
		x_axis_var = iter_vars.keys()[axes_to_plot[1]]
		
		data = load_signal_decoding_weber_law(data_flag)
		
		# Blue for non-adapted; red for adapted
		cmap = cmaps[Weber_idx]
		
		# Load CS objects for single stimuli plotting
		iter_vars_dims = []
		for iter_var in iter_vars:
			iter_vars_dims.append(len(iter_vars[iter_var]))		
		print ('Loading object list for single stimulus plot...'),
		CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
		print ('...loaded.')
	
		# Plot the bar graphs for true signal and estimates, top and 
		#  bottom; order by signal strength
		true_signal = CS_object_array[mu_dSs_to_plot, seed_Kk2_to_plot].dSs
		est_signal = CS_object_array[mu_dSs_to_plot, seed_Kk2_to_plot].dSs_est
		sorted_idxs = sp.argsort(true_signal)[::-1]
		plt.bar(sp.arange(Nn), true_signal[sorted_idxs]*(-1)**Weber_idx, 
			lw=true_signal_lw, edgecolor=true_signal_color,
			zorder=100, width=1.0, fill=False)
		plt.bar(sp.arange(Nn), est_signal[sorted_idxs]*(-1)**Weber_idx, 
			color=cmap(0.7), zorder=2, width=1.0)
	
		save_signal_estimation_fig(fig, data_flag)
		plt.show()
		
		
if __name__ == '__main__':
	data_flags = get_flags()
	plot_signal_estimation_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict(normal_eps_tuning_width=15))