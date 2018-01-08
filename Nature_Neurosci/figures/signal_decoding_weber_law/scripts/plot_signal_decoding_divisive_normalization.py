"""
Plot estimation error of inferred signal in compressed sensing 
from error objects (.npz) generated by calculate_errors.py, 
averaged over 2nd variables. Allows to plot a batch of 
errors to compare. This script plots single plots only at a time, 
for the divisive normalization plots. The first half of the 
command line data flags are for Weber-Law (plotted in red), 
the second half for non-Weber law (plotted in blue).

Current plots generated with command line arguments:


Created by Nirag Kadakia at 21:42 01-07-2018
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
from save_load_data import save_decoding_accuracy_fig, \
							load_signal_decoding_weber_law, \
							load_aggregated_object_list
from figure_plot_formats import divisive_normalization_subfigures


def plot_signal_decoding_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict()):
	
	# Define the plot indices
	assert len(data_flags) % 2 == 0, \
		"Need command line arguments to be normalization_idxs*2, " \
		"alternating Weber law and non-Weber law."
	
	# Ready the plotting window; colormaps; colors; signals to plot
	cmaps = [cm.Reds, cm.Blues]
	shades = sp.linspace(0.3, 0.7, len(data_flags)/2)
	
	# Plot success error figures
	for data_flag_idx, data_flag in enumerate(data_flags):
		
		Weber_idx = int(data_flag_idx / (len(data_flags)/2))
		normalization_idx = int(data_flag_idx % (len(data_flags)/2))
		
		# Decoding accuracy subfigures
		fig = divisive_normalization_subfigures()
		
		# Blue for non-adapted; red for adapted
		cmap = cmaps[Weber_idx]
		
		# Shade darkens for successive plots
		shade = shades[normalization_idx]
			
		list_dict = read_specs_file(data_flag)
		iter_vars = list_dict['iter_vars']
		Nn = list_dict['params']['Nn']
		iter_plot_var = iter_vars.keys()[axes_to_plot[0]]
		x_axis_var = iter_vars.keys()[axes_to_plot[1]]
		
		data = load_signal_decoding_weber_law(data_flag)
		successes = data['successes']
		
		nAxes = len(successes.shape)
		if nAxes > 2:
			successes = project_tensor(successes, iter_vars, 
										projected_variable_components,
										axes_to_plot)

		# Switch axes if necessary
		if axes_to_plot[0] > axes_to_plot[1]:    
			successes = successes.T
				
		# Plot successes, averaged over second axis of successes array
		avg_successes = sp.average(successes, axis=1)*100.0
		plt.plot(iter_vars[iter_plot_var], avg_successes, 
					color=cmap(shade), zorder=normalization_idx, lw=3.0)
	
		# Save same plot in both Weber Law and non-Weber Law folders
		for Weber_idx in range(2):
			data_flag = data_flags[data_flag_idx]
			save_decoding_accuracy_fig(fig, data_flag)
	
		
if __name__ == '__main__':
	data_flags = get_flags()
	plot_signal_decoding_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict(normal_eps_tuning_width=15))