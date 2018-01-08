"""
Plot receptor activity distributions for multiple estimations.

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
from save_load_data import save_activities_fig, \
							load_signal_decoding_weber_law, \
							load_aggregated_object_list
from figure_plot_formats import activities_subfigures


def plot_activities(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict()):
	
	# Define the plot indices
	assert len(data_flags) % 2 == 0, \
		"Need command line arguments to be normalization_idxs*2, " \
		"alternating Weber law and non-Weber law."
	
	# Ready the plotting window; colormaps; colors; signals to plot
	cmaps = [cm.Reds, cm.Blues]
	cmaps_r = [cm.Reds_r, cm.Blues_r]
	
	# Plot success error figures
	for data_flag_idx, data_flag in enumerate(data_flags):
		
		# Blue or red depends on this index
		Weber_idx = data_flag_idx % 2
		
		list_dict = read_specs_file(data_flag)
		for key in list_dict:
			exec("%s = list_dict[key]" % key)
			
		iter_plot_var = iter_vars.keys()[axes_to_plot[0]]
		x_axis_var = iter_vars.keys()[axes_to_plot[1]]
		
		data = load_signal_decoding_weber_law(data_flag)
		activities = data['activities']
		
		for iter_var_idx in range(len(iter_vars[iter_plot_var])):
			if iter_var_idx == 0:
				bins = sp.linspace(-3, 0, 5000)
				activities_data = sp.zeros((len(iter_vars[iter_plot_var]), 
										len(bins) - 1))
				
			# Histogram takes all activities for all stimuli choices, for 
			#  this particular background mean
			hist, bin_edges = sp.histogram(
								sp.log(activities[iter_var_idx, :, :])\
								/sp.log(10), bins=100)
			
			# Interpolate to identical scale for all columns.
			interp_hist = sp.interp(bins, bin_edges[:-1], hist)[:-1]
			activities_data[iter_var_idx, :] = interp_hist
		
		fig = activities_subfigures()
		plt.pcolormesh(iter_vars[iter_plot_var], 10.**bins, activities_data.T, 
					cmap=cmaps_r[Weber_idx], rasterized=True) 
		save_activities_fig(fig, data_flag)
		
	
if __name__ == '__main__':
	data_flags = get_flags()
	plot_activities(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict(normal_eps_tuning_width=15))