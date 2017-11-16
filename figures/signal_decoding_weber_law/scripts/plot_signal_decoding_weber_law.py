"""
Plot estimation error of inferred signal in compressed sensing 
from error objects (.npz) generated by calculate_errors.py, 
averaged over 2nd variables. Allows to plot a batch of 
errors to compare.

Current plots generated with command line arguments:

mu_dSs_lo_Kk2_lo_diverse_Kk=7_WL 
mu_dSs_lo_Kk2_lo_diverse_Kk=7_no-WL 
mu_dSs_lo_Kk2_vryhi_diverse_Kk=7_WL 
mu_dSs_lo_Kk2_vryhi_diverse_Kk=7_no-WL 

Created by Nirag Kadakia at 22:00 09-21-2017
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
from matplotlib.collections import PolyCollection
from utils import get_flags, project_tensor, polygon_under_graph
from load_specs import read_specs_file
from save_load_data import save_signal_decoding_weber_law_fig, \
							load_signal_decoding_weber_law, \
							load_aggregated_object_list
from figure_plot_formats import single_decoding_weber_law_plot


def plot_signal_decoding_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict()):
	
	data_idxs = len(data_flags)
	assert data_idxs % 2 == 0, \
		"Need even number of command line arguments, alternating" \
		"Weber law and non-Weber law"
	
	# Ready the plotting window
	fig, ax = single_decoding_weber_law_plot(data_idxs=data_idxs)
	cmaps = [cm.Reds, cm.Blues]
	cmaps_r = [cm.Reds_r, cm.Blues_r]
	color_cycle_errors = sp.linspace(0.25, 0.7, data_idxs/2)
	linewidths_errors = sp.linspace(6.0, 3.0, data_idxs/2)
	
	# Plot for each command line argument
	for data_idx, data_flag in enumerate(data_flags):
	
		data_flag = str(data_flag)
		list_dict = read_specs_file(data_flag)
		for key in list_dict:
			exec("%s = list_dict[key]" % key)
			
		iter_plot_var = iter_vars.keys()[axes_to_plot[0]]
		x_axis_var = iter_vars.keys()[axes_to_plot[1]]
		Nn = list_dict['params']['Nn']
		
		data = load_signal_decoding_weber_law(data_flag)
		successes = data['successes']
		epsilons = data['epsilons']
		gains = data['gains']
		
		nAxes = len(successes.shape)
		if nAxes > 2:
			successes = project_tensor(successes, 
									iter_vars, projected_variable_components,
									axes_to_plot)
			
		# Switch axes if necessary
		if axes_to_plot[0] > axes_to_plot[1]:    
			successes = successes.T
		
		# Plot average errors
		average_successes = sp.average(successes, axis=1)*100.0
		cmap = cmaps[data_idx % 2]
		color = cmap(color_cycle_errors[data_idx / 2])
		linewidth = linewidths_errors[data_idx / 2]
		ax['successes'].plot(iter_vars[iter_plot_var], average_successes, 
					color=color, linewidth=linewidth)

		# Plot gains, averaged over odorant, and binned over different signals
		#   and over receptors
		if data_idx % 2 == 0:
			for iter_var_idx in range(len(iter_vars[iter_plot_var])):
				if iter_var_idx == 0:
					min = -8.0
					max = -1.0
					bins = sp.linspace(min, max, 1000)
					gains_data = sp.zeros((len(iter_vars[iter_plot_var]), 
											len(bins) - 1))
				odorant_averaged_gains = sp.average(gains[iter_var_idx], 
													axis=-1)
				hist, bin_edges = sp.histogram(
									sp.log(odorant_averaged_gains.flatten()) \
									/sp.log(10), bins=bins)
				gains_data[iter_var_idx, :] = hist 
			ax['gains_%s' % data_idx].pcolormesh(iter_vars[iter_plot_var],
										10.**bins, gains_data.T, cmap=cmaps_r[0],
										rasterized=True) 
		elif data_idx % 2 == 1:
			for iter_var_idx in range(len(iter_vars[iter_plot_var])):
				if iter_var_idx == 0:
					min = -2.0 
					max = 1.0
					bins = sp.linspace(min, max, 200)
					gains_data = sp.zeros((len(iter_vars[iter_plot_var]), 
											len(bins) - 1))
				odorant_averaged_gains = sp.average(gains[iter_var_idx], 
													axis=-1)
				hist, bin_edges = sp.histogram(
									sp.log(odorant_averaged_gains.flatten()) \
									/sp.log(10), bins=bins)
				gains_data[iter_var_idx, :] = hist 
			ax['gains_%s' % data_idx].pcolormesh(iter_vars[iter_plot_var],
										10.**bins, gains_data.T, cmap=cmaps_r[1], 
										rasterized=True)
			
	save_signal_decoding_weber_law_fig(fig)	
	
	# Plot estimated signals
	samples_to_plot = [0, 1]
	bkgrnds_to_plot = [40, 76]
	stimuli_to_plot = [1, 2]
	for data_idx, data_flag_idx in enumerate(samples_to_plot):
		data_flag = data_flags[data_flag_idx]
		iter_vars_dims = []
		for iter_var in iter_vars:
			iter_vars_dims.append(len(iter_vars[iter_var]))		

		print ('Loading object list...'),
		CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
		print ('...loaded.')
		
		cmap = cmaps[data_idx % 2]
		color = cmap(color_cycle_errors[-1])
		
		for bkgrnd_idx, bkgrnd_val in enumerate(bkgrnds_to_plot):
			ax['est_%s' % bkgrnd_idx].bar(
				sp.arange(Nn), CS_object_array[bkgrnd_val, 
				stimuli_to_plot[bkgrnd_idx]].dSs, edgecolor='black', 
				zorder=100, width=1.0, lw=1.0, fill=False)
			ax['est_%s' % bkgrnd_idx].bar(
				sp.arange(Nn), CS_object_array[bkgrnd_val, 
				stimuli_to_plot[bkgrnd_idx]].dSs_est, color=color, 
				zorder=2+data_idx, width=1.0)
		
		save_signal_decoding_weber_law_fig(fig)	
		
if __name__ == '__main__':
	data_flags = get_flags()
	plot_signal_decoding_weber_law(data_flags, axes_to_plot=[0, 1], 
				projected_variable_components=dict(normal_eps_tuning_width=15))