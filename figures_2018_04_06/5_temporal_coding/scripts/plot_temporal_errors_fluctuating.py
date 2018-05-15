"""
Plot estimation error in temporal coding tasks in a false color plot.

Created by Nirag Kadakia at 22:00 04-16-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									save_fig
from plot_formats import fig_signal_trace

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects, \
						load_signal_trace_from_file


def plot_temporal_errors(data_flag, rates_to_plot=[0, 1], whiff_threshold=8):

	# Which background and foreground complexities to plot
	Kk_1_idx = 0
	Kk_2_idx = 2
	
	success = load_success_ratios(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	
	iter_var_names = ['temporal_adaptation_rate', 'seed_dSs']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	
	# Set Kk_1 and Kk_2 if needed	
	if len(iter_vars) > 2:
		assert iter_vars.keys()[2] == 'Kk_1'
		assert iter_vars.keys()[3] == 'Kk_2'
		success = success[:, :, Kk_1_idx, Kk_2_idx]
	
	# Signal data can be loaded from specs file -- no need to open agg objs.
	signal_file = list_dict['fixed_vars']['signal_trace_file']
	signal_data = load_signal_trace_from_file(signal_file)
	Tt = signal_data[:, 0] - signal_data[0, 0]
	multiplier = list_dict['fixed_vars']['signal_trace_multiplier']
	offset = list_dict['fixed_vars']['signal_trace_offset']
	signal = (offset + signal_data[:, 1])*multiplier
	
	if 'signal_trace_file_2' in list_dict['fixed_vars']:
		signal_file_2 = list_dict['fixed_vars']['signal_trace_file_2']
		signal_data_2 = load_signal_trace_from_file(signal_file_2)
		multiplier = list_dict['fixed_vars']['signal_trace_multiplier_2']
		offset = list_dict['fixed_vars']['signal_trace_offset_2']
		signal_2 = (offset + signal_data_2[:, 1])*multiplier
	
	# Clip array to desired section
	xlims = (0.35, 0.40)
	xlim_idxs = [int(len(Tt)*xlims[0]), int(len(Tt)*xlims[1])]
	plot_range = range(xlim_idxs[0], xlim_idxs[1])
	Tt = Tt[plot_range]
	Tt = Tt - Tt[0]
	signal = signal[plot_range]
	if 'signal_trace_file_2' in list_dict['fixed_vars']:
		signal_2 = signal_2[plot_range]
	success = success[plot_range,...]
	
	# Signal with whiffs highlighted
	whiff_hits_array = 1.*(signal > whiff_threshold)
	whiff_begs = Tt[sp.where(sp.diff(whiff_hits_array) == 1)[0]]
	whiff_ends = Tt[sp.where(sp.diff(whiff_hits_array) == -1)[0]]
	if whiff_hits_array[0] > whiff_threshold:
		whiff_begs = sp.hstack((sp.zeros(1), whiff_begs))
	if len(whiff_begs) > len(whiff_ends):
		whiff_ends = sp.hstack((whiff_ends, Tt[-1]))
	if len(whiff_begs) < len(whiff_ends):
		whiff_begs = sp.hstack((Tt[0], whiff_begs))
	
	fig = fig_signal_trace()
	fig.set_size_inches(3.5, 2.5)
	plt.plot(Tt, signal, lw=2, color=plt.cm.Greens(0.9))
	if 'signal_trace_file_2' in list_dict['fixed_vars']:
		plt.plot(Tt, signal_2, lw=3, color=plt.cm.Blues(0.9), linestyle='-')
	for nWhf in range(len(whiff_begs)):
		plt.axvspan(whiff_begs[nWhf], whiff_ends[nWhf], alpha=0.15, color='r')
	plt.xticks(sp.arange(0, 100, 1), fontsize=18)
	plt.yticks(sp.arange(0, 1000, 100), fontsize=18)
	plt.xlim(Tt[0], Tt[-1])
	plt.ylim(0, 200)
	save_fig('signal_with_whiffs', subdir=data_flag)
	
	adapt_rate_vals = iter_vars['temporal_adaptation_rate']
	avg_success = sp.average(success, axis=2)*100.
	std_success = sp.std(success, axis=2)*100.
	
	# Separate plots for each rate to visualize easier
	for iR, iRate in enumerate(rates_to_plot):
		fig = fig_signal_trace()
		fig.set_size_inches(3.5, 2.5)
		cmap = plt.cm.Greys
		color = 0.6 + 0.4*iR/(len(rates_to_plot) - 1)
		smoothed_avg = gaussian_filter(avg_success[:, iRate].T, sigma=2)
		plt.plot(Tt, smoothed_avg, color=cmap(color), lw=2)
		
		# Plot whiff regions in each estimation plot
		for nWhf in range(len(whiff_begs)):
			plt.axvspan(whiff_begs[nWhf], whiff_ends[nWhf], alpha=0.15, 
						color='r')
		plt.xticks(sp.arange(0, 100, 1), fontsize=18)
		plt.yticks(sp.arange(0, 101, 50), fontsize=18)
		plt.xlim(Tt[0], Tt[-1])
		plt.ylim(0, 101)
		save_fig('errors_rates=%s' % iR, subdir=data_flag)
		
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)