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


def plot_temporal_errors(data_flag, rates_to_plot = [0, 1], whiff_threshold=50):
	"""
	"""
	
	success = load_success_ratios(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	
	assert len(iter_vars) == 2, "Need 2 iter_vars"
	iter_var_names = ['temporal_adaptation_rate', 'seed_dSs']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	
	# Signal data can be loaded from specs file -- no need to open agg objs.
	signal_file = list_dict['fixed_vars']['signal_trace_file']
	signal_data = load_signal_trace_from_file(signal_file)
	Tt = signal_data[:, 0] - signal_data[0, 0]
	signal = signal_data[:, 1]
	
	# Plot signal alone
	fig = fig_signal_trace()
	plt.ylim(0, 1000)
	plt.xlim(Tt[0], Tt[-1])
	plt.plot(Tt, signal, lw = 3, color='k')
	save_fig('signal', subdir=data_flag)
	
	# Signal with whiffs highlighted
	whiff_hits_array = 1.*(signal > whiff_threshold)
	whiff_begs = Tt[sp.where(sp.diff(whiff_hits_array) == 1)[0]]
	whiff_ends = Tt[sp.where(sp.diff(whiff_hits_array) == -1)[0]]
		
	fig = fig_signal_trace()
	plt.ylim(0, 1000)
	plt.xlim(Tt[0], Tt[-1])
	plt.plot(Tt, signal, lw=3, color='k')
	for nWhf in range(len(whiff_begs)):
		plt.axvspan(whiff_begs[nWhf], whiff_ends[nWhf], alpha=0.2, color='r')
	save_fig('signal_with_whiffs', subdir=data_flag)
	
	adapt_rate_vals = iter_vars['temporal_adaptation_rate']
	avg_success = sp.average(success, axis=2)*100.
	std_success = sp.std(success, axis=2)*100.
	
	# Separate plots for each rate to visualize easier
	for iR, iRate in enumerate(rates_to_plot):
		fig = fig_signal_trace()
		color = 0.5 + 0.4*iR/(len(rates_to_plot) - 1)
		smoothed_avg = gaussian_filter(avg_success[:, iRate].T, sigma=2)
		plt.plot(Tt, smoothed_avg, color=plt.cm.Greys(color), lw=3)
		
		# Plot whiff regions in each estimation plot
		for nWhf in range(len(whiff_begs)):
			plt.axvspan(whiff_begs[nWhf], whiff_ends[nWhf], alpha=0.2, color='r')
		plt.xlim(Tt[0], Tt[-1])
		plt.yticks(sp.arange(0, 101, 50), fontsize=18)
		save_fig('errors_rate=%s' % iRate, subdir=data_flag)
		
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)