"""
Calculate success ratios for CS batch runs, now with 
runs that have 3rd iter_var is Kk, to see how accuracy
depends on Kk.

Created by Nirag Kadakia at 23:35 02-14-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
sys.path.append('../src')
from utils import get_flag
from load_specs import read_specs_file
import matplotlib.pyplot as plt
from figure_plot_formats import signal_decoding_vs_Kk_subfigures
from save_load_data import load_aggregated_object_list, \
							load_signal_decoding_weber_law, \
							save_decoding_accuracy_vs_Kk_fig


def plot_signal_decoding_vs_Kk(data_flag, nonzero_bounds=[0.75, 1.25], 
								zero_bound=1./10., 
								threshold_pct_nonzero=100.0, 
								threshold_pct_zero=100.0):

	

	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	Nn = list_dict['params']['Nn']
	
	# Decoding accuracy subfigures
	fig = signal_decoding_vs_Kk_subfigures()

	assert len(iter_vars) == 3, "Need 3 iter_vars"
	iter_plot_var = iter_vars.keys()[0]
	x_axis_var = iter_vars.keys()[1]
	Kk_axis_var = iter_vars.keys()[2]
	
	data = load_signal_decoding_weber_law(data_flag)
	successes = data['successes']
	
	# Plot successes, averaged over second axis of successes array
	avg_successes = sp.average(successes, axis=1)*100.0
	plt.pcolormesh(avg_successes.T, cmap=plt.cm.hot, rasterized=True, 
					vmin=-5, vmax=110)
	plt.colorbar()
	save_decoding_accuracy_vs_Kk_fig(fig, data_flag)
	
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_signal_decoding_vs_Kk(data_flag)