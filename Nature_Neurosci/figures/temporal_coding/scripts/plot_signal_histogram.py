"""
Plot histogram of temporal data (under construction)

Created by Nirag Kadakia at 17:00 01-24-2018
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


def plot_signal_histogram(data_flag):
	"""
	UNDER CONSTRUCTION
	"""
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	# Get time axes, select colors from blue (non-adaptive) to red (adapted)
	Tt = data['Tt'] - data['Tt'][0]
	
	bins = sp.linspace(-2, -0.2, 50)
	
	#signal_hist = sp.histogram(data['signal'], bins=bins)
	hist = sp.histogram(sp.log(data['signal'])/sp.log(10), bins=bins, normed=True)
	hist_bins = (hist[1][1:] + hist[1][:-1])/2.0
	hist_dx = hist[1][1] - hist[1][0]
	hist_vals = hist[0]
	
	cum_hist = sp.cumsum(hist_vals)*hist_dx
	
	# Plot signal
	#fig = signal_trace_subfigures(xlims)
	plt.plot(10**hist_bins, hist_vals)
	#plt.xlim(-2, 0)
	plt.xscale('log')
	#plt.yscale('log')
	plt.plot(10**hist_bins, cum_hist)
	plt.show()
	

if __name__ == '__main__':
	data_flags = sys.argv[1]
	plot_signal_histogram(data_flags)
	