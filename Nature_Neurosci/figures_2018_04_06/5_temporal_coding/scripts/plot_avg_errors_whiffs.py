"""
Plot average errors over whiffs for different data flags.

Created by Nirag Kadakia at 10:40 05-06-2018
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
from save_load_figure_data import load_success_ratios, save_fig
from plot_formats import fig_avg_whiff_errors

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file
from load_data import load_signal_trace_from_file


def plot_avg_errors(data_flags, whiff_threshold=8):
	"""
	Will only plot 2 rates (indices 0 and 1 for temporal_adaptation_rate
	index in data file).
	"""

	# This should come from the data flags themselves
	background_intensities = [10, 50, 100, 300]
	
	for iFlag, data_flag in enumerate(data_flags):
		
		success = load_success_ratios(data_flag)
		list_dict = read_specs_file(data_flag)
		iter_vars_dims = []
		for iter_var in list_dict['iter_vars']:
			iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
		iter_vars = list_dict['iter_vars']
		
		iter_var_names = ['temporal_adaptation_rate', 'seed_dSs', 'Kk_1', 'Kk_2']
		for iName, name in enumerate(iter_var_names):
			assert iter_vars.keys()[iName] == name, "%sth variable "\
				"must have name %s" % (iName, name)
		
		try: 
			avg_whiff_errors
		except:
			error_array_dims = (len(data_flags), ) + (iter_vars_dims[0], 
								iter_vars_dims[2], iter_vars_dims[3])
			avg_whiff_errors = sp.zeros(error_array_dims)
			
		# Signal data can be loaded from specs file -- no need to open agg objs.
		signal_file = list_dict['fixed_vars']['signal_trace_file']
		signal_data = load_signal_trace_from_file(signal_file)
		Tt = sp.arange(len(signal_data))
		multiplier = list_dict['fixed_vars']['signal_trace_multiplier']
		offset = list_dict['fixed_vars']['signal_trace_offset']
		signal = (offset + signal_data[:, 1])*multiplier
		
		# Clip array to desired range
		xlims = (0.0, 0.4)
		xlim_idxs = [int(len(Tt)*xlims[0]), int(len(Tt)*xlims[1])]
		plot_range = range(xlim_idxs[0], xlim_idxs[1])
		Tt = Tt[plot_range]
		Tt = Tt - Tt[0]
		signal = signal[plot_range]
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
		
		avg_success = sp.average(success, axis=2)*100.
		cum_success = 0
		num_pts = 0
		
		# Grab the data
		for nWhf in range(len(whiff_begs)):
			whf_range = range(whiff_begs[nWhf], whiff_ends[nWhf])
			cum_success += sp.sum(avg_success[whf_range,...], axis=0)
			num_pts += len(whf_range)
		whf_success = 1.*cum_success/num_pts
		avg_whiff_errors[iFlag,...] = whf_success
	
	# x-axis of heatmap is background strength; y- is background complexity
	x_range = background_intensities
	y_range = range(iter_vars_dims[3])
	X, Y = sp.meshgrid(x_range, y_range)
	
	# Separate figure for each foreground/background complexity pair
	for Kk_1 in range(iter_vars_dims[2]):
		for Kk_2 in range(iter_vars_dims[2]):
			fig = fig_avg_whiff_errors()
			#plt.fill_between(x_range, 0, avg_whiff_errors[:, 1, Kk_1, Kk_2], 
			#	color=plt.cm.Greys(0.7))
			#plt.fill_between(x_range, 0, avg_whiff_errors[:, 0, Kk_1, Kk_2], 
			#	color=plt.cm.Greys(0.4))
			plt.plot(x_range, avg_whiff_errors[:, 1, Kk_1, Kk_2], 
						color=plt.cm.Greys(0.7), lw=4)
			plt.plot(x_range, avg_whiff_errors[:, 0, Kk_1, Kk_2], 
						color=plt.cm.Greys(0.4), lw=4)
			plt.yticks([0, 50, 100])
			plt.xscale('log')
			plt.ylim(0, 100)
			save_fig('%s_Kk_1=%s_Kk_2=%s' % (x_range, Kk_1, Kk_2), 
						subdir='avg_whiff_errors')
	
if __name__ == '__main__':
	data_flags = get_flags()
	plot_avg_errors(data_flags)