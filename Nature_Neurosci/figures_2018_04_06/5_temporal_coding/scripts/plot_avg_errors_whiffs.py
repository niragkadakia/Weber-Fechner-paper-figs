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
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									save_fig
from plot_formats import fig_signal_trace

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects, \
						load_signal_trace_from_file


def plot_avg_errors(data_flags, rates_to_plot=[0, 1], whiff_threshold=8):
	"""
	"""

	avg_whiff_errors = sp.zeros((len(data_flags), 2))
	Kk_1_idx = 0
	Kk_2_idx = 2
	
	for iFlag, data_flag in enumerate(data_flags):
		
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
		
		# Clip array
		xlims = (0, 1.0)
		xlim_idxs = [int(len(Tt)*xlims[0]), int(len(Tt)*xlims[1])]
		plot_range = range(xlim_idxs[0], xlim_idxs[1])
		Tt = Tt[plot_range]
		Tt = sp.arange(len(Tt))
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
		
		# Get success rate during whiff times
		avg_success = sp.average(success, axis=2)*100.
		cum_success = 0
		num_pts = 0
		for nWhf in range(len(whiff_begs)):
			whf_range = range(whiff_begs[nWhf], whiff_ends[nWhf])
			cum_success += sp.sum(avg_success[whf_range, :], axis=0)
			num_pts += len(whf_range)
		whf_success = 1.*cum_success/num_pts
		avg_whiff_errors[iFlag, :] = whf_success

	for iFlag in range(len(data_flags)):
		x_range = [iFlag + 0.75, iFlag + 1.25]
		plt.bar(x_range, avg_whiff_errors[iFlag, :], width=0.5)
	plt.show()
	
	
if __name__ == '__main__':
	data_flags = get_flags()
	plot_avg_errors(data_flags)