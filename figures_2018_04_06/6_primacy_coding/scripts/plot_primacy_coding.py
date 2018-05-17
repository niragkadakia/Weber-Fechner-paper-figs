"""
Analysis of primacy coding estimations. 

Created by Nirag Kadakia at 10:16 05-10-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from plot_formats import fig_primacy_errors_num_active, \
							fig_primacy_active_moment, \
							fig_primacy_errors_time
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									load_activities, save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_signal_trace_from_file


def primacy_coding(data_flag, active_rate=10, 
					intensity_idxs_to_plot=range(0, 5)):
	
	binary_errors = load_binary_errors(data_flag)
	errors_nonzero = binary_errors['errors_nonzero']
	errors_zero = binary_errors['errors_zero']
	Yy = load_activities(data_flag)
	
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	iter_var_names = ['temporal_adaptation_rate', 'seed_dSs', 'Kk', 
						'signal_trace_multiplier']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	Mm = list_dict['params']['Mm']
	
	# Signal array has shape (timepoints, iter_vars_dims[3])
	signal_file = list_dict['fixed_vars']['signal_trace_file']
	signal_data = load_signal_trace_from_file(signal_file)
	signal_array = list_dict['iter_vars']['signal_trace_multiplier']
	offset = list_dict['fixed_vars']['signal_trace_offset']
	Tt = signal_data[:, 0] - signal_data[0, 0]
	
	# Separate plots for adaptation rate, Kk, and signal step intensity
	iter_array = [iter_vars_dims[0], iter_vars_dims[2], iter_vars_dims[3]]
	it = sp.nditer(sp.zeros(iter_array), flags=['multi_index'])
	
	while not it.finished:
	
		rate, Kk, intensity = it.multi_index
		signal = (offset + signal_data[:, 1])*signal_array[intensity]
		
		# Binary -- receptors active? Has shape (timepoints, signal IDs, Mm)
		active_array = (Yy[:, rate, :, Kk, intensity, :]
								> active_rate)
		
		# Num. of active receptors; has shape (timepoints, signal IDs)
		num_activated = sp.sum(active_array, axis=-1)
		errors = 1./100*errors_nonzero[:, rate, :, Kk, intensity]*\
					errors_zero[:, rate, :, Kk, intensity]
		
		# list of Mm lists; each holds errors when iM receptors are active
		error_lists = [[] for x in xrange(Mm)]
				
		# For each signal, append the errors to the Mm lists in error_lists
		for idSs in range(iter_vars_dims[1]):
			active_error_arr = sp.vstack((num_activated[:, idSs].T, 
											errors[:, idSs].T)).T
			for iM in range(Mm):
				errors_iM = errors[sp.where(active_error_arr == iM)[0], idSs]
				error_lists[iM] = error_lists[iM] + errors_iM.tolist()
		
		# Plot errors versus number of active receptors
		fig = fig_primacy_errors_num_active()
		errors_avg = sp.zeros(Mm)
		for iM in range(Mm):
			errors_avg[iM] = sp.mean(error_lists[iM])
		plt.plot(range(1, Mm + 1), errors_avg, marker='o', color='k', lw=2)
		title='rate=%s, Kk=%s, intensity=%s' % (rate, Kk, intensity)
		plt.xlim(0, Mm)
		save_fig('primacy_errors_vs_num_active_rate=%s_Kk=%s_intensity=%s' 
					% (rate, Kk, intensity), subdir=data_flag)
	
		it.iternext()
	
	
	# Plot time shift to 50% fidelity, as a function of stimulus strength
	iter_array = [iter_vars_dims[0], iter_vars_dims[2]]
	it = sp.nditer(sp.zeros(iter_array), flags=['multi_index'])
	cmap = plt.cm.Reds
	
	while not it.finished:
		
		time_to_fifty = sp.zeros(len(intensity_idxs_to_plot))
		max_signal = sp.zeros(len(intensity_idxs_to_plot))
		
		for iI, intensity in enumerate(intensity_idxs_to_plot):
			rate, Kk = it.multi_index
			signal = (offset + signal_data[:, 1])*signal_array[intensity]
			errors = 1./100*errors_nonzero[:, rate, :, Kk, intensity]*\
						errors_zero[:, rate, :, Kk, intensity]
			max_signal[iI] = max(signal)
			time_to_fifty[iI] = Tt[sp.where(errors > 50)[0][0]]
			
		fig = fig_primacy_active_moment()
		plt.plot(max_signal, -1e3*(time_to_fifty - max(time_to_fifty)), 
					lw=2, color='k')
		plt.xlim(0, 10.1)
		plt.ylim(0, 30)
		save_fig('primacy_errors_vs_time_rate=%s_Kk=%s_intensities=%s' 
					% (rate, Kk, intensity_idxs_to_plot), subdir=data_flag)
			
		it.iternext()
		
		
		#color_val = 0.5 + 0.4*intensity/(len(intensity_idxs_to_plot) - 1)
			#plt.plot(Tt, sp.mean(errors, axis=1), color=cmap(color_val), lw=2)
			#plt.xlim(2.35, 2.6)
			
		
if __name__ == '__main__':
	data_flag = get_flag()
	primacy_coding(data_flag)