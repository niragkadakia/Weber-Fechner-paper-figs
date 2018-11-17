"""
Plots of decoding accuracy versus number of active ORNs

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
from plot_formats import fig_primacy_errors_num_active
from save_load_figure_data import load_binary_errors, \
									load_activities, save_fig_no_whtspc

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file


def accuracy_vs_num_active(data_flag, active_rate=15):
	
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
		assert list(iter_vars.keys())[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	Mm = list_dict['params']['Mm']
	num_rates = iter_vars_dims[0]
	num_odors = iter_vars_dims[1]
	num_Kks = iter_vars_dims[2]
	num_intensities = iter_vars_dims[3]
	
	# Separate plots for adaptation rate, Kk, and signal step intensity
	iter_array_dims = [num_rates, num_Kks, num_intensities]
	it = sp.nditer(sp.zeros(iter_array_dims), flags=['multi_index'])
	
	while not it.finished:
	
		rate, Kk, intensity = it.multi_index
		
		# Are ORNs active (1) or not (0)? Shape = (timepoints, signal IDs, Mm)
		active_array = (Yy[:, rate, :, Kk, intensity, :] > active_rate)
		
		# Num. of active receptors; has shape (timepoints, signal IDs)
		num_activated = sp.sum(active_array, axis=-1)
		errors = (errors_nonzero[:, rate, :, Kk, intensity]*\
					errors_zero[:, rate, :, Kk, intensity]*1./100)
		
		# list of Mm lists; each holds errors when iM receptors are active
		error_lists = [[] for x in sp.arange(Mm)]
				
		# For each signal, append the errors to the Mm lists in error_lists
		for idSs in range(num_odors):
			active_error_arr = sp.vstack((num_activated[:, idSs].T, 
											errors[:, idSs].T)).T
			
			# Find errors at all times at which a given # (iM) of ORNs active
			for iM in range(Mm):
				errors_iM = errors[sp.where(active_error_arr == iM)[0], idSs]
				error_lists[iM] = error_lists[iM] + errors_iM.tolist()
		
		# Plot errors versus number of active receptors
		fig = fig_primacy_errors_num_active()
		errors_avg = sp.zeros(Mm)
		for iM in range(Mm):
			errors_avg[iM] = sp.mean(error_lists[iM])
		color_val = 0.85*Kk/(num_Kks - 1)
		plt.plot(range(1, Mm + 1), errors_avg, marker='o', lw=2,
					color=plt.cm.viridis(color_val))
		plt.xlim(0, Mm)
		save_fig_no_whtspc('primacy_errors_vs_num_active_rate=%s_Kk=%s_intensity=%s' 
					% (rate, Kk, intensity), subdir=data_flag, 
					no_ax=False)
	
		it.iternext()
	
		
if __name__ == '__main__':
	data_flag = get_flag()
	accuracy_vs_num_active(data_flag)