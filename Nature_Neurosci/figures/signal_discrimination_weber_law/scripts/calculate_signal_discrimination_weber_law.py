"""
Calculate success ratios for CS batch runs, now regarding each CS
estimation as successful or not. 

Created by Nirag Kadakia at 11:00 10-06-2017
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
from analysis import binary_success, binary_errors, binary_errors_dual_odor
from save_load_data import load_aggregated_object_list, \
							save_signal_discrimination_weber_law


def calculate_signal_discrimination_weber_law(data_flags, 
										nonzero_bounds=[0.8, 1.2], 
										zero_bound=1./10., 
										threshold_pct_nonzero=75.0, 
										threshold_pct_zero=75.0):
	"""
	"""
	
	for data_flag in data_flags:
	
		list_dict = read_specs_file(data_flag)
		for key in list_dict:
			exec("%s = list_dict[key]" % key)

		iter_vars_dims = []
		for iter_var in iter_vars:
			iter_vars_dims.append(len(iter_vars[iter_var]))		
		
		print ('Loading object list...'),
		CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
		print ('...loaded.')

		#Code not written for Kk_split = 0
		assert CS_object_array[0, 0].Kk_split != 0, "Need nonzero Kk_split"
		
		# Data structures 
		errors_zero = sp.zeros(iter_vars_dims)
		errors_nonzero_2 = sp.zeros(iter_vars_dims)
		errors_nonzero = sp.zeros(iter_vars_dims)
		successes_2 = sp.zeros(iter_vars_dims)
		successes = sp.zeros(iter_vars_dims)
		
		# Calculate binary errors
		it = sp.nditer(sp.zeros(iter_vars_dims), flags=['multi_index'])	
		while not it.finished:
			errors = binary_errors_dual_odor(CS_object_array[it.multi_index], 
									nonzero_bounds=nonzero_bounds,
									zero_bound=zero_bound)
			
			errors_nonzero[it.multi_index] = errors['errors_nonzero']
			errors_nonzero_2[it.multi_index] = errors['errors_nonzero_2']
			errors_zero[it.multi_index] = errors['errors_zero']
			it.iternext()
		
		# Calculate success ratios from binary errors
		it = sp.nditer(sp.zeros(iter_vars_dims), flags = ['multi_index'])
		while not it.finished:
			successes_2[it.multi_index] = binary_success(
						errors_nonzero_2[it.multi_index], 
						errors_zero[it.multi_index], 
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
			successes[it.multi_index] = binary_success(
						errors_nonzero[it.multi_index], 
						errors_zero[it.multi_index], 
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
			it.iternext()
			
		save_signal_discrimination_weber_law(successes, successes_2, data_flag)
	
if __name__ == '__main__':
	data_flags = sys.argv[1:]
	calculate_signal_discrimination_weber_law(data_flags)