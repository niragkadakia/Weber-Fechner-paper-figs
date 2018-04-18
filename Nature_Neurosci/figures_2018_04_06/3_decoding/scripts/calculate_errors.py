"""
Calculate binary errors and success ratios for CS batch runs, now regarding 
each CS estimation as successful or not. 

Created by Nirag Kadakia at 11:00 04-06-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
sys.path.append('../../shared_src')
from save_load_figure_data import save_binary_errors, save_MSE_errors, \
									save_success_ratios, save_odor_ID_errors

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from analysis import binary_success, binary_errors, MSE_errors
from load_data import load_aggregated_object_list


def calculate_errors(data_flags, nonzero_bounds=[0.7, 1.3], zero_bound=1./10.,
						threshold_pct_nonzero=100.0, threshold_pct_zero=100.0):
	"""
	Calcualte successful odor estimations based on binary errors.
	"""
		
	for data_flag in data_flags:
	
		list_dict = read_specs_file(data_flag)
		iter_vars_dims = []
		for iter_var in list_dict['iter_vars']:
			iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
		
		print ('Loading object list...'),
		CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
		print ('...loaded.')

		# Data structures 
		binary_dict = dict()
		binary_dict['nonzero'] = sp.zeros(iter_vars_dims)
		binary_dict['zero'] = sp.zeros(iter_vars_dims)
		successes = sp.zeros(iter_vars_dims)
		
		# Calculate binary errors
		it = sp.nditer(sp.zeros(iter_vars_dims), flags=['multi_index'])	
		while not it.finished:
			
			binary_error = binary_errors(CS_object_array[it.multi_index], 
									nonzero_bounds=nonzero_bounds,
									zero_bound=zero_bound)
			binary_dict['nonzero'][it.multi_index] = \
						binary_error['errors_nonzero']
			binary_dict['zero'][it.multi_index] = binary_error['errors_zero']
			
			successes[it.multi_index] = binary_success(
						binary_error['errors_nonzero'], 
						binary_error['errors_zero'],
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
			it.iternext()
		
		# Save binary errors (for identify v intensity) and pct success
		save_binary_errors(binary_dict['nonzero'], binary_dict['zero'], 
							data_flag)
		save_success_ratios(successes, data_flag)
	
	
if __name__ == '__main__':
	data_flags = sys.argv[1:]
	calculate_errors(data_flags)