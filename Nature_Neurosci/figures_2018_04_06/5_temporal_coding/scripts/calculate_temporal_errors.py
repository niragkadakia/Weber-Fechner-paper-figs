"""
Calculate estimation error of inferred signal in compressed sensing 
decoding for a full signal trace in time. Success ratios are calculated
as a function of time. Saved to objects as a distinct file particular
to these figures for the paper, temporal_coding_figures_errors.pklz.

Note that "load_data" must be in src/ (in addition to save_load_data),
since it is called by four_state_receptor_CS for signal loading. 
However, all i/o functions for figure generation are in src/save_load_data.


Created by Nirag Kadakia at 10:00 04-17-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_binary_errors

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects
from analysis import binary_errors_temporal_run


def calculate_temporal_errors(data_flags, nonzero_bounds=[0.7, 1.3], 
							zero_bound=1./10., threshold_pct_nonzero=100, 
							threshold_pct_zero=100):
	"""
	Calculate temporal errors from the aggregated data
	"""
	
	# Get full iterated variable indices from specss
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(iter_vars[iter_var]))		
	it = sp.nditer(sp.zeros(iter_vars_dims), flags = ['multi_index'])	
	
	# To hold saved CS data for each timepoint and for each iter_var_idx
	data = load_aggregated_temporal_objects(data_flag)
	nT = data['dSs'].shape[0]
	Nn = data['dSs'].shape[-1]
	init_objs = sp.reshape(data['init_objs'], iter_vars_dims)
	
	# Data structures to hold error data
	error_dims = tuple(iter_vars_dims) + (nT, )
	errors_nonzero = sp.zeros(error_dims)
	errors_zero = sp.zeros(error_dims)
	
	while not it.finished:
	
		init_CS_object = init_objs[it.multi_index]
		
		dSs_idxs = [slice(None)]*len(data['dSs'].shape)
		dSs_idxs[1:-1] = it.multi_index 
		dSs = data['dSs'][dSs_idxs]
		dSs_est = data['dSs_est'][dSs_idxs]
		
		error_idxs = [slice(None)]*len(errors_nonzero.shape)
		error_idxs[:-1] = it.multi_index
		
		binary_error = binary_errors_temporal_run(init_CS_object, dSs, dSs_est, 
						nonzero_bounds=nonzero_bounds, zero_bound=zero_bound)
		errors_nonzero[error_idxs] = binary_error['errors_nonzero']
		errors_zero[error_idxs] = binary_error['errors_zero']
	
		it.iternext()

	save_binary_errors(errors_nonzero, errors_zero, data_flag)
	
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	calculate_temporal_errors(data_flag)
