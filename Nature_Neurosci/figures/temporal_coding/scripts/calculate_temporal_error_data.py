"""
Calculate estimation error of inferred signal in compressed sensing 
decoding for a full signal trace in time. Success ratios are calculated
as a function of time. Saved to objects as a distinct file particular
to these figures for the paper, temporal_coding_figures_errors.pklz.

Note that "load_data" must be in src/ (in addition to save_load_data),
since it is called by four_state_receptor_CS for signal loading. 
However, all i/o functions for figure generation are in src/save_load_data.


Created by Nirag Kadakia at 10:00 01-23-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
sys.path.append('../src')
from load_specs import read_specs_file
from save_load_data import load_objects, save_temporal_errors
from analysis import binary_errors, binary_success


def calculate_temporal_success(data_flags, nonzero_bounds=[0.75, 1.25], 
							zero_bound=1./30., threshold_pct_nonzero=100, 
							threshold_pct_zero=100):
	"""
	Compile all the data from all temporal traces for all iterated variables.
	Saves zero/nonzero errors, successes, and epsilons for each point
	along the signal trace. 
	"""
	
	# Get full iterated variable indices from specss
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(iter_vars[iter_var]))		
	it = sp.nditer(sp.zeros(iter_vars_dims), flags = ['multi_index'])	
	
	# To hold data for each timepoint, for each iter_var_idx
	data = dict()
	
	# Grab data from initial point; not needed in every iteration
	CS_init_array = load_objects(list(it.multi_index), data_flag)
	nT = len(CS_init_array)
	data['Tt'] = CS_init_array[0].signal_trace_Tt
	data['signal'] = CS_init_array[0].signal_trace
	if 'Kk_split' in list_dict['params'].keys():
		if list_dict['params']['Kk_split'] != 0:
			signal_2 = CS_init_array[0].signal_trace_2
			data['signal_2'] = signal_2
	
	# Data structures to hold stats at every iteration
	array_shape = sp.hstack((nT, iter_vars_dims))
	data['success_ratios'] = sp.zeros(array_shape)
	data['nonzero_errors'] = sp.zeros(array_shape)
	data['zero_errors'] = sp.zeros(array_shape)
	data['avg_eps'] = sp.zeros(array_shape)
	data['avg_dYy'] = sp.zeros(array_shape)
	data['avg_Yy'] = sp.zeros(array_shape)
	if 'Kk_split' in list_dict['params'].keys():
		if list_dict['params']['Kk_split'] != 0:
			data['success_ratios_2'] = sp.zeros(array_shape)
			data['nonzero_errors_2'] = sp.zeros(array_shape)
	
	# Non-scalar, array-shape variables 
	array_shape_dSs = sp.hstack((nT, iter_vars_dims, list_dict['params']['Nn']))
	data['dSs_est'] = sp.zeros(array_shape_dSs)
	data['dSs'] = sp.zeros(array_shape_dSs)
	
	array_shape_Mm = sp.hstack((nT, iter_vars_dims, list_dict['params']['Mm']))
	data['dYy'] = sp.zeros(array_shape_Mm)
	data['Yy'] = sp.zeros(array_shape_Mm)
	data['epsilons'] = sp.zeros(array_shape_Mm)
	data['adaptation_rates'] = sp.zeros(array_shape_Mm)
	
	while not it.finished:
		
		# Load temporal data trace for this particular iterated variable index
		print 'Loading index:', it.multi_index
		temporal_CS_array = load_objects(list(it.multi_index), data_flag)
		
		# Grab all the errors and stats, timepoint-by-timepoint
		for iT, dt in enumerate(data['Tt']):
		
			full_idx = (iT, ) + it.multi_index
			
			if temporal_CS_array[iT].Kk_split == 0:
				errors = binary_errors(temporal_CS_array[iT], 
									nonzero_bounds=nonzero_bounds,
									zero_bound=zero_bound)
				success = binary_success(
						errors['errors_nonzero'], errors['errors_zero'], 
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
			elif temporal_CS_array[iT].Kk_split > 0:
				errors = binary_errors_dual_odor(temporal_CS_array[iT], 
									nonzero_bounds=nonzero_bounds,
									zero_bound=zero_bound)
				success = binary_success(
						errors['errors_nonzero'], errors['errors_zero'], 
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
				success_2 = binary_success(
						errors['errors_nonzero_2'], errors['errors_zero'],
						threshold_pct_nonzero=threshold_pct_nonzero,
						threshold_pct_zero=threshold_pct_zero)
				
			data['nonzero_errors'][full_idx] = errors['errors_nonzero']
			data['zero_errors'][full_idx] = errors['errors_zero']
			data['success_ratios'][full_idx] = success
			data['avg_eps'][full_idx] = sp.average(temporal_CS_array[iT].eps)
			data['avg_dYy'][full_idx] = sp.average(temporal_CS_array[iT].dYy)
			data['avg_Yy'][full_idx] = sp.average(temporal_CS_array[iT].Yy)
			data['dSs_est'][full_idx] = temporal_CS_array[iT].dSs_est
			data['dSs'][full_idx] = temporal_CS_array[iT].dSs
			data['dYy'][full_idx] = temporal_CS_array[iT].dYy
			data['Yy'][full_idx] = temporal_CS_array[iT].Yy
			data['epsilons'][full_idx] = temporal_CS_array[iT].eps
			if temporal_CS_array[iT].Kk_split > 0:
				data['nonzero_errors_2'][full_idx] = errors['errors_nonzero_2']
				data['success_ratios_2'][full_idx] = success_2

			if 'temporal_adaptation_rate_sigma' in list_dict['fixed_vars'].keys():
				if list_dict['fixed_vars']['temporal_adaptation_rate_sigma'] != 0:
					data['adaptation_rates'][full_idx] = \
						temporal_CS_array[iT].temporal_adaptation_rate_vector
			
	
		it.iternext()
				
	save_temporal_errors(data, data_flag)
		
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	calculate_temporal_success(data_flag)