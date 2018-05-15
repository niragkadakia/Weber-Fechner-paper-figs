"""
Calculate tuning curves for a particularr specs file

Two figures generated with separate command line arguments
figure_tuning_curves_adaptation_Kk2
figure_tuning_curves_no_adaptation_Kk2


Created by Nirag Kadakia at 23:00 10-04-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import os
sys.path.append('../src')
from utils import get_flag
from load_specs import read_specs_file, parse_iterated_vars, \
						parse_relative_vars
from utils import merge_two_dicts, get_flag
from encode_CS import single_encode_CS
from save_load_data import save_tuning_curve


def calculate_tuning_curves(data_flag):

	list_dict = read_specs_file(data_flag)
	for key in list_dict:
		exec("%s = list_dict[key]" % key)
	
	# Get the iterated variable dimensions
	iter_vars_dims = []
	for iter_var in iter_vars:
		iter_vars_dims.append(len(iter_vars[iter_var]))		
	it = sp.nditer(sp.zeros(iter_vars_dims), flags = ['multi_index'])	
	
	# Set up array to hold tuning curve curves
	tuning_curve = sp.zeros((iter_vars_dims[0], iter_vars_dims[1], 
									params['Nn'], params['Mm']))
	
	# Set array to hold epsilons, Kk2, and activities
	epsilons = sp.zeros((iter_vars_dims[0], iter_vars_dims[1], params['Mm']))
	Kk2s = sp.zeros((iter_vars_dims[0], iter_vars_dims[1], 
						params['Mm'], params['Nn']))
	activities = sp.zeros((iter_vars_dims[0], iter_vars_dims[1], params['Mm']))
	gains = sp.zeros((iter_vars_dims[0], iter_vars_dims[1], 
						params['Mm'], params['Nn']))
	
	# Iterate tuning curve calculation over all iterable variables 
	while not it.finished:
		iter_var_idxs = it.multi_index
		
		vars_to_pass = dict()
		vars_to_pass = parse_iterated_vars(iter_vars, iter_var_idxs, 
											vars_to_pass)
		vars_to_pass = parse_relative_vars(rel_vars, iter_vars, vars_to_pass)
		vars_to_pass = merge_two_dicts(vars_to_pass, fixed_vars)
		vars_to_pass = merge_two_dicts(vars_to_pass, params)
		
		# Calculate tuning curve
		for iN in range(vars_to_pass['Nn']):
			vars_to_pass['manual_dSs_idxs'] = sp.array([iN])
			obj = single_encode_CS(vars_to_pass, run_specs)
			tuning_curve[iter_var_idxs[0], iter_var_idxs[1], iN, :] = obj.Yy
		
		epsilons[it.multi_index] = obj.eps
		Kk2s[it.multi_index] = obj.Kk2
		activities[it.multi_index] = obj.Yy
		gains[it.multi_index] = obj.Rr
		
		it.iternext()
	
	save_tuning_curve(tuning_curve, epsilons, Kk2s, activities, gains, data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	calculate_tuning_curves(data_flag)
	