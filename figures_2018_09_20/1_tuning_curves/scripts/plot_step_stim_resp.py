"""
Plot tge firing rate response to a short pulse stimulus.


Created by Nirag Kadakia at 11:00 10-01-2018
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

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir, scripts_dir
sys.path.append(src_dir())
from utils import merge_two_dicts, get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from temporal_CS_run import temporal_CS_run

def plot_step_stim_resp(data_flag, Ss0_to_plot=1e1, seed_to_plot=0):
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	Nn = list_dict['params']['Nn']
	Mm = list_dict['params']['Mm']
	
	corr_iter_vars = ['seed_dSs', 'Kk_1', 'Kk_2']
	for iK, key in enumerate(corr_iter_vars):
		assert list(iter_vars.keys())[iK] == key, "iK iter_var must be %s" % key
	vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
	obj = four_state_receptor_CS(**vars_to_pass)
	
	temporal_CS_run(data_flag, iter_var_idxs,
					mu_dSs_offset=0, mu_dSs_multiplier=1./3., 
					sigma_dSs_offset=0, sigma_dSs_multiplier=1./9., 
					signal_window=None, save_data=True, 
					decode=True)
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	iter_var_idxs = sys.argv[2:]
	plot_step_stim_resp(data_flag)