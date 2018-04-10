"""
Plot the tuning curves 


Created by Nirag Kadakia at 11:00 04-06-2017
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
from save_load_figure_data import save_success_ratios, save_fig
from plot_formats import fig_tuning_curve, fig_Kk2

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import merge_two_dicts, get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from encode_CS import single_encode_CS


def plot_tuning_curves(data_flag, Ss0_to_plot=10.0, seed_to_plot=9):
	"""
	Plot tuning curves at a given stimulus intensity level.
	"""
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	Nn = list_dict['params']['Nn']
	Mm = list_dict['params']['Mm']
	tuning_curve = sp.zeros((Nn, Mm))
	
	# Only need to specify the mu_Ss0 and seed_Kk2 index. Kk index not needed.
	assert iter_vars.keys()[0] == 'mu_Ss0', "Need mu_Ss0 as first iter var"
	iter_var_idxs = sp.zeros(len(iter_vars)).astype(int)
	iter_var_idxs[0] = (abs(iter_vars['mu_Ss0'] - Ss0_to_plot)).argmin()
	iter_var_idxs[1] = seed_to_plot	
	
	# Calculate tuning curve
	for iN in range(Nn):
		
		vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
		
		# Overwrite signal to single component at this iN
		vars_to_pass['manual_dSs_idxs'] = sp.array([iN])
		
		# Remove threshold (only needed to limit decoding)
		vars_to_pass['NL_threshold'] = 0
		obj = four_state_receptor_CS(**vars_to_pass)
		obj = single_encode_CS(obj, list_dict['run_specs'])
		tuning_curve[iN, :] = obj.Yy
	
	# Plot the activity matrix
	fig = fig_Kk2()
	ax = plt.gca()
	im = ax.imshow(tuning_curve.T, cmap=plt.cm.bone_r, aspect=1.0, 
						vmin=-1, vmax=300)
	cbar = plt.colorbar(im, fraction=0.07, pad=0.04, orientation="horizontal")
	cbar.set_ticks([0, 100, 200, 300])
	cbar.ax.tick_params(labelsize=12) 
	save_fig('act_matrix_seed=%s' % seed_to_plot, subdir=data_flag)
	
	# Plot sorted individual tuning curve
	for iM in range(Mm):
		fig = fig_tuning_curve()
		tuning_curve_sorted = sp.sort(tuning_curve[:, iM])
		tuning_to_plot = sp.hstack((tuning_curve_sorted[::2], 
									tuning_curve_sorted[1::2][::-1]))
		plt.bar(sp.arange(-Nn/2, Nn/2), tuning_to_plot, 
					color='0.1', width=1)
		save_fig('seed=%s_iM=%s' % (seed_to_plot, iM), subdir=data_flag)
	
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_tuning_curves(data_flag)