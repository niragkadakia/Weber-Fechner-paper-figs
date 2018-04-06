"""
Plot tuning curve matrix for all signal concentrations

spec files used to create data (data_flags):
figure_tuning_curves_adaptation_Kk2
figure_tuning_curves_no_adaptation_Kk2

Current figure (Figure_tuning_curves.svg) uses odor_seed 12 and 22

Created by Nirag Kadakia at 17:30 02-14-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
from scipy.ndimage.filters import gaussian_filter
import sys
sys.path.append('../src')
import matplotlib.pyplot as plt
params = {'text.usetex': False, 'mathtext.fontset': 'dejavusans'}
plt.rcParams.update(params)
from matplotlib import cm
from utils import get_flag
from load_specs import read_specs_file, parse_iterated_vars, \
						parse_relative_vars
from utils import get_flags, merge_two_dicts
from save_load_data import load_tuning_curve, save_tuning_curve_matrix_fig, \
							save_tuning_curve_sparse_odors_matrix_fig, \
							save_tuning_curve_std_fig
from figure_plot_formats import tuning_curve_matrix_subfigure, \
								tuning_curve_std_subfigures
from encode_CS import single_encode_CS


def plot_tuning_curve_matrix(data_flag):
	
	# Which diversity of the system to plot?
	sigma_Kk2_idx = 2
	
	# Wha signal indices to plot?
	mu_Ss0_idxs = range(0, 15)
	
	# Receptors to plot
	receptors_to_plot = [4, 14, 24]
	
	# Odor seeds and components
	odor_seeds = sp.arange(100)
	num_odor_comps = 6
	
	# Which standard deviation figures to highlight
	highlight_colors = [cm.Greens, cm.Blues, cm.Oranges]
	
	list_dict = read_specs_file(data_flag)
	for key in list_dict:
		exec("%s = list_dict[key]" % key)
	
	
	
	#############################
	### Tuning to sparse odor ###
	#############################
	
	std_devs = sp.zeros((len(mu_Ss0_idxs), len(receptors_to_plot)))
	
	for iM, receptor_idx in enumerate(receptors_to_plot):
	
		fig = tuning_curve_matrix_subfigure()
		data_to_plot = sp.zeros((len(mu_Ss0_idxs), len(odor_seeds)))
		
		for iSeed, odor_seed in enumerate(odor_seeds):
			
			# Loop over different background stimuli; plot at successive dt
			for idSs, mu_dSs_idx in enumerate(mu_Ss0_idxs):
					
				iter_var_idxs = [mu_dSs_idx, sigma_Kk2_idx]
				vars_to_pass = dict()
				vars_to_pass = parse_iterated_vars(iter_vars, 
													iter_var_idxs, 
													vars_to_pass)
				vars_to_pass = parse_relative_vars(rel_vars, iter_vars, 
													vars_to_pass)
				vars_to_pass = merge_two_dicts(vars_to_pass, fixed_vars)
				vars_to_pass = merge_two_dicts(vars_to_pass, params)
				
				# Choose the random stimulus of num_odor_comps components
				sp.random.seed(odor_seed)
				dSs_idxs = sp.random.choice(sp.arange(params['Nn']), 
											size=num_odor_comps, 
											replace=False)
				
				# Call object and get activity; plot in appropriate region
				vars_to_pass['manual_dSs_idxs'] = dSs_idxs
				obj = single_encode_CS(vars_to_pass, run_specs)
				data_to_plot[idSs, iSeed] = obj.Yy[receptor_idx]
		plt.pcolormesh(data_to_plot.T, vmin=0, vmax=1, cmap=plt.cm.hot, 
						rasterized=True)
		plt.colorbar()
		save_tuning_curve_sparse_odors_matrix_fig(fig, sigma_Kk2_idx, 
													receptor_idx, data_flag)
		
		std_devs[:, iM] = sp.std(data_to_plot, axis=1)

		
		
	#############################
	###  Dev of tuning curve  ###
	#############################
		
	fig = tuning_curve_std_subfigures()
	for iM in range(len(receptors_to_plot)):
		plt.plot(range(len(mu_Ss0_idxs)), std_devs[:, iM], 
				color=highlight_colors[iM](0.8), alpha=0.65, lw=2.5)
	
	save_tuning_curve_std_fig(fig, sigma_Kk2_idx, receptors_to_plot, data_flag)
	
	

	#############################
	###  Full tuning curve    ###
	#############################
	
	tuning_curve_data = load_tuning_curve(data_flag)
	
	for iM, receptor_idx in enumerate(receptors_to_plot):
		fig = tuning_curve_matrix_subfigure()
		data = tuning_curve_data['tuning_curve'][mu_Ss0_idxs, 
				sigma_Kk2_idx, :, receptor_idx]
		plt.pcolormesh(data.T, cmap=plt.cm.hot, rasterized=True)
		plt.colorbar()
		save_tuning_curve_matrix_fig(fig, sigma_Kk2_idx, 
										receptor_idx, data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_tuning_curve_matrix(data_flag)
	