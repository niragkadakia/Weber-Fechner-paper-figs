"""
Plot tuning curves for a particularr specs file

spec files used to create data (data_flags):
figure_tuning_curves_adaptation_Kk2
figure_tuning_curves_no_adaptation_Kk2

Current figure (Figure_tuning_curves.svg) uses odor_seed 12 and 22

Created by Nirag Kadakia at 12:00 12-07-2017
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import get_flag
from load_specs import read_specs_file, parse_iterated_vars, \
						parse_relative_vars
from utils import get_flags, merge_two_dicts
from encode_CS import single_encode_CS
from save_load_data import load_tuning_curve, save_tuning_curve_fig, \
							save_Kk2_fig, save_firing_rate_fig, \
							save_firing_rate_stimulus_fig
from figure_plot_formats import tuning_curve_subfigures, \
								Kk2_subfigure, firing_rate_subfigure, \
								firing_rate_stimulus_subfigure


def plot_tuning_curve_subfigures(data_flag, plot_tuning_curves=True, 
									plot_Kk2=True, plot_firing_rate=True):
	
	# Which level of background stimulus and diversity?
	mu_dSs_idx = 2
	sigma_Kk2_idx = 0
	
	# How many odorants to plot in the Kk2 matrix?
	num_odorants_Kk2 = 100
	
	# Which background level indices to plot for adaptation plots?
	mu_dSs_idxs = [0, 3, 6]
	
	# Choose the sparsity of the signal for the adaptation plots.
	adaptation_Kk = 20
	
	# What are the random number seeds for picking odor components?
	odor_seeds = range(40)
	
	# Which of the (num_figs) to highlight, and with what colors? (incr order!)
	highlight_figs = [14, 24, 32]
	highlight_colors = [cm.Blues, cm.Oranges, cm.Greens]
	
	list_dict = read_specs_file(data_flag)
	for key in list_dict:
		exec("%s = list_dict[key]" % key)
	
	tuning_curve_data = load_tuning_curve(data_flag)
	tuning_curve = tuning_curve_data['tuning_curve']
	Kk2s = tuning_curve_data['Kk2s']

		
	#############################
	###          Kk2          ###
	#############################
	
	if plot_Kk2 == True:
		
		fig = Kk2_subfigure()
		ax = fig.add_subplot(111)
		img1 = ax.imshow(sp.log(Kk2s[mu_dSs_idx, sigma_Kk2_idx, :, \
				:num_odorants_Kk2].T)/sp.log(10), aspect=2.0, cmap='hot_r', 
				rasterized=True, vmin=-3.5, vmax=0, interpolation="nearest")
		
		# Draw arrows to indicate representative plots
		highlight_idx = 0
		for fig_num, iM in enumerate(range(params['Mm'])):
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx](0.5)
				highlight_idx += 1
				ax.annotate('', fontsize=20, xy=(iM, 1), xycoords='data', 
					xytext=(0, 30), textcoords='offset points', 
					arrowprops=dict(arrowstyle="->", lw=2.5, color=color))
			else:
				continue
			
		
		# matplotlib weirdness in making colorbar correct size.
		cb = fig.colorbar(img1, ax=ax, ticks=sp.arange(-5, 1, 1))
		ax.set_aspect('auto')
		cb.ax.tick_params(labelsize=14)
		save_Kk2_fig(fig, sigma_Kk2_idx, data_flag)
		
	
	#############################
	###     Tuning Curves     ###
	#############################
	
	if plot_tuning_curves == True:
		
		# Plot the given receptor idxs
		highlight_idx = 0
		for fig_num, iM in enumerate(range(params['Mm'])):
			
			fig = tuning_curve_subfigures()
			
			# Colors for highlighted figures versus non-highlighted figures
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx](0.8)
				highlight_idx += 1
				lw = 3.
			else:
				color = cm.Greys(0.9)
				lw = 2.
				
			# Plot every even value on left; odd on right to show symmetry
			LHS_curve = tuning_curve[mu_dSs_idx, sigma_Kk2_idx, ::2, iM]
			LHS_curve = sp.sort(LHS_curve)
			RHS_curve = tuning_curve[mu_dSs_idx, sigma_Kk2_idx, 1:-1:2, iM]
			RHS_curve = sp.sort(RHS_curve)[::-1]
			
			# Connect right side and left side
			RHS_curve = sp.hstack((LHS_curve[-1], sp.sort(RHS_curve)[::-1]))
			
			plt.bar(sp.arange(params['Nn']/2), LHS_curve, color=color, 
								width=1.2)
			plt.bar(sp.arange(params['Nn']/2 - 1, params['Nn'] - 1), 
								RHS_curve, color=color, width=1.2)
			
			save_tuning_curve_fig(fig, mu_dSs_idx, sigma_Kk2_idx, 
									fig_num, data_flag)
			

	#############################
	###     Firing rates      ###
	#############################
	
	if plot_firing_rate == True:
		
		# Plot firing rate for distinct stimuli
		for odor_seed in odor_seeds:
			fig = firing_rate_subfigure()
			highlight_idx = 0	
			for fig_num, iM in enumerate(range(params['Mm'])):
			
				# Indicate receptors colored blue, orange, green; rest grey
				if fig_num in highlight_figs:
					color = highlight_colors[highlight_idx](0.8)
					lw = 2.0
					highlight_idx += 1
					zorder=iM
				else:
					color = ('0.8')
					lw = 1.75
					zorder=-iM
				
				# Loop over different background stimuli; plot at successive dt
				for idSs, mu_dSs_idx in enumerate(mu_dSs_idxs):
						
					iter_var_idxs = [mu_dSs_idx, sigma_Kk2_idx]
					vars_to_pass = dict()
					vars_to_pass = parse_iterated_vars(iter_vars, 
														iter_var_idxs, 
														vars_to_pass)
					vars_to_pass = parse_relative_vars(rel_vars, iter_vars, 
														vars_to_pass)
					vars_to_pass = merge_two_dicts(vars_to_pass, fixed_vars)
					vars_to_pass = merge_two_dicts(vars_to_pass, params)
					
					# Choose the random stimulus of adaptation_Kk components
					sp.random.seed(odor_seed)
					dSs_idxs = sp.random.choice(sp.arange(params['Nn']), 
												size=adaptation_Kk, 
												replace=False)
					
					# Call object and get activity; plot in appropriate region
					vars_to_pass['manual_dSs_idxs'] = dSs_idxs
					obj = single_encode_CS(vars_to_pass, run_specs)
					stim_beg = (2*idSs + 1)*2
					stim_end = (2*idSs + 2)*2
					x_range = [stim_beg, stim_end]
					response = [obj.Yy[iM]]*2
					plt.plot(x_range, response, color=color, lw=lw, zorder=zorder)
				
			plt.xlim(1, 13)
			save_firing_rate_fig(fig, sigma_Kk2_idx, mu_dSs_idxs, 
									odor_seed, data_flag)
			
			# Plot stimulus vector as a heatmap; roll to double size of bar
			fig = firing_rate_stimulus_subfigure()
			shifted_dSs = obj.dSs + sp.roll(obj.dSs, 1)
			plt.imshow([shifted_dSs], cmap=cm.Greys, aspect=20, 
						vmin=-1, interpolation='nearest')
			save_firing_rate_stimulus_fig(fig, sigma_Kk2_idx, mu_dSs_idxs,
											odor_seed, data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_tuning_curve_subfigures(data_flag)
	