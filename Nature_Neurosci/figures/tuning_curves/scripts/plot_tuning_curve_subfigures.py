"""
Plot tuning curves for a particularr specs file

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
import matplotlib.pyplot as plt
from matplotlib import cm


def plot_tuning_curve_subfigures(data_flag, plot_tuning_curves=True, 
									plot_Kk2=True, plot_firing_rate=True):
	
	# Which level of background stimulus and diversity, and how many figures?
	mu_dSs_idx = 3
	sigma_Kk2_idx = 1
	num_figs = 9
	
	# Which background level indices to plot for adaptation plots?
	mu_dSs_idxs = [0, 8, 13]
	
	# Which sparsity for adaptation plots?
	adaptation_Kk = 6
	
	# What are the random number seeds for picking odor components?
	odor_seeds = range(40)
	
	# Which receptors to highlight; choose these by ordering Kk2
	highlight_figs = [2, 3, 7]
	highlight_colors = [cm.Blues, cm.Oranges, cm.Greens]
	
	list_dict = read_specs_file(data_flag)
	for key in list_dict:
		exec("%s = list_dict[key]" % key)
	
	tuning_curve_data = load_tuning_curve(data_flag)
	tuning_curve = tuning_curve_data['tuning_curve']
	Kk2s = tuning_curve_data['Kk2s']
	
	# All receptors to plot
	iMs_to_plot = sp.linspace(6, params['Mm'] - 1, num_figs, dtype='int')


	#############################
	###     Tuning Curves     ###
	#############################
	
	if plot_tuning_curves == True:
		
		# Plot the given receptor idxs
		highlight_idx = 0
		for fig_num, iM in enumerate(iMs_to_plot):
			
			fig = tuning_curve_subfigures()
			
			# Colors for highlighted figures versus non-highlighted figures
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx](0.5)
				highlight_idx += 1
				lw = 3.
			else:
				color = cm.Greys(0.8)
				lw = 2.
				
			# Plot every even value on left; odd on right to show symmetry
			LHS_curve = tuning_curve[mu_dSs_idx, sigma_Kk2_idx, ::2, iM]
			LHS_curve = sp.sort(LHS_curve)
			RHS_curve = tuning_curve[mu_dSs_idx, sigma_Kk2_idx, 1:-1:2, iM]
			RHS_curve = sp.sort(RHS_curve)[::-1]
			
			# Connect right side and left side
			RHS_curve = sp.hstack((LHS_curve[-1], sp.sort(RHS_curve)[::-1]))
			
			plt.plot(sp.arange(params['Nn']/2), LHS_curve, color=color, lw=lw)
			plt.plot(sp.arange(params['Nn']/2 - 1, params['Nn'] - 1), 
								RHS_curve, color=color, lw=lw)
			
			save_tuning_curve_fig(fig, mu_dSs_idx, sigma_Kk2_idx, 
									fig_num, data_flag)
			
		
	#############################
	###          Kk2          ###
	#############################
	
	if plot_Kk2 == True:
		
		fig = Kk2_subfigure()
		ax = fig.add_subplot(111)
		plt.imshow(sp.log(Kk2s[mu_dSs_idx, sigma_Kk2_idx, :, :].T), aspect=0.2, 
					cmap='bone', rasterized=True, vmin=-8.5, vmax=-3)
		
		# Draw arrows to indicate representative plots
		highlight_idx = 0
		for fig_num, iM in enumerate(iMs_to_plot):
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx](0.5)
				highlight_idx += 1
				ax.annotate('', fontsize=20, xy=(iM, 2), xycoords='data', 
					xytext=(0, 30), textcoords='offset points', 
					arrowprops=dict(arrowstyle="->", lw=2.5, color=color))
			else:
				continue
			
		cb = plt.colorbar(ticks=[-8, -6, -4])
		cb.ax.tick_params(labelsize=14)
		save_Kk2_fig(fig, sigma_Kk2_idx, data_flag)
	
	
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

				# Array to hold full time trace; mu_dSs incremented
				num_mu_dSs_to_plot = len(mu_dSs_idxs)
				response = sp.zeros(num_mu_dSs_to_plot*100)
				
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
					len_stim = int(len(response)/7)
					stim_beg = (2*idSs + 1)*len_stim
					stim_end = (2*idSs + 2)*len_stim
					response[stim_beg: stim_end] = obj.Yy[iM]
					
					# Set xlimits based on length of first stimulus
					if idSs == 0:
						plt.xlim(stim_beg/2.0, len(response) - stim_beg/2.0)
						
					# Red region when odor is on; only call once per plot
					if iM == 0:
						plt.axvspan(stim_beg, stim_end, color=
									cm.Reds(0.25*(idSs + 1)/len(mu_dSs_idxs)),
									zorder=-1000)
				
			
				# Smooth edges to look more natural
				response = gaussian_filter(response, sigma=2)
				plt.plot(response, color=color, lw=lw, zorder=zorder)
				
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
	