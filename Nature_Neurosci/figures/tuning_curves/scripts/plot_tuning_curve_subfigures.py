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
import sys
sys.path.append('../src')
from utils import get_flag
from load_specs import read_specs_file, parse_iterated_vars, \
						parse_relative_vars
from utils import get_flags, merge_two_dicts
from save_load_data import load_tuning_curve, save_tuning_curve_fig, \
							save_Kk2_fig
from figure_plot_formats import tuning_curve_subfigures, \
								Kk2_subfigure
import matplotlib.pyplot as plt
from matplotlib import cm


def plot_tuning_curve_subfigures(data_flag, plot_tuning_curves=True, 
									plot_Kk2=True, plot_adaptation=True):
	
	# Which level of background stimulus and diversity, and how many figures?
	mu_dSs_idx = 3
	sigma_Kk2_idx = 1
	num_figs = 9
	
	# Which receptors to highlight; choose these by ordering Kk2
	highlight_figs = [2, 3, 7]
	highlight_colors = [cm.Blues(0.5), cm.Oranges(0.5), cm.Greens(0.5)]
	highlight_idx = 0
	
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
		for fig_num, iM in enumerate(iMs_to_plot):
			
			fig = tuning_curve_subfigures()
			
			# Colors for highlighted figures versus non-highlighted figures
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx]
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
		ax = fig.add_subplot(111)#, autoscale_on=False)
		plt.imshow(sp.log(Kk2s[mu_dSs_idx, sigma_Kk2_idx, :, :].T), aspect=0.2, 
					cmap='bone', rasterized=True, vmin=-8.5, vmax=-3)
		
		highlight_idx = 0
		for fig_num, iM in enumerate(iMs_to_plot):
			if fig_num in highlight_figs:
				color = highlight_colors[highlight_idx]
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
	###       Adaptive        ###
	#############################
	
	if plot_adaptation == True:
		pass
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_tuning_curve_subfigures(data_flag)
	