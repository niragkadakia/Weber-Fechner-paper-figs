"""
Plot tuning curves for a particularr specs file

Created by Nirag Kadakia at 23:00 10-04-2017
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
from encode_CS import single_encode_CS
from save_load_data import load_tuning_curve, save_tuning_curve_fig
from figure_plot_formats import tuning_curve_plot_Kk2
import matplotlib.pyplot as plt
from matplotlib import cm


def plot_tuning_curves(data_flags):
	
	# First entries are for mu_dSs, second are for Kk2 diversity
	plot_vars = [[0, 9, 19], [8, 15, 19]]
	cmaps = [[cm.Greys, cm.Purples, cm.Blues], 
				[cm.Greys, cm.Purples, cm.Blues]]
				
	for data_idx, data_flag in enumerate(data_flags):

		list_dict = read_specs_file(data_flag)
		for key in list_dict:
			exec("%s = list_dict[key]" % key)
		
		tuning_curve_data = load_tuning_curve(data_flag)
		tuning_curve = tuning_curve_data['tuning_curve']
		epsilons = tuning_curve_data['epsilons']
		Kk2s = tuning_curve_data['Kk2s']
		
		if data_idx == 0: 
			fig, plot_dims, axes_tuning, axes_Kk2, axes_signal = \
					tuning_curve_plot_Kk2(plot_vars, iter_vars, params)
		
		for idx, idx_var in enumerate(plot_vars[0]):
			for idy, idy_var in enumerate(plot_vars[1]):
			
				colors = cmaps[data_idx][idx](sp.linspace(0.75, 0.3, 
												params['Mm']))
				
				for iM in range(params['Mm']):
					axes_tuning[idx, idy].plot(
						sp.arange(params['Nn']/2), 
						sp.sort(tuning_curve[idx_var, idy_var, ::2, iM]), 
						color=colors[iM], linewidth = 0.7,
						zorder = params['Mm'] - iM)
					axes_tuning[idx, idy].plot(
						sp.arange(params['Nn']/2 - 1, params['Nn']-1), 
						sp.sort(tuning_curve[idx_var, idy_var, 1::2, iM])[::-1],
						color=colors[iM], linewidth = 0.7,
						zorder = params['Mm'] - iM)
				
				
				if idx == 0:
					sorted_idxs = sp.argsort(sp.std(Kk2s[0, idy_var, : , :], axis = 1))
					axes_Kk2[idy].imshow(Kk2s[0, idy_var, sorted_idxs, :].T, 
											aspect=0.3, cmap='bone', 
											rasterized=True)
				
	save_tuning_curve_fig(fig, data_flag)
	
if __name__ == '__main__':
	data_flags = get_flags()
	plot_tuning_curves(data_flags)
	