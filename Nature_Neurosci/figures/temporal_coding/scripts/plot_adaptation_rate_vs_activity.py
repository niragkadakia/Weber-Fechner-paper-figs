"""
Plot adaptation rate versus mean activity level.

Created by Nirag Kadakia at 22:30 01-28-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
sys.path.append('../src')
import matplotlib.pyplot as plt
from load_specs import read_specs_file
from save_load_data import load_temporal_errors, save_plot_adapt_rate_vs_act_fig
from figure_plot_formats import adapt_rate_vs_act


def plot_adapt_rate_vs_act(data_flag,  iT=0, iter_var_idx_to_plot=0, 
							avg_var_idx_to_plot=0):
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	Yys = data['Yy'][iT, iter_var_idx_to_plot, avg_var_idx_to_plot]
	adaptation_rates = data['adaptation_rates'][iT, 
					iter_var_idx_to_plot, avg_var_idx_to_plot]
	Yys_ordered_idxs = sp.argsort(Yys)
	adaptation_rates_ordered = adaptation_rates[Yys_ordered_idxs]
	Yys_ordered = Yys[Yys_ordered_idxs]
	
	colors = plt.cm.Reds(sp.linspace(0.25, 0.75, len(adaptation_rates)))
	fig = adapt_rate_vs_act(xvals=Yys, yvals=adaptation_rates)
	plt.scatter(Yys_ordered, adaptation_rates_ordered, s=20, color=colors)
	save_plot_adapt_rate_vs_act_fig(fig, data_flag)

	
if __name__ == '__main__':
	data_flags = sys.argv[1]
	plot_adapt_rate_vs_act(data_flags, iT=95, iter_var_idx_to_plot=8)