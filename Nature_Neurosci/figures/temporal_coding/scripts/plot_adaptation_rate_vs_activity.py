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
from save_load_data import load_temporal_errors


def plot_adapt_rate_vs_act(data_flag,  iT=100, iter_var_idx_to_plot=7, 
							avg_var_idxs_to_plot=[0, 1, 5]):
	"""
	"""
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	
	for avg_var_idx in avg_var_idxs_to_plot:
		adaptation_rates = data['adaptation_rates'][iT, 
						iter_var_idx_to_plot, avg_var_idx]
		Yys = data['Yy'][iT, iter_var_idx_to_plot, avg_var_idx]
		dYys = data['dYy'][iT, iter_var_idx_to_plot, avg_var_idx]
	
		plt.scatter(Yys, adaptation_rates)
		plt.xlim(min(Yys)*0.8, max(Yys)*1.2)
		plt.ylim(min(adaptation_rates)*0.8, max(adaptation_rates)*1.2)
		plt.yscale('log')	
		
	plt.show()
	
	
if __name__ == '__main__':
	data_flags = sys.argv[1]
	plot_adapt_rate_vs_act(data_flags)