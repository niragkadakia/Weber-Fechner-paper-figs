"""
Plot estimation error in temporal coding tasks in a false color plot.

Created by Nirag Kadakia at 22:00 04-16-2018
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
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									save_fig
from plot_formats import fig_temporal_heatmap

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects, load_signal_trace_from_file


def plot_temporal_errors(data_flag, nT_to_plot=85, Kk_idx=2, 
							rate_idxs=range(6)):
	
	success = load_success_ratios(data_flag)
	list_dict = read_specs_file(data_flag)
	signal_file = list_dict['fixed_vars']['signal_trace_file']
	signal_trace = load_signal_trace_from_file(signal_file)
	
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	
	assert len(iter_vars) == 3, "Need 3 iter_vars"
	iter_var_names = ['temporal_adaptation_rate', 'seed_Kk2', 'Kk']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)

	Tt = signal_trace[:nT_to_plot, 0]
	rates = iter_vars['temporal_adaptation_rate'][rate_idxs]
	X, Y = sp.meshgrid(Tt, rates)
	
	avg_success = sp.average(success, axis=2)
	Z = avg_success[:nT_to_plot, rate_idxs, Kk_idx]
	
	fig = fig_temporal_heatmap()
	plt.pcolormesh(X, Y, Z.T, cmap='hot', vmin=-0.05, vmax=1.05)
	plt.yscale('log')
	save_fig('temporal_heatmap', subdir=data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)