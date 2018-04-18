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
from save_load_figure_data import load_binary_errors, load_success_ratios, save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects


def plot_temporal_errors(data_flag, complexities_to_plot=range(1, 6)):
	"""
	"""
	
	data_dict = load_binary_errors(data_flag)
	success = load_success_ratios(data_flag)
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	
	assert len(iter_vars) == 3, "Need 3 iter_vars"
	iter_var_names = ['temporal_adaptation_rate', 'seed_Kk2', 'Kk']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	adapt_rate_vals = iter_vars['temporal_adaptation_rate']
	
	zero = data_dict['errors_zero']
	avg_zero = sp.average(zero, axis=1)

	rates_to_plot = [0, 1, 2]
	Kk_to_plot = 1
	
	plt.subplot(311)
	Tt = range(len(data_dict['errors_zero'][:, 0, 0, 0]))
	for rate in rates_to_plot:
		color = 1.*rate/14+0.3
		plt.plot(avg_zero[:, rate, Kk_to_plot].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	
	
	nonzero = data_dict['errors_nonzero']
	avg_nonzero = sp.average(nonzero, axis=2)
	
	plt.subplot(312)
	for rate in rates_to_plot:
		color = 1.*rate/14+0.3
		plt.plot(avg_nonzero[:, rate, Kk_to_plot].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	
	avg_success = sp.average(success, axis=2)
	
	plt.subplot(313)
	for rate in rates_to_plot:
		color = 1.*rate/14+0.3
		plt.plot(avg_success[:, rate, Kk_to_plot].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	
	plt.show()
	
	
	plt.imshow(avg_success[:, rates_to_plot, Kk_to_plot].T, aspect=40, cmap='hot', vmin=-0.05, vmax=1.05)
	plt.colorbar()
	plt.show()
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)