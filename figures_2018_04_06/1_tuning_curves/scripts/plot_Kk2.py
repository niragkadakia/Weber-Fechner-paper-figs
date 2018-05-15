"""
Plot the binding matrices and their distribution in a histogram.


Created by Nirag Kadakia at 9:00 04-08-2018
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
from plot_formats import fig_Kk2, fig_Kk2_hist

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_Kk2(data_flag, seed_to_plot=9):
	"""
	Plot binding matrix and histogram.
	"""
	
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	assert list_dict['iter_vars'].keys()[1] == 'seed_Kk2', 'Second iter_var '\
		'must be seed_Kk2'
	assert len(iter_vars_dims) == 3, 'Need 3 iter_vars'	
	
	print ('Loading object list...'),
	CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
	print ('...loaded.')

	Kk2 = CS_object_array[0, seed_to_plot, 0].Kk2
	
	fig = fig_Kk2()
	plt.imshow(sp.log(Kk2)/sp.log(10), cmap=plt.cm.afmhot_r, aspect=1.5)
	save_fig('Kk2_seed=%s' % seed_to_plot, subdir=data_flag)
	
	fig = fig_Kk2_hist()
	hist, bins = sp.histogram(sp.log(sp.ndarray.flatten(Kk2))/sp.log(10), bins = 50, normed=True)
	plt.plot(bins[:-1], sp.log((bins[1:] - bins[:-1])*sp.cumsum(hist))/sp.log(10), color='k', lw=5)
	
	# Convert to log plot
	log_bins = (bins[:-1] + bins[1:])/2.0
	log_vals = sp.log((bins[1:] - bins[:-1])*sp.cumsum(hist))/sp.log(10)
	exp_bins = 10**log_bins
	exp_vals = 10**(log_vals)
	plt.plot(exp_bins, exp_vals, color='k', lw=5)
	plt.xscale('log')
	plt.yscale('log')
	save_fig('Kk2_hist_seed=%s' % seed_to_plot, subdir=data_flag)
	
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_Kk2(data_flag)