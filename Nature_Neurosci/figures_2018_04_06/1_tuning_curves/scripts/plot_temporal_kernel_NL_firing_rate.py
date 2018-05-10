"""
Plot the temporal filter for visualization


Created by Nirag Kadakia at 14:37 04-09-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.ndimage.filters import gaussian_filter
sys.path.append('../../shared_src')
from plot_formats import fig_temporal_kernel
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_temporal_kernel_NL_firing_rate(data_flag):
	"""
	Plot a trace of the temporal kernel.
	"""
	
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	
	print ('Loading object list...'),
	CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
	print ('...loaded.')
	obj = CS_object_array[tuple([0]*len(iter_vars_dims))]
	
	kernel_Tt = sp.arange(0, 0.075, 1e-3)
	NL_Tt = sp.arange(-1, 1.6, 1e-3)
	temporal_kernel = obj.kernel_scale*(gamma.pdf(kernel_Tt, 2, 
						scale=obj.kernel_tau_1) - obj.kernel_alpha*\
						gamma.pdf(kernel_Tt, 3, scale=obj.kernel_tau_2))
	NL = NL_Tt*(NL_Tt > 0)
	firing_rate = sp.exp(-NL)*(NL_Tt > 0)
	firing_rate = gaussian_filter(firing_rate, sigma=50)
	
	fig = fig_temporal_kernel()
	ax = plt.axes(frameon=False)
	ax.axis('off')
	plt.plot(kernel_Tt, temporal_kernel, lw=7, color='k')
	save_fig('temporal_kernel', subdir=data_flag)
	
	fig = fig_temporal_kernel()
	ax = plt.axes(frameon=False)
	ax.axis('off')
	plt.plot(NL_Tt, NL, lw=7, color='k')
	save_fig('NL', subdir=data_flag)
		
	fig = fig_temporal_kernel()
	ax = plt.axes(frameon=False)
	plt.plot(firing_rate, lw=7, color='k')
	ax.axis('off')
	save_fig('r(t)', subdir=data_flag)
	
	
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_temporal_kernel_NL_firing_rate(data_flag)