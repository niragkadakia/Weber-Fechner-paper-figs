"""
Plot estimation of odor signal intensity and identity.

Created by Nirag Kadakia at 12:40 04-11-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig
from plot_formats import fig_signal_estimation

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_estimated_signal(data_flag, mu_Ss0_to_plot=28, Kk_to_plot=2, 
							conc_shift=0):
	"""
	Args:
		data_flag: string; run identifier.
		mu_Ss0_to_plot: int; index of the signal intensity to plot
		Kk_to_plot: int; index of the Kk sparsity to plot
		
	"""
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	
	assert len(iter_vars) == 3, "Need 3 iter_vars"
	
	print ('Loading object list...'),
	CS_object_array = load_aggregated_object_list(iter_vars_dims, data_flag)
	print ('...loaded.')
	
	mu_Ss0_var = iter_vars.keys()[0]
	Kk2_seed_var = iter_vars.keys()[1]
	Kk_var = iter_vars.keys()[2]
	
	color_lo = plt.cm.Purples(0.85)
	color_hi = plt.cm.Greens(0.5)
	
	for iSeed in range(iter_vars_dims[1]):
		obj = CS_object_array[mu_Ss0_to_plot, iSeed, Kk_to_plot]
		obj.dSs = obj.dSs*10.**conc_shift
		obj.dSs_est = obj.dSs_est*10.**conc_shift
		idxs_sorted = sp.argsort(obj.dSs)[::-1]
		nonzero_idxs = idxs_sorted[:obj.Kk]
		zero_idxs = idxs_sorted[obj.Kk:]
		
		fig = fig_signal_estimation()
		gridspec.GridSpec(1, 7)
		
		# Nonzero components; ordered
		ax1 = plt.subplot2grid((1, 7), (0, 0))
		plt.ylim(-max(obj.dSs)*1.2, max(obj.dSs)*1.2)
		plt.xticks([])
		ax1.tick_params(labelsize=15)
		plt.fill_between(sp.arange(obj.Kk), 0, obj.dSs[nonzero_idxs], 
							color=color_lo)
		plt.fill_between(sp.arange(obj.Kk), 0, obj.dSs_est[nonzero_idxs], 
							color=color_hi, alpha=0.75)
		
		# Zero components
		ax2 = plt.subplot2grid((1, 7), (0, 1), colspan=5)
		plt.fill_between(sp.arange(obj.Nn - obj.Kk), 0, 
							sp.sort(obj.dSs_est[zero_idxs])[::-1], 
							color=color_hi, alpha=0.75)
		plt.ylim(-min(obj.dSs[nonzero_idxs])/2., min(obj.dSs[nonzero_idxs])/2.)
		ax2.yaxis.tick_right()
		plt.xticks([])
		ax2.tick_params(labelsize=15)
		save_fig('est_signal_[%s,%s,%s]' % (mu_Ss0_to_plot, iSeed, 
									Kk_to_plot), subdir=data_flag, 
									tight_layout=False)
				
				
if __name__ == '__main__':
	data_flag = get_flag()
	plot_estimated_signal(data_flag)