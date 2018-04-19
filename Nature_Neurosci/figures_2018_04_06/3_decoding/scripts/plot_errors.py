"""
Plot estimation of odor signal intensity and identity, and the 
full errors which is the product of these.

Created by Nirag Kadakia at 12:40 04-11-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig
from plot_formats import fig_errors_vs_Kk

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_errors(data_flag, zero_thresh=0.1, nonzero_thresh=[0.7, 1.3],
							row_placement=[[1, 1], [1, 1], [1, 1]]):
	"""
	Args:
		data_flag: string; run identifier.
		conc_shift: integer; shift of concentration for plotting.
		zero_thresh: multiplier of mu_dSs for which to consider a false
			negative or false positive for the estimated odor signal.
		nonzero_thresh: 2-element list; range for which |s_est/s| is 
			considered correctly intensity decoded
		row_placement: 3-element list of 2-element lists. Entries are 
			whether to draw tick labels for x and y-axes for intensity, 
			identity, and full errors plots.
		
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
	Kk_var = iter_vars.keys()[2]
	x = sp.log(iter_vars[mu_Ss0_var])/sp.log(10) + conc_shift
	y = iter_vars[Kk_var]
	X, Y = sp.meshgrid(x, y)
		
	zero_errors = sp.zeros((iter_vars_dims[0], iter_vars_dims[2]))
	nonzero_errors = sp.zeros((iter_vars_dims[0], iter_vars_dims[2]))
	full_errors = sp.zeros((iter_vars_dims[0], iter_vars_dims[2]))
	
	for iSs0 in range(iter_vars_dims[0]):
		for iKk in range(iter_vars_dims[2]):
			
			zero_errors_all_seeds = []
			nonzero_errors_all_seeds = []
			full_all_seeds = []
			
			for iSeed in range(iter_vars_dims[1]):
				obj = CS_object_array[iSs0, iSeed, iKk]
				
				# Get zero and nonzero indices of signal
				all_idxs = sp.arange(obj.Nn)
				nonzero_idxs = obj.idxs[0]
				mask = sp.ones(obj.Nn, dtype=bool)
				mask[nonzero_idxs] = False
				zero_idxs = all_idxs[mask]
				
				# Count number of nonzero components within error
				scaled_est = sp.absolute(obj.dSs_est[nonzero_idxs]/
											obj.dSs[nonzero_idxs])
				num_good_nonzeros = 100.*sp.sum((scaled_est 
									<= nonzero_thresh[1])*(scaled_est 
									> nonzero_thresh[0]))/obj.Kk
				
				# Count number of zero components that are below min threshold
				# and number of nonzero components above threshold
				num_good_zeros = sp.sum(abs(obj.dSs_est[zero_idxs]) <= 
								abs(zero_thresh*obj.mu_dSs)) + \
								sp.sum(abs(obj.dSs_est[nonzero_idxs]) >= 
								abs(zero_thresh*obj.mu_dSs))
				
				# Binary success for identity and intensity, separately
				nonzero_errors_all_seeds.append(num_good_nonzeros == 100)
				zero_errors_all_seeds.append(num_good_zeros == 100)
				full_all_seeds.append((num_good_nonzeros == 100)*
										(num_good_zeros == 100))
				
			zero_errors[iSs0, iKk] = 100.*sp.average(zero_errors_all_seeds)
			nonzero_errors[iSs0, iKk] = 100.*sp.average(nonzero_errors_all_seeds)
			full_errors[iSs0, iKk] = 100.*sp.average(full_all_seeds)
		
	errors = dict()
	errors['nonzero'] = nonzero_errors
	errors['zero'] = zero_errors
	errors['full'] = full_errors
	
	vminmax = [0, 100]
	ticks = [0, 50, 100]
	tick_labels = ['0', '50', '100']
	row_placement = row_placement
	
	for key in errors.keys():
	
		fig = fig_errors_vs_Kk()
		
		# Whether or not to put axes
		if key == 'nonzero':
			if row_placement[0][0] == 0:
				plt.xticks([])
			if row_placement[0][1] == 0:
				plt.yticks([])
		if key == 'zero':
			if row_placement[1][0] == 0:
				plt.xticks([])
			if row_placement[1][1] == 0:
				plt.yticks([])
		if key == 'full':
			if row_placement[2][0] == 0:
				plt.xticks([])
			if row_placement[2][1] == 0:
				plt.yticks([])
		
		plt.xlim(0, 4)
		plt.ylim(1, 7)
		plt.pcolormesh(X, Y, errors[key].T, cmap=plt.cm.hot, 
						rasterized=True, shading='gouraud', 
						vmin=vminmax[0], vmax=vminmax[1])
		save_fig('%s_errors' % key, subdir=data_flag)
		
		# Separate figure for colorbar
		fig = plt.figure()
		fig.set_size_inches(1, 5)
		ax1 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
		norm = mpl.colors.Normalize(vmin=vminmax[0], vmax=vminmax[1])
		cbar = mpl.colorbar.ColorbarBase(ax1, cmap=plt.cm.hot,
                                norm=norm, ticks=ticks,
                                orientation='vertical')
		cbar.ax.tick_params(labelsize=25)
		cbar.ax.set_yticklabels(tick_labels)
		save_fig('%s_errors_colorbar' % key, subdir=data_flag, 
					tight_layout=False)
		
		
				
if __name__ == '__main__':
	data_flag = get_flag()
	plot_errors(data_flag)