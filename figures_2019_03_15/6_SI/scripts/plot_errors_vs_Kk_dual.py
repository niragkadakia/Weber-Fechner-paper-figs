"""
Plot estimation of odor signal intensity and identity, and the 
full avg_successes which is the product of these. This uses data generated by 
calculate_avg_successes.py.


Created by Nirag Kadakia at 12:40 04-11-2018
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
from save_load_figure_data import save_fig, save_fig_no_whtspc, \
								  load_success_ratios
from plot_formats import fig_errors_vs_Kk

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_avg_successes_vs_Kk(data_flag, zero_thresh=0.1, nonzero_thresh=[0.7, 1.3],
							row_placement=[[1, 1], [1, 1], [1, 1]]):
	"""
	Args:
		data_flag: string; run identifier.
		zero_thresh: multiplier of mu_dSs for which to consider a false
			negative or false positive for the estimated odor signal.
		nonzero_thresh: 2-element list; range for which |s_est/s| is 
			considered correctly intensity decoded
		row_placement: 3-element list of 2-element lists. Entries are 
			whether to draw tick labels for x and y-axes for intensity, 
			identity, and full avg_successes plots.
		
	"""
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	
	
	assert len(iter_vars) == 4, "Need 3 iter_vars"
	iter_var_names = ['mu_Ss0', 'seed_dSs', 'Kk_1', 'Kk_2']
	for iName, name in enumerate(iter_var_names):
		assert list(iter_vars.keys())[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	mu_Ss0_vals = iter_vars['mu_Ss0']
	Kk_1_vals = iter_vars['Kk_1']
	Kk_2_vals = iter_vars['Kk_2']
	Nn = list_dict['params']['Nn']
	
	x = mu_Ss0_vals
	y = Kk_1_vals
	X, Y = sp.meshgrid(x, y)
	successes = load_success_ratios(data_flag)
	avg_successes = 100*sp.average(successes, axis=1)
	
	vminmax = [0, 100]
	ticks = [0, 50, 100]
	tick_labels = ['0', '50', '100']
	
	if row_placement[2][0] == 0:
		plt.xticks([])
	if row_placement[2][1] == 0:
		plt.yticks([])
	
	for iKk_2, Kk_2 in enumerate(Kk_2_vals):
		fig = fig_errors_vs_Kk()
		plt.xlim(10**0, 10**4)
		plt.ylim(1, 5)
		plt.xscale('log')
		plt.pcolormesh(X, Y, avg_successes[...,iKk_2].T, cmap=plt.cm.hot,
						rasterized=True, #shading='gouraud', 
						vmin=vminmax[0], vmax=vminmax[1])
		save_fig_no_whtspc('avg_successes_bkgrnd_complexity=%s' % Kk_2, \
							subdir=data_flag)
		
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
		save_fig('avg_successes_colorbar_bkgrnd_complexity=%s' % Kk_2, 
					subdir=data_flag, tight_layout=False)
		
	
			
if __name__ == '__main__':
	data_flag = get_flag()
	plot_avg_successes_vs_Kk(data_flag)