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
from save_load_figure_data import load_temporal_errors, save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file


def plot_temporal_errors(data_flag, complexities_to_plot=range(1, 6)):
	"""
	"""
	
	data_dict = load_temporal_errors(data_flag)
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
	
	successes = data_dict['zero_errors']
	avg_successes = sp.average(successes, axis=2)
	
	plt.subplot(311)
	Tt = data_dict['Tt']
	for rate in sp.arange(2, 8):
		color = 1.*rate/14+0.3
		plt.plot(avg_successes[:, rate, 3].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	
	
	successes = data_dict['nonzero_errors']
	avg_successes = sp.average(successes, axis=2)
	
	Tt = data_dict['Tt']
	plt.subplot(312)
	for rate in sp.arange(3, 7):
		color = 1.*rate/14+0.3
		plt.plot(avg_successes[:, rate, 3].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	
	successes = data_dict['success_ratios']
	avg_successes = sp.average(successes, axis=2)
	
	Tt = data_dict['Tt']
	plt.subplot(313)
	for rate in sp.arange(2, 8):
		color = 1.*rate/14+0.3
		plt.plot(avg_successes[:, rate, 2].T, color='%.2f' % color)#, aspect=10, cmap=plt.cm.hot)#, shading='gouraud')
	#plt.colorbar()
	
	
	plt.show()
	
	
	
	print successes.shape
	
	
	# Each plot is foreground complexity versus foreground intensity
	x = adapt_rate_vals*10**conc_shift
	y = Kk_1_vals
	X, Y = sp.meshgrid(x, y)
	
	# Distinct plots for each background complexity
	bkgrnd_complexities = Kk_2_vals
	
	# Plot successes, averaged over various odor identities
	# Distinct plot for each background complexity
	for iKk2, Kk2 in enumerate(bkgrnd_complexities):
	
		if Kk2 not in complexities_to_plot:
			continue
	
		fig = fig_errors_vs_Kk()
		avg_successes = sp.average(successes[:, :, :, iKk2], axis=1)
		plt.pcolormesh(X, Y, avg_successes.T, cmap=plt.cm.hot, rasterized=True,
						shading='gouraud', vmin=-0.05, vmax=1.05)
		plt.xscale('log')
		plt.xlim(10**(-8), 10**(-4))
		plt.ylim(1, 5)
		if xticks == False:
			plt.xticks([])
		if yticks == False:
			plt.yticks([])
		if plot_colorbar == True:
			cbar = plt.colorbar()
			cbar.ax.tick_params(labelsize=15) 
		save_fig('errors_vs_Kk_bkgrnd_complexity=%s' % Kk2, subdir=data_flag)
		
		
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)