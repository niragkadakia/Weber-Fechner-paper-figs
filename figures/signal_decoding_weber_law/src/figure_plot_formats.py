"""
Functions for generating plot formats for figures

Created by Nirag Kadakia at 10:13 10-05-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import matplotlib
from matplotlib import cm
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from local_methods import def_data_dir


DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')

					
def single_decoding_weber_law_plot(data_idxs=2):
	""" 
	Generate the figure frame for tuning curve plots
	
	Args: 
		plot_vars: numpy array with 2 columns giving 
			iterated variable indices to plot.
		params: dictionary of parameters pulled from specs file
				
	Returns:
		fig: The figure object.
	"""

	decode_plot_size = 8
	num_gain_plots = data_idxs / 2
	gain_plot_size = 4
	est_plot_height = 4
	num_est_plots = 2
	
	fig = plt.figure()
	grid_width = decode_plot_size + 1 + gain_plot_size*num_gain_plots
	grid_height = decode_plot_size + 1 + est_plot_height*num_est_plots
	gs = gridspec.GridSpec(grid_height, grid_width)
	gs.update(wspace=1.3, hspace=1.0)
	fig.set_size_inches(12, 12)
	
	ax = dict()
	
	# Success error plots
	ax['successes'] = plt.subplot(gs[0:decode_plot_size, 0:decode_plot_size])
	plt.xscale('log')
	plt.xticks(10.**sp.arange(-5, 5),  fontsize=14)
	plt.yticks([0, 50, 100], fontsize=14)
	ax['successes'].set_xlabel(r'Background signal', fontsize=18)
	ax['successes'].set_ylabel(r'Correctly decoded signals (\%)', fontsize=18)
	ax['successes'].set_xlim([1e-1, 1e2])
	
	# Gain plots
	for data_idx in sp.arange(0, data_idxs):
		if data_idx % 2 == 0: 
			ax['gains_%s' % data_idx] = \
				plt.subplot(gs[0:gain_plot_size, decode_plot_size + \
							1 + (data_idx/2)*gain_plot_size: decode_plot_size + \
							1 + gain_plot_size + (data_idx/2)*gain_plot_size])
		elif data_idx % 2 == 1:
			ax['gains_%s' % data_idx] = \
				plt.subplot(gs[gain_plot_size:2*gain_plot_size, decode_plot_size + \
							1 + (data_idx/2)*gain_plot_size: decode_plot_size  + \
							1 + gain_plot_size + (data_idx/2)*gain_plot_size])
	
		ax['gains_%s' % data_idx].set_xscale('log')
		ax['gains_%s' % data_idx].set_yscale('log')
		ax['gains_%s' % data_idx].set_xlim(1e-2, 1e0)
		ax['gains_%s' % data_idx].set_xticks([])
		ax['gains_%s' % data_idx].set_yticks([])
		
		if data_idx / 2 == 0:
			ax['gains_%s' % data_idx].set_yticks(10.**sp.arange(-5, 0))
			ax['gains_%s' % data_idx].set_xticks([1e-2, 1e-1, 1e0])	
			ax['gains_%s' % data_idx].tick_params(labelsize=14)

		if data_idx % 2 == 0:
			ax['gains_%s' % data_idx].set_ylim(3e-3, 4e-2)
		elif data_idx % 2 == 1:
			ax['gains_%s' % data_idx].set_ylim(1e-2, 2e-1)
		ax['gains_%s' % data_idx].set_xlim(0.98e-2, 1.1e-1)
			
	# Sample estimation plot
	for est_idx in range(num_est_plots):
		ax['est_%s' % est_idx] = \
			plt.subplot(gs[decode_plot_size + 1 + est_idx*est_plot_height:
				decode_plot_size + 1 + (est_idx + 1)*est_plot_height, 
				0:grid_width])
		ax['est_%s' % est_idx].set_xticks([])
		ax['est_%s' % est_idx].set_yticks([])
	
	return fig, ax

