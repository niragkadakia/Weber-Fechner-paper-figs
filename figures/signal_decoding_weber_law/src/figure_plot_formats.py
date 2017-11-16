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
import matplotlib.pyplot as plt
plt.rcParams['text.latex.preamble'] = [
		r'\usepackage{helvet}',
		r'\usepackage{sansmath}', # math --> arial
		r'\sansmath']             # <-- actually tell tex to use it!
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
	num_gain_plots = data_idxs/2
	gain_plot_size = 4
	est_plot_upper_pad = 3
	est_plot_height = 4
	num_est_plots = 2
	
	
	fig = plt.figure()
	grid_width = decode_plot_size + 1 + gain_plot_size*num_gain_plots
	grid_height = decode_plot_size + est_plot_upper_pad + \
					est_plot_height*num_est_plots
	gs = gridspec.GridSpec(grid_height, grid_width)
	gs.update(wspace=6.5, hspace=3.0)
	fig.set_size_inches(13, 12)
	
	# Some font and label sizes for later use
	tick_label_size_small = 20
	tick_label_size_large = 22
	axis_label_size = 28
	x_axis_label_pad = 13
	y_axis_successes_labelpad = 10
	x_tick_label_pad = 11
	
	ax = dict()
	
	# Success error plots
	ax['successes'] = plt.subplot(gs[0:decode_plot_size, 0:decode_plot_size])
	plt.xscale('log')
	plt.xticks(10.**sp.arange(-5, 5),  fontsize=tick_label_size_small)
	plt.yticks([0, 50, 100], fontsize=tick_label_size_large)
	ax['successes'].tick_params(axis='x', which='major', pad=x_tick_label_pad)
	ax['successes'].set_xlabel(r'Background signal', fontsize=axis_label_size, 
								labelpad=x_axis_label_pad)
	ax['successes'].set_ylabel(r'Correctly decoded''\n''signals ($\%$)', 
								fontsize=axis_label_size, 
								labelpad=y_axis_successes_labelpad, 
								multialignment='center')
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
		
		# Gain plots ticks labels and axis ranges
		if data_idx / 2 == 0:
			ax['gains_%s' % data_idx].set_yticks(10.**sp.arange(-5, 0))
			ax['gains_%s' % data_idx].tick_params(labelsize=tick_label_size_large)
		if data_idx % 2 == 0:
			ax['gains_%s' % data_idx].set_ylim(3e-3, 4e-2)
		elif data_idx % 2 == 1:
			ax['gains_%s' % data_idx].set_ylim(1e-2, 2e-1)
			ax['gains_%s' % data_idx].set_xticks([1e-2, 1e-1, 1e0])	
			ax['gains_%s' % data_idx].tick_params(labelsize=tick_label_size_small)
		ax['gains_%s' % data_idx].set_xlim(0.98e-2, 1.1e-1)
		ax['gains_%s' % data_idx].tick_params(axis='x', which='major', 
												pad=x_tick_label_pad)
			
		# Gaain axis label in center of range of plots
		if num_gain_plots % 2 == 0:
			if data_idx == data_idxs/2 - 1:
				ax['gains_%s' % data_idx].xaxis.set_label_coords(1.1, -0.38)
				ax['gains_%s' % data_idx].set_xlabel(r'Background signal', 
											fontsize=axis_label_size, 
											labelpad=x_axis_label_pad)
		elif num_gain_plots % 2 == 1:
			if data_idx == data_idxs/2:
				ax['gains_%s' % data_idx].set_xlabel(r'Background signal', 
											fontsize=axis_label_size, 
											labelpad=x_axis_label_pad)
		if data_idx == data_idxs - 2:
			ax['gains_%s' % data_idx].yaxis.set_label_position('right')
			ax['gains_%s' % data_idx].set_ylabel(
						r'Receptor gain''\n''$\Delta a / \Delta s$', 
						fontsize=axis_label_size, rotation=270)
			ax['gains_%s' % data_idx].yaxis.set_label_coords(1.50, -0.1)
				
	# Sample estimation plot
	for est_idx in range(num_est_plots):
		ax['est_%s' % est_idx] = \
			plt.subplot(gs[decode_plot_size + est_plot_upper_pad
						+ est_idx*est_plot_height:
						decode_plot_size + est_plot_upper_pad 
						+ (est_idx + 1)*est_plot_height, 0:grid_width])
		ax['est_%s' % est_idx].set_yticks([])
		
		# Sample estimation plots labels only on outside
		if est_idx == num_est_plots/2:
			ax['est_%s' % est_idx].set_xticks(sp.arange(0, 100, 10))
			ax['est_%s' % est_idx].yaxis.set_label_coords(-0.03, 1.1)
			ax['est_%s' % est_idx].set_xlabel(r'Odorant identity', 
												fontsize=axis_label_size,
												labelpad=x_axis_label_pad)
			ax['est_%s' % est_idx].set_ylabel(r'Odorant intensity', 
												fontsize=axis_label_size) 
			ax['est_%s' % est_idx].tick_params(labelsize=tick_label_size_large, 
												pad=x_tick_label_pad)
		else:
			ax['est_%s' % est_idx].set_xticks([])
		ax['est_%s' % est_idx].set_xlim(35, 99)
		
	return fig, ax

