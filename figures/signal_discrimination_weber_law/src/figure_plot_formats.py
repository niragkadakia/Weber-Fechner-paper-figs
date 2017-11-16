"""
Functions for generating plot formats for figures

Created by Nirag Kadakia at 10:13 10-30-2017
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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from local_methods import def_data_dir


DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')

					
def signal_discrimination_weber_law_plot(Kk_split_idxs=None):
	""" 
	Generate the figure frame for discrimination task plots.
	"""

	signal_plot_height = 9
	decode_plot_size = signal_plot_height
	signal_plot_width = int(decode_plot_size*5./3.)
	
	fig = plt.figure()
	grid_width = decode_plot_size + 1 + signal_plot_width
	grid_height = (decode_plot_size+1)*(Kk_split_idxs)
	gs = gridspec.GridSpec(grid_height, grid_width)
	gs.update(wspace=1, hspace=0.5)
	fig.set_size_inches(12, 5*Kk_split_idxs)
	
	# Labels and tick sizes for later
	success_x_tick_pad = 12
	success_x_label_pad = 12
	success_y_label_pad = 16
	signal_x_label_pad = 15
	signal_y_label_pad = 42
	tick_label_size = 24
	axis_label_size = 28
	
	ax = dict()
	
	for Kk_split_idx in range(Kk_split_idxs):
		
		# Success error plots
		plot_y1 = (decode_plot_size + 1)*Kk_split_idx 
		plot_y2 = (decode_plot_size + 1)*Kk_split_idx  + decode_plot_size
		plot_x1 = 0
		plot_x2 = decode_plot_size
		
		ax['successes_%s' % Kk_split_idx] = \
			plt.subplot(gs[plot_y1:plot_y2, plot_x1:plot_x2])
		plt.xscale('log')
		plt.yticks([0, 50, 100], fontsize=tick_label_size)
		plt.xticks([])
		
		# X-labels only on bottom / set xlim for all though
		if Kk_split_idx == Kk_split_idxs - 1:
			ax['successes_%s' % Kk_split_idx].\
				set_xlabel(r'Background odor''\n''strength', 
							fontsize=axis_label_size,
							labelpad=success_x_label_pad)
			plt.xticks(10.**sp.arange(-5, 5),  fontsize=tick_label_size)
		ax['successes_%s' % Kk_split_idx].tick_params(axis='x', 
							which='major', pad=success_x_tick_pad)
		ax['successes_%s' % Kk_split_idx].set_xlim([1e-1, 1e1])
		
		# Signal plots		
		plot_y1 = (decode_plot_size + 1)*Kk_split_idx 
		plot_y2 = (decode_plot_size + 1)*Kk_split_idx \
					+ signal_plot_height
		plot_x1 = decode_plot_size + 1
		plot_x2 = decode_plot_size + signal_plot_width + 1
		
		ax['signal_%s' % Kk_split_idx] = \
			plt.subplot(gs[plot_y1:plot_y2, plot_x1:plot_x2])
		plt.yticks([])
		if Kk_split_idx == Kk_split_idxs - 1:
			plt.xticks([20, 40, 60, 80], fontsize=tick_label_size)
			ax['signal_%s' % Kk_split_idx].tick_params(axis='x', pad=10)
		else:
			plt.xticks([])
		plt.xlim(18, 70)
		
		# Signal insert (zoom-in) plots 
		insert_xlims = [30, 50]
		insert_ylims = [-0.14, 0.14]
		ax['signal_insert_%s' % Kk_split_idx] = \
			zoomed_inset_axes(ax['signal_%s' % Kk_split_idx], 2.0, loc=2)
		ax['signal_insert_%s' % Kk_split_idx].\
			set_xlim(insert_xlims[0], insert_xlims[1])
		ax['signal_insert_%s' % Kk_split_idx].\
			set_ylim(insert_ylims[0], insert_ylims[1])
		ax['signal_insert_%s' % Kk_split_idx].set_xticks([])
		ax['signal_insert_%s' % Kk_split_idx].set_yticks([])
		ax['signal_insert_%s' % Kk_split_idx].patch.set_alpha(0.75)
		ax['signal_insert_%s' % Kk_split_idx].set_facecolor('0.96')
		
		# X-labels only on bottom / set xlim for all though
		if Kk_split_idx == Kk_split_idxs - 1:
			ax['signal_%s' % Kk_split_idx].\
				set_xlabel(r'Odorant identity', fontsize=axis_label_size, 
							labelpad=signal_x_label_pad)
		
		# All Y-labels: only in middle if odd number of Kk_split_idxs
		if Kk_split_idxs % 2 == 0:
			ax['successes_%s' % Kk_split_idx].\
					set_ylabel(r'Correctly decoded signals (%)', 
					labelpad=success_y_label_pad, fontsize=axis_label_size)
			ax['signal_%s' % Kk_split_idx].\
				yaxis.set_label_position("right")
			ax['signal_%s' % Kk_split_idx].\
				set_ylabel(r'Concentration', fontsize=axis_label_size, 
							rotation=270, labelpad=signal_y_label_pad)
		else:
			if Kk_split_idx == Kk_split_idxs / 2:
				ax['successes_%s' % Kk_split_idx].\
					set_ylabel(r'Correctly decoded signals (%)', 
					labelpad=success_y_label_pad, fontsize=axis_label_size)
				ax['signal_%s' % Kk_split_idx].\
					yaxis.set_label_position("right")
				ax['signal_%s' % Kk_split_idx].\
					set_ylabel(r'Concentration', fontsize=axis_label_size, 
								rotation=270, labelpad=signal_y_label_pad)
		
		
	return fig, ax