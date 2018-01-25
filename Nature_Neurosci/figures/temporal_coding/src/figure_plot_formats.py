"""
Functions for generating plot formats for figures

Created by Nirag Kadakia at 15:03 01-24-2018
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
params = {'text.usetex': False, 'mathtext.fontset': 'dejavusans'}
plt.rcParams.update(params)
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from local_methods import def_data_dir


DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')

					
def signal_trace_subfigures(xlims):
	"""
	"""

	plot_height = 3
	plot_width = 25*(xlims[1] - xlims[0])
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)	
	plt.xticks(fontsize=tick_label_size)
	plt.yticks(sp.arange(0, 5, 0.2), fontsize=tick_label_size)
	
	return fig
	
def epsilon_trace_subfigures(xlims):
	"""
	"""

	plot_height = 3
	plot_width = 25*(xlims[1] - xlims[0])
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)	
	plt.xticks(fontsize=tick_label_size)
	plt.yticks(sp.arange(-10, 10, 1), fontsize=tick_label_size)
	
	return fig
	
def errors_trace_subfigures(xlims):
	"""
	"""

	plot_height = 3
	plot_width = 25*(xlims[1] - xlims[0])
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)	
	plt.xticks(fontsize=tick_label_size)
	plt.yticks(sp.arange(0, 150, 25), fontsize=tick_label_size)
	
	plt.ylim(-2, 102)
	
	return fig
	
def est_signal_zeros_subfigures(ylims):
	"""
	"""
	
	plot_height = 3
	plot_width = 3
	tick_label_size = 16
	yticks = [ylims[0], 0, ylims[1]]
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)	
	plt.xticks([0])
	plt.yticks(yticks, fontsize=tick_label_size)
	plt.ylim(ylims[0], ylims[1])
		
	return fig
	
def est_signal_nonzeros_subfigures(ylims):
	"""
	"""
	
	plot_height = 3
	plot_width = 3
	tick_label_size = 16
	yticks = [ylims[0], 0, ylims[1]]
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)	
	plt.xticks([0])
	plt.yticks(yticks, fontsize=tick_label_size)
	plt.ylim(ylims[0], ylims[1])
		
	return fig
	