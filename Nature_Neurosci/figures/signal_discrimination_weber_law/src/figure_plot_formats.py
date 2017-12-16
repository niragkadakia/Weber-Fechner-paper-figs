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
#matplotlib.rcParams['font.sans-serif'] = "Helvetica"
#matplotlib.rcParams['font.family'] = "sans-serif"
#rc('text', usetex=True)
import matplotlib.pyplot as plt
params = {'text.usetex': False, 'mathtext.fontset': 'dejavusans'}
plt.rcParams.update(params)
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from local_methods import def_data_dir


DATA_DIR = def_data_dir()

					
def decoding_accuracy_subfigures():
	""" 
	Generate the figure frame for the accuracy plots. Axis labels not 
	included here; add in inkscape.
	"""

	plot_size = 4
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_size, plot_size)	
	plt.xscale('log')
	plt.xticks(10.**sp.arange(-5, 5), fontsize=tick_label_size)
	plt.yticks([0, 50, 100], fontsize=tick_label_size)
	plt.xlim(1e-1, 1e2)
	
	return fig
	
	
def sample_estimation_subfigures(plot_xlims=[20, 60], inset_xlims=[0, 100],
									inset_ylims=[-0.14, 0.14], 
									inset_loc=1, inset_mag=3.0):
	""" 
	Generate the figure frame for the accuracy plots. Axis labels not 
	included here; add in inkscape.
	
	Args:
		plot_xlims: 2-element list; lower and upper x-axis limits of main plot
		inset_xlims: 2-element list; lower and upper x-axis limits of inset
		inset_ylims: 2-element list; lower and upper y-axis limits of inset
		inset_loc: integer 1-4; location of inset
		inset_mag: float; magnification of inset
	"""

	plot_height = 3
	plot_width = 5
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_width, plot_height)
	ax = plt.subplot(111)
	
	plt.xticks([0, 20, 40, 60, 80, 100], fontsize=tick_label_size)
	plt.yticks([])
	plt.xlim(plot_xlims[0], plot_xlims[1])

	# Signal insert (zoom-in) plots 
	ax_insert = zoomed_inset_axes(ax, inset_mag, loc=inset_loc)
	ax_insert.set_xlim(inset_xlims[0], inset_xlims[1])
	ax_insert.set_ylim(inset_ylims[0], inset_ylims[1])
	ax_insert.set_xticks([])
	ax_insert.set_yticks([])
	ax_insert.patch.set_alpha(0.80)
	ax_insert.set_facecolor('0.97')
	
	return fig, ax, ax_insert