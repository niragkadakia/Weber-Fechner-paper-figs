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
params = {'text.usetex': False, 'mathtext.fontset': 'dejavusans'}
plt.rcParams.update(params)
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from local_methods import def_data_dir


DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')

					
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
	plt.xticks([10**-2, 10**0, 10**2], fontsize=tick_label_size)
	plt.yticks([0, 50, 100], fontsize=tick_label_size)
	plt.xlim(1e-2, 1e2)
	plt.ylim(-2, 102)
	
	return fig

	
def divisive_normalization_subfigures():
	""" 
	Generate the figure frame for the accuracy plots. Axis labels not 
	included here; add in inkscape.
	"""

	plot_size = 3
	tick_label_size = 16
	
	fig = plt.figure()
	fig.set_size_inches(plot_size, plot_size)	
	plt.xscale('log')
	plt.xticks([10**-2, 10**0, 10**2], fontsize=tick_label_size)
	plt.yticks([0, 50, 100], fontsize=tick_label_size)
	plt.xlim(1e-2, 1e2)
	plt.ylim(-2, 102)
	
	return fig