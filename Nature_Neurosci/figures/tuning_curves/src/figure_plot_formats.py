"""
Functions for generating plot formats for figures

Created by Nirag Kadakia at 11:13 12-07-2017
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
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
rc('text', usetex=True)
import matplotlib.pyplot as plt
from local_methods import def_data_dir


DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')


def tuning_curve_subfigures():
	""" 
	Generate the figure frame for tuning curve plots
	
	Returns:
		fig: The figure object.
	"""

	fig = plt.figure()
	fig.set_size_inches(2, 1.5)
	plt.ylim(-0.08, 1.0)
	plt.yticks([])
	plt.xticks([])

	return fig

	
def Kk2_subfigure():
	""" 
	Generate the figure frame for tuning curve plots
	
	Returns:
		fig: The figure object.
	"""

	fig = plt.figure()
	fig.set_size_inches(4, 4)
	plt.yticks([])
	plt.xticks([])
	
	return fig

	
def firing_rate_subfigure():
	""" 
	Generate the figure frame for tuning curve plots
	
	Returns:
		fig: The figure object.
	"""

	fig = plt.figure()
	fig.set_size_inches(4, 1.5)
	plt.xticks([])
	plt.yticks([0, 1.0], fontsize=20)
	plt.ylim(-0.025, 1.05)
	
	return fig
	
	
def firing_rate_stimulus_subfigure():
	""" 
	Generate the figure frame for tuning curve plots
	
	Returns:
		fig: The figure object.
	"""

	fig = plt.figure()
	plt.yticks([])
	plt.xticks([])
	
	return fig