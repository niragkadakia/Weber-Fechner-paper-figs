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


def fig_errors_vs_Kk():
	"""
	Plot decoding error as a function of odor stimulus intensity
	and odor complexity, in a heatmap.
	"""
		
	fig = plt.figure()
	
	fig.set_size_inches(4, 3)
	plt.xticks(sp.arange(-10, 10), fontsize=20)
	plt.yticks(sp.arange(1, 15, 2), fontsize=20)
	
	return fig

def fig_tuning_curve():
	"""
	Plot tuning curves
	"""
		
	fig = plt.figure()
	
	fig.set_size_inches(5, 3)
	plt.xticks([])
	plt.yticks([])
	plt.ylim(0, 285)		
		
	return fig
	
def fig_Kk2():
	"""
	Plot binding matrices
	"""
		
	fig = plt.figure()
	
	plt.xticks([])
	plt.yticks([])
	
	return fig
	
def fig_Kk2_hist():
	"""
	Plot binding constnats histogram
	"""
	
	fig = plt.figure()
	
	fig.set_size_inches(3, 2.5)
	plt.ylim(-3, 0)
	plt.xlim(-3, 0)
	plt.xticks(sp.arange(-4,1), fontsize=15)
	plt.yticks(sp.arange(-4, 1), fontsize=15)
	plt.xlim(-3.2, 0.2)
	plt.ylim(-3.1, 0.5)
	
	return fig
	
def fig_temporal_kernel():
	"""
	Plot binding matrices
	"""
		
	fig = plt.figure()
	fig.set_size_inches(2.5, 1.5)
	plt.xticks([])
	plt.yticks([])
	
	return fig
	
def fig_signal_estimation():
	"""
	Plot estimated signals, ordered.
	"""
		
	fig = plt.figure()
	fig.set_size_inches(4, 2.5)
	
	return fig
	