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
	plt.xticks(sp.arange(-10, 10), fontsize=16)
	plt.yticks(sp.arange(1, 15, 2), fontsize=16)
	
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
	xticks = sp.arange(-6, 0, 2)
	plt.xticks(10.**xticks, fontsize=15)
	plt.yticks(10.**xticks, fontsize=15)
	plt.xlim(10**(-3.2), 10**0.2)
	plt.ylim(10**(-3.1), 10**0.5)
	
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
	
def fig_temporal_heatmap():
	"""
	Heatmap of estimation with time as x-axis 
	and adaptation rate (or something else) as y-axis.
	"""
	
	fig = plt.figure()
	fig.set_size_inches(5.5, 2)
	
	plt.xticks(sp.arange(-10, 10, 0.5), fontsize=18)
	plt.yticks(sp.arange(1, 15, 2), fontsize=18)
	
	return fig
	
def fig_signal_trace():
	"""
	Figure of temporal signal trace.
	"""
	
	fig = plt.figure()
	fig.set_size_inches(9, 2)
	
	plt.xticks(sp.arange(0, 10), fontsize=18)
	plt.yticks(sp.arange(0, 5000, 500), fontsize=18)
	
	return fig
	
def fig_MI_trace():
	"""
	Figure of trace of mutual information.
	"""
	
	fig = plt.figure()
	fig.set_size_inches(3, 3)
	
	plt.xscale('log')
	plt.xticks([1, 100, 10000], fontsize=18)
	plt.yticks(range(0, 100, 2), fontsize=18)
	
	return fig
	
def fig_tnse():
	"""
	Figure of 2D tsne projection of high-D data
	"""
	
	fig = plt.figure()
	fig.set_size_inches(5, 4)
	
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	
	return fig
	
def fig_avg_whiff_errors():
	"""
	Figure of plots comparing errors during whiffs
	"""
	
	fig = plt.figure()
	fig.set_size_inches(4, 3)
	
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	
	return fig
	
def fig_primacy_errors_num_active():
	"""
	Figure of plots comparing errors during whiffs
	"""
	
	fig = plt.figure()
	fig.set_size_inches(3.5, 2.5)
	
	plt.xticks(fontsize=16)
	plt.yticks([0, 50, 100], fontsize=16)
	plt.ylim(-1, 101)
		
	return fig

def fig_primacy_errors_time():
	"""
	Figure of plots comparing time to active for different intensities
	"""
	
	fig = plt.figure()
	fig.set_size_inches(3, 3)
	
	plt.xscale('log')
	plt.xticks([1, 3, 10, 30, 100], ['1x', '3x', '10x', '30x', '100x'], fontsize=16)
	plt.yticks(fontsize=16)
		
	return fig
	
def fig_primacy_schematic():
	"""
	Figure of plots to schematize primacy coding 
	"""
	
	fig = plt.figure()
	fig.set_size_inches(3, 3)
	
	plt.xticks([])
	plt.yticks([])
		
	return fig

def primacy_min_active_for_accurate_heatmap():
	"""
	Plot heatmap of minimum number of active receptors
	to correctly decode a signal, as a function of 
	odor complexity and intensity
	"""
		
	fig = plt.figure()
	
	fig.set_size_inches(4, 3)
	plt.xscale('log')
	plt.xticks([1, 3, 10, 30, 100], ['1x', '3x', '10x', '30x', '100x'], fontsize=16)
	plt.yticks(sp.arange(1, 15, 2), fontsize=16)
	
	return fig
	
def fig_classification_accuracy():
	"""
	Figure of plots to compare classification accuracy of odors in full network 
	"""
	
	fig = plt.figure()
	fig.set_size_inches(4, 4)
		
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
		
	plt.xscale('log')
	plt.ylim(0.5, 1.002)
	plt.xlim(1, 5000)
		
	return fig
	
def primacy_set_overlap():
	"""
	Figure showing % overlap of primacy sets.
	"""
	
	fig = plt.figure()
	fig.set_size_inches(3.5, 3.5)
	
	return fig