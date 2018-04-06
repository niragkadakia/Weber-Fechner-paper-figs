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
	
	fig.set_size_inches(5, 4)
	plt.xticks(sp.arange(-10, 0), fontsize=13)
	plt.yticks(sp.arange(1, 20), fontsize=13)
	
	return fig

def fig_tuning_curve():
	"""
	Plot tuning curves
	"""
		
	fig = plt.figure()
	
	fig.set_size_inches(5, 4)
	plt.xticks([])
	plt.yticks([])
	
	return fig