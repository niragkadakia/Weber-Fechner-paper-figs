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
from local_methods import def_data_dir

DATA_DIR = def_data_dir()
VAR_STRINGS = dict(mu_Ss0 = '$\langle s_0 \\rangle$')

					
def tuning_curve_plot_epsilon(plot_vars, iter_vars, params):
	""" 
	Generate the figure frame for tuning curve plots
	
	Args: 
		plot_vars: numpy array with 2 columns giving 
			iterated variable indices to plot.
		params: dictionary of parameters pulled from specs file
				
	Returns:
		fig: The figure object.
	"""
	
	plot_dims = [len(plot_vars[0]), len(plot_vars[1])]
	
	assert (len(plot_vars[0]) % 2 == 1), \
		'Best to have an odd number of columns'
	assert (len(plot_vars[1]) % 2 == 1), \
		'Best to have an odd number of rows'
	
	fig = plt.figure()
	fig.set_size_inches(8, 8)
	axes_tuning = sp.zeros(plot_dims, dtype='object')
	axes_Kk2 = sp.zeros(plot_dims[1], dtype='object')
	axes_signal = sp.zeros(plot_dims, dtype='object')
	
	gs = gridspec.GridSpec(4*plot_dims[0], 4*plot_dims[1] +  5)
	
	for idx in range(plot_dims[0]):
			for idy in range(plot_dims[1]):
				
				# Tuning curve frames
				for iM in range(params['Mm']):
					axes_tuning[idx, idy] = \
						plt.subplot(gs[idy*4:(idy + 1)*4, idx*4:(idx + 1)*4])
					plt.ylim(0, 1.05)
					plt.yticks([])
					plt.xticks([])
				
				# Ylabel in center plot
				if idx == 0:
					if idy == len(plot_vars[1]) / 2:
						axes_tuning[idx, idy].set_ylabel(r'activity', 
														fontsize=19)
					plt.yticks([0, 1], fontsize = 15)
					
				# Xlabel in center plot
				if idy == plot_dims[1] - 1:
					if idx == len(plot_vars[0]) / 2:
						axes_tuning[idx, idy].set_xlabel(r'Odorant' \
														'$(i = 1,...,100)$', \
														fontsize=15)
				
				
				
				# Epsilon frame
				if idx == 0:
					axes_Kk2[idy] = plt.subplot(gs[idy*4:(idy + 1)*4, 
										4*plot_dims[0]+1:plot_dims[0]*4 + 4])
					
					
					axes_Kk2[idy].yaxis.set_label_position('right')
					plt.yticks(sp.arange(0, 50, 10))
					plt.ylim(3, 25)
					axes_Kk2[idy].yaxis.tick_left()
					axes_Kk2[idy].tick_params(axis='y', labelsize=14)
					plt.xticks([])
					
					# Ylabel in center plot
					if idy == len(plot_vars[1]) / 2:
						axes_Kk2[idy].set_ylabel(r'$\epsilon^\rho_' \
												'{\\textup {u}}$', 
												fontsize=19, rotation=0, 
												labelpad=16)
					
					if idy == plot_dims[1] - 1:
						axes_Kk2[idy].set_xlabel(r' Receptor', 
													fontsize=15)
	
				# Side text
				if idx == 0:
					if idy == 0:
						axes_tuning[idx, idy].text(0, 1.20, 
								r'{\textbf{\textit{Increasing signal' \
									' background}}} $\\rightarrow$', 
									fontsize=16, color='green')
				
					if idy == plot_dims[1] - 1:
						axes_tuning[idx, idy].text(-160, 2.1, 
								r'{\textbf{\textit{Increasing receptor' \
									' diversity}}} $\\rightarrow$', 
									fontsize=18, rotation=90, color='green')
								
	
	return fig, plot_dims, axes_tuning, axes_Kk2, axes_signal


def tuning_curve_plot_Kk2(plot_vars, iter_vars, params):
	""" 
	Generate the figure frame for tuning curve plots
	
	Args: 
		plot_vars: numpy array with 2 columns giving 
			iterated variable indices to plot.
		params: dictionary of parameters pulled from specs file
				
	Returns:
		fig: The figure object.
	"""
	
	plot_dims = [len(plot_vars[0]), len(plot_vars[1])]
	
	assert (len(plot_vars[0]) % 2 == 1), \
		'Best to have an odd number of columns'
	assert (len(plot_vars[1]) % 2 == 1), \
		'Best to have an odd number of rows'
	
	fig = plt.figure()
	fig.set_size_inches(8, 8)
	axes_tuning = sp.zeros(plot_dims, dtype='object')
	axes_Kk2 = sp.zeros(plot_dims[1], dtype='object')
	axes_signal = sp.zeros(plot_dims, dtype='object')
	
	gs = gridspec.GridSpec(4*plot_dims[0], 4*plot_dims[1] +  5)
	gs.update(wspace=1.0, hspace=1.0)
	
	for idx in range(plot_dims[0]):
			for idy in range(plot_dims[1]):
				
				# Tuning curve frames
				for iM in range(params['Mm']):
					axes_tuning[idx, idy] = \
						plt.subplot(gs[idy*4:(idy + 1)*4, idx*4:(idx + 1)*4])
					plt.ylim(0, 1.05)
					plt.yticks([])
					plt.xticks([])
				
				# Ylabel in center plot
				if idx == 0:
					if idy == len(plot_vars[1]) / 2:
						axes_tuning[idx, idy].set_ylabel(r'Activity', 
														fontsize=17)
					plt.yticks([0, 1], fontsize = 15)
					
				# Xlabel in center plot
				if idy == plot_dims[1] - 1:
					if idx == len(plot_vars[0]) / 2:
						axes_tuning[idx, idy].set_xlabel(r'Odorant' \
														'$(i = 1,...,200)$', \
														fontsize=17, labelpad=10)
				
				
				
				# Kk2 frame
				if idx == 0:
					axes_Kk2[idy] = plt.subplot(gs[idy*4:(idy + 1)*4, 
										4*plot_dims[0]:plot_dims[0]*4 + 5])
					axes_Kk2[idy].yaxis.set_label_position('right')
					plt.xticks([])
					plt.yticks([])
					
					# Ylabel in center plot
					if idy == len(plot_vars[1]) / 2:
						axes_Kk2[idy].set_ylabel(r'Odorant $(i = 1,...,200)$',
												fontsize=17, rotation=270, 
												labelpad=26)
					
					if idy == plot_dims[1] - 1:
						axes_Kk2[idy].set_xlabel(r' Receptor', 
													fontsize=17, labelpad=10)
	
				# Side text
				if idx == 0:
					if idy == 0:
						axes_tuning[idx, idy].text(0, 1.20, 
								r'{\textbf{\textit{Increasing signal' \
									' background}}} $\\rightarrow$', 
									fontsize=16, color='green')
				
					if idy == plot_dims[1] - 1:
						axes_tuning[idx, idy].text(-160, 2.1, 
								r'{\textbf{\textit{Increasing receptor' \
									' diversity}}} $\\rightarrow$', 
									fontsize=18, rotation=90, color='green')
								
	
	return fig, plot_dims, axes_tuning, axes_Kk2, axes_signal

