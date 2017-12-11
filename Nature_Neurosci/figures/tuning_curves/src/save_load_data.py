"""
Functions for loading data for paper figures.

Created by Nirag Kadakia at 20:30 10-04-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import cPickle
import shelve
import gzip
import os
from local_methods import def_data_dir, def_analysis_dir
import matplotlib.pyplot as plt

DATA_DIR = def_data_dir()
ANALYSIS_DIR = def_analysis_dir()
					
def save_tuning_curve(tuning_curve, epsilons, Kk2s, activities, data_flag):
	"""
	Save objects saved by save_objects in the save_data module. 

	Args:
		tuning_curves: numpy array holding response data.
		data_flag: Data identifier for loading and saving.
	"""

	out_dir = '%s/objects/%s' % (DATA_DIR, data_flag)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	filename = '%s/tuning_curve' % out_dir
	sp.savez(filename, tuning_curve=tuning_curve, epsilons=epsilons, 
				Kk2s=Kk2s, activities=activities)
	
	print ('Tuning curve data saved to %s' % filename) 
					
					
def load_tuning_curve(data_flag):
	"""
	Load tuning response curve data, saved as npz file.
	
	Args: 
		data_flag: Data identifier for loading and saving.
	
	Returns:
		tuning_curves: Dictionary of keyed items and their respective values.
	"""
	
	out_dir = '%s/objects/%s' % (DATA_DIR, data_flag)
	tuning_curve_data = sp.load('%s/tuning_curve.npz' % out_dir)
	
	return tuning_curve_data

	
def save_tuning_curve_fig(fig, mu_dSs_idx, sigma_Kk2_idx, fig_num, data_flag):
	"""
	Save tuning curve subfigures.
	
	Args:
		fig: figure object to save
		mu_dSs_idx: index of mu_dSs (among range of first iter_var in specs)
		sigma_Kk2_idx: index of Kk2 standard deviation (in range of second 
						iter_var in specs)
		fig_num: number of tuning curve
		data_flag: data identifier
	"""
	
	out_dir = '%s/figures/tuning_curves/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'tuning_curve_mu_dSs_idx=%s,sigma_Kk2_idx=%s,fig_num=%s' \
				% (mu_dSs_idx, sigma_Kk2_idx, fig_num)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	

def save_Kk2_fig(fig, sigma_Kk2_idx, data_flag):
	"""
	Save Kk2 matrix subfigure.
	
	Args:
		fig: figure object to save
		sigma_Kk2_idx: index of Kk2 standard deviation (in range of second 
						iter_var in specs)
		data_flag: data identifier
	"""
	
	out_dir = '%s/figures/tuning_curves/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'Kk2_sigma_Kk2_idx=%s' % sigma_Kk2_idx
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	

def save_firing_rate_fig(fig, sigma_Kk2_idx, mu_dSs_idxs, 
							odor_seed, data_flag):

	"""
	Save firing rate subfigures.
	
	Args:
		fig: figure object to save
		sigma_Kk2_idx: index of Kk2 standard deviation (in range of second 
						iter_var in specs)
		mu_dSs_idxs: indices of mu_dSs for which stimulus is plotted 
						(among range of first iter_var in specs)
		odor_seed: random number seed for odor stimulus plotted in this figure
		data_flag: data identifier
	"""
	
	out_dir = '%s/figures/tuning_curves/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'firing_rate_sigma_Kk2_idx=%s,mu_dSs_idxs=%s,odor_seed=%s,' \
				% (sigma_Kk2_idx, mu_dSs_idxs, odor_seed)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	
def save_firing_rate_stimulus_fig(fig, sigma_Kk2_idx, mu_dSs_idxs, 
							odor_seed, data_flag):

	"""
	Save firing rate stimuli figures.
	
	Args:
		fig: figure object to save
		sigma_Kk2_idx: index of Kk2 standard deviation (in range of second 
						iter_var in specs)
		mu_dSs_idxs: indices of mu_dSs for which stimulus is plotted 
						(among range of first iter_var in specs)
		odor_seed: random number seed for odor stimulus plotted in this figure
		data_flag: data identifier
	"""
	
	out_dir = '%s/figures/tuning_curves/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'firing_rate_stimulus_sigma_Kk2_idx=%s,mu_dSs_idxs=%s,' \
				'odor_seed=%s,' % (sigma_Kk2_idx, mu_dSs_idxs, odor_seed)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')