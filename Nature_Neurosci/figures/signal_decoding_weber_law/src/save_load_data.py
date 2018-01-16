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
from local_methods import def_data_dir
from local_methods import def_analysis_dir
import matplotlib.pyplot as plt

DATA_DIR = def_data_dir()
ANALYSIS_DIR = def_analysis_dir()		

def load_aggregated_object_list(iter_vars_dims, data_flag):
	"""
	Load objects saved by aggregated_objects module; in correct array
	reflecting the shape of the iterated variables.

	Args:
		data_flag: Data identifier for loading and saving.
	
	Returns:
		agg_obj: Array of objects
	"""

	filename = '%s/objects/%s/aggregated_objects.pklz' % (DATA_DIR, data_flag)
	with gzip.open(filename, 'rb') as f:
		CS_object_list = cPickle.load(f)
	
	CS_object_array = sp.reshape(CS_object_list, iter_vars_dims)

	return CS_object_array
	
def save_signal_decoding_weber_law(data, data_flag):
	"""
	Save list of successes based on decoding error of CS
	objects.
	
	Args:
		successes: numpy array of number of binary data for
				success (1) or not success (0), for full CS
				object array.
		gains: numpy array of gains
		epsilons: numpy aray of free energy 
		activities: numpy array for dYy
		Kk2s: numpy array for Kk2s
		data_flag: Data identifier for loading and saving.
	"""
	
	out_dir = '%s/analysis/%s' % (DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/signal_decoding_weber_law.npy' % out_dir
	sp.save(filename, data)
	print ('\nDecoding error data file saved to %s' % filename)

def load_signal_decoding_weber_law(data_flag):
	"""
	Load list of successes, epsilons, and gains based on decoding 
	error of CS objects.
	
	Args:
		data_flag: Data identifier for loading and saving.
		
	Returns:
		data: dictionary with successes (numpy object of success ratio 
				data), gains, and epsilons.
	"""
	
	filename = '%s/analysis/%s/signal_decoding_weber_law.npy' \
				% (DATA_DIR, data_flag)
	data = sp.load(filename).item()
	
	return data
	
def save_decoding_accuracy_fig(fig, data_flag):

	"""
	Save decoding accuracy subfigures.
	"""
	
	out_dir = '%s/figures/signal_decoding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'decoding_accuracy' 
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
	
def save_activities_fig(fig, data_flag):

	"""
	Save activities subfigures.
	"""
	
	out_dir = '%s/figures/signal_decoding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'activities' 
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
	plt.close()
	
def save_Kk2_fig(fig, data_flag):

	"""
	Save activities subfigures.
	"""
	
	out_dir = '%s/figures/signal_decoding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'Kk2' 
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
	
def save_signal_estimation_zeros_fig(fig, data_flag, mu_dSs, seed_Kk2):

	"""
	Save activities subfigures.
	"""
	
	out_dir = '%s/figures/signal_decoding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'signal_estimation_zeros_%s,%s' % (mu_dSs, seed_Kk2) 
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
	
	
def save_signal_estimation_nonzeros_fig(fig, data_flag, mu_dSs, seed_Kk2):

	"""
	Save activities subfigures.
	"""
	
	out_dir = '%s/figures/signal_decoding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'signal_estimation_nonzeros_%s,%s' % (mu_dSs, seed_Kk2) 
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
	
	
	
	
