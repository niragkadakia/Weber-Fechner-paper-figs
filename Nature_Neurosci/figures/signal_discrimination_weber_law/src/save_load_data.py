"""
Functions for loading data for paper figures.

Created by Nirag Kadakia at 12:00 10-30-2017
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
	
def save_signal_discrimination_weber_law(successes,  successes_2, data_flag):
	"""
	Save list of successes based on decoding error of CS
	objects.
	
	Args:
		successes: numpy array of number of binary data for
				success (1) or not success (0), for full CS
				object array.
		successes_2: numpy array of number of binary data for
				success (1) or not success (0), for full CS
				object array, for fluctuating odor 2, indices
				given by .idxs_2
		data_flag: Data identifier for loading and saving.
	"""
	
	out_dir = '%s/analysis/%s' % (DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/signal_discrimination_weber_law.npz' % out_dir
	sp.savez(filename, successes=successes, successes_2=successes_2)
	print ('\nDecoding error data file saved to %s' % filename)

	
def load_signal_discrimination_weber_law(data_flag):
	"""
	Load list of successes, epsilons, and gains based on decoding 
	error of CS objects.
	
	Args:
		data_flag: Data identifier for loading and saving.
		
	Returns:
		data: dictionary with successes (numpy object of success ratio 
				data), gains, and epsilons.
	"""
	
	filename = '%s/analysis/%s/signal_discrimination_weber_law.npz' \
				% (DATA_DIR, data_flag)
	data = sp.load(filename)
	
	return data
	
	
def save_discrimination_accuracy_fig(fig, data_flag):

	"""
	Save firing rate subfigures.
	
	TODO
	"""
	
	out_dir = '%s/figures/signal_discrimination/%s' % (ANALYSIS_DIR, data_flag)
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
	
	
def save_sample_estimation_fig(fig, data_flag):

	"""
	Save firing rate subfigures.
	
	TODO
	"""
	
	out_dir = '%s/figures/signal_discrimination/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'sample_estimation' 
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