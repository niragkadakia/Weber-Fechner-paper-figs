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
					
def save_tuning_curve(tuning_curve, epsilons, Kk2s, data_flag):
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
	sp.savez(filename, tuning_curve=tuning_curve, epsilons=epsilons, Kk2s=Kk2s)
	
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

	
def save_tuning_curve_fig(fig, data_flag):
	"""
	Save tuning curve figure.
	"""
	
	out_dir = '%s/figures' % (ANALYSIS_DIR)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/tuning_curves/%s.png' % (out_dir, data_flag)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '../%s.png' % data_flag
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '../%s.svg' % data_flag
	plt.savefig(filename, bbox_inches = 'tight')