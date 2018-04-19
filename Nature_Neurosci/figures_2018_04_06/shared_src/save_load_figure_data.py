"""
Functions for saving data for later analysation

Created by Nirag Kadakia at 23:30 08-02-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import os
import matplotlib.pyplot as plt
from local_methods import def_data_dir, def_figure_analysis_data_dir
import gzip, cPickle


DATA_DIR = def_data_dir()
FIGURE_ANALYSIS_DATA_DIR = def_figure_analysis_data_dir()


def save_MSE_errors(errors_nonzero, errors_zero, data_flag):
	"""
	Save decoding error from array of CS objects as numpy object.

	Args:
		errors: Error array to be saved
		data_flag: Data identifier for saving and loading.
	"""

	out_dir = '%s/%s' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/MSE_errors.npz' % out_dir
	sp.savez(filename, errors_nonzero=errors_nonzero, errors_zero=errors_zero)
	print ('\nSignal errors file saved to %s' % filename)
			
def save_binary_errors(errors_nonzero, errors_zero, data_flag, 
						errors_zero_2=None, errors_nonzero_2=None):
	"""
	Save decoding error from array of CS objects as numpy object, 
	above or below certain threshold for nonzero and zero components.

	Args:
		errors_nonzero_components: Error array to be saved, for 
			nonzero components
		errors_zero_components: Error array to be saved, for 
			zero components
		data_flag: Data identifier for saving and loading.
	"""

	out_dir = '%s/%s' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/binary_errors.npz' % out_dir
	
	if (errors_nonzero_2 is not None) and (errors_zero_2 is not None):
		sp.savez(filename, errors_nonzero=errors_nonzero, 
					errors_zero=errors_zero, errors_nonzero_2=errors_nonzero_2,
					errors_zero_2=errors_zero_2)
	else:
		sp.savez(filename, errors_nonzero=errors_nonzero, 
					errors_zero=errors_zero)
	print ('\nSignal errors file saved to %s' % filename)
	
def save_success_ratios(successes, data_flag):
	"""
	Save list of successes based on decoding error of CS
	objects.
	
	Args:
		successes: numpy array of number of binary data for
					success (1) or not success (0), for full CS
					object array.
		data_flag: Data identifier for loading and saving.
	"""
	
	out_dir = '%s/%s' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/successes.npz' % out_dir
	sp.savez(filename, successes=successes)
	print ('\nSignal binary successes file saved to %s' % filename)
	
def save_odor_ID_errors(errors, data_flag):
	"""
	Save decoding error of odor ID; by thresholding nonzeros above small value 
	and zeros below small value.
	
	Args:
		errors: Error array to be saved.
		data_flag: Data identifier for saving and loading.
	"""

	out_dir = '%s/%s' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/odor_ID_errors.npz' % out_dir
	
	sp.savez(filename, odor_ID_errors=errors)
	print ('\nOdor ID errors file saved to %s' % filename)
	
def load_binary_errors(data_flag):
	"""
	Load .npz file containing error data.

	Args:
		data_flag: Data identifier for loading and saving.
	
	Returns:
		errors: Numpy object containing error data.
	"""
	
	filename = '%s/%s/binary_errors.npz' \
				% (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	errors = dict()
	infile = sp.load(filename)
	errors['errors_nonzero'] = infile['errors_nonzero']
	errors['errors_zero'] = infile['errors_zero']
	
	return errors	
	
def load_MSE_errors(data_flag):
	"""
	Load .npz file containing error data.

	Args:
		data_flag: Data identifier for loading and saving.
	
	Returns:
		errors: Numpy object containing error data.
	"""

	filename = '%s/%s/MSE_errors.npz' \
				% (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	errors = dict()
	infile = sp.load(filename)
	errors['errors_nonzero'] = infile['errors_nonzero']
	errors['errors_zero'] = infile['errors_zero']
	
	return errors

def load_success_ratios(data_flag):
	"""
	Load list of successes based on decoding error of CS
	objects.
	
	Args:
		data_flag: Data identifier for loading and saving.
		
	Returns:
		successes: numpy object of success ratio data
	"""
	
	filename = '%s/%s/successes.npz' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	infile = sp.load(filename)
	successes = infile['successes']
	
	return successes

def load_odor_ID_errors(data_flag):
	"""
	Save decoding error of odor ID; by thresholding nonzeros above small value 
	and zeros below small value.
	
	Args:
		errors: Error array to be saved.
		data_flag: Data identifier for saving and loading.
	"""

	out_dir = '%s/%s' % (FIGURE_ANALYSIS_DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/odor_ID_errors.npz' % out_dir
	
	
	infile = sp.load(filename)
	odor_ID_errors = infile['odor_ID_errors']
	
	return odor_ID_errors
			
def save_fig(fig_name, subdir=None, clear_plot=True, tight_layout=True):
	"""
	"""
	
	if dir is None:
		out_dir = '../subfigures'
	else:
		out_dir = '../subfigures/%s' % subdir
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	if tight_layout == True:
		plt.tight_layout()
	filename = '%s/%s' % (out_dir, fig_name)
	plt.savefig('%s.svg' % filename, bbox_inches = 'tight')
	plt.savefig('%s.png' % filename, bbox_inches = 'tight')
	
	if clear_plot == True:
		plt.close()