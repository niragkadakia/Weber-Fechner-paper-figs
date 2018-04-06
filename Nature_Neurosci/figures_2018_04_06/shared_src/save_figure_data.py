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
import cPickle
import shelve
import gzip
import os
import time
try:
	import matplotlib.pyplot as plt
except:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
from local_methods import def_data_dir, def_figure_analysis_data_dir


DATA_DIR = def_data_dir()
SAVE_FIGURE_ANALYSIS_DATA = def_figure_analysis_data_dir()

	
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
	
	out_dir = '%s/%s' % (SAVE_FIGURE_ANALYSIS_DATA, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	filename = '%s/successes.npz' % out_dir
	sp.savez(filename, successes=successes)
	print ('\nSignal binary successes file saved to %s' % filename)
	