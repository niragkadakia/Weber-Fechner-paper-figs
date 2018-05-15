"""
Analysis of primacy coding idea

Created by Nirag Kadakia at 10:16 05-10-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									load_activities, save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_signal_trace_from_file


def primacy_coding(data_flag):
	"""
	"""

	success = load_success_ratios(data_flag)
	binary_errors = load_binary_errors(data_flag)
	errors_nonzero = binary_errors['errors_nonzero']
	errors_zero = binary_errors['errors_zero']
	Yy = load_activities(data_flag)
	
	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	iter_var_names = ['temporal_adaptation_rate', 'seed_dSs', 'Kk', 
						'signal_trace_multiplier']
	for iName, name in enumerate(iter_var_names):
		assert iter_vars.keys()[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	
	cmap = plt.cm.Reds
	print Yy.shape
	for rate in range(iter_vars_dims[0]):
		for Kk in [1]:#range(iter_vars_dims[2]):
			for signal_intensity in [2]:#in range(iter_vars_dims[3]):
				
				
				# TODO , this is just one dSs
				activated_receptors = (Yy[:, rate, :, Kk, signal_intensity, :] > 7.)
				pct_activated = sp.mean(activated_receptors, axis=-1)
				errors = errors_nonzero[:, rate, :, Kk, signal_intensity]
				
				title='%s, %s, %s' % (rate, Kk, signal_intensity)
				plt.title(title)
				plt.xlim(0, 1)
				for idSs in range(iter_vars_dims[1]):
					plt.subplot(5, 4, idSs + 1)
					color = cmap(0.4 + 0.5*idSs/(iter_vars_dims[1]))
					data = sp.vstack((pct_activated[:, idSs].T, errors[:, idSs].T)).T
					sorted_data = data[sp.argsort(data[:, 0])]
					plt.plot(sorted_data[:, 0], sorted_data[:, 1], color=color)
					plt.xlim(0, 1)
					plt.xticks([])
					plt.yticks([])
				plt.show()
	
				
				title='%s, %s, %s' % (rate, Kk, signal_intensity)
				plt.title(title)
				activities = Yy[:, rate, :, Kk, signal_intensity, :]
				times = [491, 495, 500, 505, 509]
				for iT, time in enumerate(times):
					color = plt.cm.seismic(0.2 + 0.6*iT/(len(times) - 1))
					act = sp.zeros(35)
					for idSs in range(iter_vars_dims[1]):
						act += activated_receptors[time, idSs]
					plt.plot(sp.sort(act), color=color)
					plt.ylim(0, 20)
					
				plt.show()
	
				
				
if __name__ == '__main__':
	data_flag = get_flag()
	primacy_coding(data_flag)