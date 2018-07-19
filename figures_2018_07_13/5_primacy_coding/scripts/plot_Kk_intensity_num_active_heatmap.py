"""
Plot heatmap of minimum active neurons for x% accuracy, 
as a function of Kk and signal intensity.

Created by Nirag Kadakia at 10:16 07-16-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from plot_formats import primacy_min_active_for_accurate_heatmap
from save_load_figure_data import load_binary_errors, \
									load_activities, save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file


def accuracy_vs_num_active(data_flag, active_rate=20, min_accurate_pct=75):
	
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
		assert list(iter_vars.keys())[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	Mm = list_dict['params']['Mm']
	num_rates = iter_vars_dims[0]
	num_odors = iter_vars_dims[1]
	num_Kks = iter_vars_dims[2]
	num_intensities = iter_vars_dims[3]
	
	# Separate plots for adaptation rate, Kk, and signal step intensity
	iter_array_dims = [num_rates, num_Kks, num_intensities]
	it = sp.nditer(sp.zeros(iter_array_dims), flags=['multi_index'])
	errors_shape = (num_rates, num_Kks, num_intensities)
	total_errors = sp.zeros(errors_shape)
	
	print (total_errors.shape)
	
	while not it.finished:
	
		rate, Kk, intensity = it.multi_index
		
		# Are ORNs active (1) or not (0)? Shape = (timepoints, signal IDs, Mm)
		active_array = (Yy[:, rate, :, Kk, intensity, :] > active_rate)
		
		# Num. of active receptors; has shape (timepoints, signal IDs)
		num_activated = sp.sum(active_array, axis=-1)
		errors = (errors_nonzero[:, rate, :, Kk, intensity]*\
					errors_zero[:, rate, :, Kk, intensity]*1./100)
		
		# list of Mm lists; each holds errors when iM receptors are active
		error_lists = [[] for x in sp.arange(Mm)]
				
		# For each signal, append the errors to the Mm lists in error_lists
		for idSs in range(num_odors):
			active_error_arr = sp.vstack((num_activated[:, idSs].T, 
											errors[:, idSs].T)).T
			
			# Find errors at all times at which a given # (iM) of ORNs active
			for iM in range(Mm):
				errors_iM = errors[sp.where(active_error_arr == iM)[0], idSs]
				error_lists[iM] = error_lists[iM] + errors_iM.tolist()
		
		# Get minimum # of receptors to get to min_accurate_pct accuracy
		errors_avg = sp.zeros(Mm)
		for iM in range(Mm):
			errors_avg[iM] = sp.mean(error_lists[iM])
		min_Mm = sp.where(errors_avg > min_accurate_pct)[0] + 1
		if len(min_Mm) == 0:
			min_Mm = [Mm]
		total_errors[it.multi_index] = min_Mm[0]
		
		it.iternext()
	
	# Plot heatmap
	fig = primacy_min_active_for_accurate_heatmap()
	X, Y = sp.meshgrid(range(num_Kks), range(num_intensities)[1:])
	plt.pcolormesh(Y, X, total_errors[0, :, 1:].T, cmap=plt.cm.hot)
	#plt.ylim(range(num_intensities)[1], range(num_intensities)[-1])
	save_fig('min_num_active_for_accurate_vs_Kk_intensity_heatmap', 
				subdir=data_flag)
	
	# Separate figure for colorbar
	fig = plt.figure()
	fig.set_size_inches(1, 5)
	ax1 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
	norm = mpl.colors.Normalize(vmin=0, vmax=Mm)
	ticks = range(0, Mm + 1, 10)
	cbar = mpl.colorbar.ColorbarBase(ax1, cmap=plt.cm.hot,
							norm=norm, ticks=ticks,
							orientation='vertical')
	cbar.ax.tick_params(labelsize=25)
	#cbar.ax.set_yticklabels(tick_labels)
	save_fig('min_num_active_colorbar', subdir=data_flag, 
				tight_layout=False)
	
if __name__ == '__main__':
	data_flag = get_flag()
	accuracy_vs_num_active(data_flag)