"""
Plot the percent of correctly decoded on a heatmap, as a function
of active ORNS (x) and odor complexity Kk (y)

Created by Nirag Kadakia at 09:55 05-17-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from plot_formats import fig_primacy_errors_time
from save_load_figure_data import load_binary_errors, load_activities, \
									save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_signal_trace_from_file


def active_heatmap(data_flag, decoded_pct=75, Kk_idxs_to_plot=[0, 2, 4, 6, 8],
					intensity_idxs_to_plot=range(1, 9)):
	
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
	
	signal_file = list_dict['fixed_vars']['signal_trace_file']
	signal_data = load_signal_trace_from_file(signal_file)
	signal_array = list_dict['iter_vars']['signal_trace_multiplier']
	offset = list_dict['fixed_vars']['signal_trace_offset']
	signal_intensities = max(offset + signal_data[:, 1])*\
							signal_array[intensity_idxs_to_plot]
	Tt = signal_data[:, 0] - signal_data[0, 0]
	
		
	# Separate figure for each rate
	for iRate in range(num_rates):
		
		# If decoding accuracy not reached, throw that signal out
		nans = 0
		
		fig = fig_primacy_errors_time()
	
		for iKk, Kk in enumerate(Kk_idxs_to_plot):
			
			decoded_times = sp.zeros((len(intensity_idxs_to_plot), num_odors))
			latencies = sp.zeros(decoded_times.shape)
			
			# Errors for each intensity
			for iI, intensity_idx in enumerate(intensity_idxs_to_plot):
				errors = 1./100*errors_nonzero[:, iRate, :, iKk, intensity_idx]*\
							errors_zero[:, iRate, :, iKk, intensity_idx]
				for idSs in range(num_odors):
					active_times = sp.where(errors[:, idSs] >= decoded_pct)[0]
					if len(active_times) > 0:
						decoded_times[iI, idSs] = Tt[active_times[0]]
					else:
						nans += 1
						decoded_times[iI, idSs] = sp.nan
			
			# Latencies for each intensity, for each odor, then avgd over odors
			for iI in range(len(intensity_idxs_to_plot)):
				latencies[iI, :] = decoded_times[iI, :] - decoded_times[0, :]
			avg_latencies = sp.nanmean(latencies, axis=1)
			
			color = plt.cm.viridis(0.85*iKk/(len(Kk_idxs_to_plot) - 1))
			plt.plot(signal_intensities, 1e3*avg_latencies, color=color, lw=2)

		print 'Thrown out indices = %s out of %s' % (nans, (iKk + 1)*num_odors*(iI + 1))
		save_fig('primacy_errors_vs_time_rate=%s_intensities=%s' 
					% (iRate, intensity_idxs_to_plot), subdir=data_flag)
		
		
if __name__ == '__main__':
	data_flag = get_flag()
	active_heatmap(data_flag)