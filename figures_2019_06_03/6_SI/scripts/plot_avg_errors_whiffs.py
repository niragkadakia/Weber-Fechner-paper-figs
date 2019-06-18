"""
Plot average errors over whiffs for different data flags.

Created by Nirag Kadakia at 10:40 05-06-2018
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
from save_load_figure_data import load_tf_success_ratios, save_fig_no_whtspc
from save_load_figure_data import load_success_ratios
from plot_formats import fig_avg_whiff_errors

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file
from load_data import load_signal_trace_from_file


def plot_avg_errors(data_flags, whf_thresh=4,
					int_window_rate_mults=[0.1, 1, 5], 
					bck_intensities=sp.array([1, 15, 40, 100]),
					plot_x_shift=1):
					
	"""
	Will only plot 2 rates (indices 0 and 1 for temporal_adaptation_rate
	index in data file). Will consider various forgetting times.
	"""

	# This should come from the data flags themselves
	# i.e. data_flags are ...
	#	jul18_power_law_temporal_bkgrnd_step_EA_2_mult=1_broad_tA
	#	jul18_power_law_temporal_bkgrnd_step_EA_2_mult=15_broad_tA
	#	jul18_power_law_temporal_bkgrnd_step_EA_2_mult=40_broad_tA
	#	jul18_power_law_temporal_bkgrnd_step_EA_2_mult=100_broad_tA
	
	
	for iMult, int_window_rate_mult in enumerate(int_window_rate_mults):
	
		for iFlag, data_flag in enumerate(data_flags):
			
			success = load_tf_success_ratios(data_flag, int_window_rate_mult)
			list_dict = read_specs_file(data_flag)
			iter_vars_dims = []
			for iter_var in list_dict['iter_vars']:
				iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
			iter_vars = list_dict['iter_vars']
			
			iter_var_names = ['temporal_adaptation_rate', 'seed_dSs', 'Kk_1', 'Kk_2']
			for iName, name in enumerate(iter_var_names):
				assert list(iter_vars.keys())[iName] == name, "%sth variable "\
					"must have name %s" % (iName, name)
			
			try: 
				avg_whf_errs
			except:
				err_arr_dims = (len(data_flags), len(int_window_rate_mults)) \
										+ (iter_vars_dims[0], iter_vars_dims[2], 
										iter_vars_dims[3])
				avg_whf_errs = sp.zeros(err_arr_dims)
				
			# Signal data can be loaded from specs file -- no need to open agg objs.
			signal_file = list_dict['fixed_vars']['signal_trace_file']
			signal_data = load_signal_trace_from_file(signal_file)
			Tt = sp.arange(len(signal_data))
			multiplier = list_dict['fixed_vars']['signal_trace_multiplier']
			offset = list_dict['fixed_vars']['signal_trace_offset']
			signal = (offset + signal_data[:, 1])*multiplier
			
			# Clip array to desired range
			xlims = (0, 1)
			xlim_idxs = [int(len(Tt)*xlims[0]), int(len(Tt)*xlims[1])]
			plot_range = range(xlim_idxs[0], xlim_idxs[1])
			Tt = Tt[plot_range]
			Tt = Tt - Tt[0]
			dt = Tt[1] - Tt[0]
			signal = signal[plot_range]
			success = success[plot_range,...]
			
			# Get the whiff occurences
			whf_bin = 1.*(signal > whf_thresh)
			whf_begs = Tt[sp.where(sp.diff(whf_bin) == 1)[0]]
			whf_ends = Tt[sp.where(sp.diff(whf_bin) == -1)[0]]
			if whf_bin[0] == 1:
				whf_begs = sp.hstack((sp.zeros(1), whf_begs))
			if len(whf_begs) > len(whf_ends):
				whf_ends = sp.hstack((whf_ends, Tt[-1]))
			if len(whf_begs) < len(whf_ends):
				whf_begs = sp.hstack((Tt[0], whf_begs))
			whf_begs = whf_begs.astype(int)
			whf_ends = whf_ends.astype(int)
			
			cum_success = 0
			num_whfs = 0
			
			# For each whiff -- see max decoding accuracy during whiff
			for nWhf in range(len(whf_begs)):
				whf_range = range(whf_begs[nWhf], whf_ends[nWhf])
				success[whf_range,...]
				cum_success += sp.amax(success[whf_range,...], axis=0)
				num_whfs += 1
			whf_success = 1.*cum_success/num_whfs
			avg_whf_succ = sp.average(whf_success, axis=1)*100.0
			avg_whf_errs[iFlag, iMult, ...] = avg_whf_succ
			print ("Number of whiffs = %s" % num_whfs)
			
		x_range = bck_intensities*plot_x_shift
	
	# Separate figure for each foreground/background complexity pair
	for Kk_1 in range(iter_vars_dims[2]):
		for Kk_2 in range(iter_vars_dims[3]):
			
			#fig = fig_avg_whiff_errors()
			fig = plt.figure()
			for iMult, int_window_rate_mult in enumerate(int_window_rate_mults):
				
				# Different values of color for each forgetting time
				color_val = 0.3	 + iMult*0.55/(len(int_window_rate_mults) - 1)
				
				# Fast adaptation in red, slow adaptation in green
				plt.plot(x_range, avg_whf_errs[:, iMult, 1, Kk_1, Kk_2],
							color=plt.cm.Reds(color_val), lw=5)
				plt.plot(x_range, avg_whf_errs[:, iMult, 0, Kk_1, Kk_2], 
							color=plt.cm.Blues(color_val), lw=5)
				plt.yticks([0, 50, 100])
				plt.xscale('log')
				plt.xlim(0.8, 115)
				plt.ylim(-5, 105)
			save_fig_no_whtspc('%s_Kk_1=%s_Kk_2=%s' % (bck_intensities, Kk_1, Kk_2),
						subdir='avg_whf_errs', no_ax=False)
		
if __name__ == '__main__':
	data_flags = get_flags()
	plot_avg_errors(data_flags)