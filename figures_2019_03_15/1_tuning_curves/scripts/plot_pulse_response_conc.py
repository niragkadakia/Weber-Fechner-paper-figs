"""
Plot illustrateive firing rate responses to step stimuli
for a single ORN and single OR at different concentrations

Created by Nirag Kadakia at 13:00 04-25-2019
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import os
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from plot_formats import fig_tuning_curves_pulse, fig_tuning_curves_norm
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from encode_CS import single_encode_CS

dt = 1e-2
pulse_beg = 0.25
pulse_end = 0.75
pulse_rnd = 0.015
wind_len = 1.0
num_odors = 1
#pulse_heights = [0, 0.7, 4, 12, 40, 80, 160, 300, 600]
#pulse_heights = [0, 0.5, 1, 2, 3, 7, 15, 30, 60]
pulse_heights = [0, 0.7, 1.5, 4, 8, 15, 30, 60, 120]
ORN_to_plot = 15

data_flag = get_flag()
iter_var_idxs = list(map(int, sys.argv[2:]))
list_dict = read_specs_file(data_flag)

fig = fig_tuning_curves_pulse()
		
for iP, pulse_height in enumerate(pulse_heights):
	vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
	obj = four_state_receptor_CS(**vars_to_pass)
	obj.temporal_run = True

	# Set signal manually: step from  pulse_beg to pulse_end, rounded by pulse_rnd
	Tt = sp.arange(0, wind_len, dt)
	signal = pulse_height/2*(sp.tanh(((Tt - pulse_beg)/pulse_rnd)) + 1)*\
				(1 - 0.5*(sp.tanh(((Tt - pulse_end)/pulse_rnd)) + 1))
	signal += 1e-5
	obj.signal_trace_Tt = Tt
	obj.signal_trace = signal

	Yy_list = []
	for iT in range(len(Tt)):

		# Set estimation dSs values from signal trace -- mu_dSs must just 
		# be non-negative to get sparse_idxs
		obj.mu_Ss0 = obj.signal_trace[iT]
		obj.mu_dSs = obj.signal_trace[iT]*1e-5
		obj.sigma_dSs = 0

		# Encode / decode fully first time; then just update eps and responses
		if iT == 0:
			obj = single_encode_CS(obj, list_dict['run_specs'])

			# Spread adaptation rates over the system if sigma is set
			if obj.temporal_adaptation_rate_sigma != 0:
				obj.set_ordered_temporal_adaptation_rate()
		else:
			obj.set_sparse_signals()
			obj.set_temporal_adapted_epsilon()
			obj.set_measured_activity()
		Yy_list.append(obj.Yy[ORN_to_plot])
		
	col_val = 0 + 0.8*iP/(len(pulse_heights) - 1)
	lw = 2
	color=plt.cm.inferno(col_val)
	plt.plot(Tt, Yy_list, color=color, lw=lw, alpha=0.7)

plt.ylim(-5, 300)
fig_name = '%s_int=%s_iORN=%s' % (iter_var_idxs, pulse_height, ORN_to_plot)
save_fig(fig_name, subdir=data_flag, clear_plot=True)