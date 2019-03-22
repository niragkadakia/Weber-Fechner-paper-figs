"""
Plot illustrateive firing rate responses to step stimuli. 

Created by Nirag Kadakia at 13:00 09-20-2018
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

dt = 1e-3
pulse_beg = 0.25
pulse_end = 0.75
pulse_height = 300
pulse_rnd = 0.005
wind_len = 1.0
num_odors = 20

data_flag = get_flag()
iter_var_idxs = list(map(int, sys.argv[2:]))
list_dict = read_specs_file(data_flag)
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

# For each odor identity, generate time-traces of ORN responses
obj_arr = [[] for i in range(num_odors)]
Yy_arr = [[] for i in range(num_odors)]
for seed in range(num_odors):
		
	# Set new odor identity and run
	print (seed)
	obj.seed_dSs = seed
	
	obj_list = []
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
			
		Yy_list.append(obj.Yy)

	Yy_arr[seed] = Yy_list

# Plot, for each receptor, all normalized responses
for iM in range(obj.Mm):
	num_plots = 0
	fig = fig_tuning_curves_norm()
	for seed in range(num_odors):
		Yy_list = Yy_arr[seed]
		act = []
		for Yy in Yy_list:
			act.append(Yy[iM])
		
		# Only can normalize if response isn't too low
		if max(act) < 1: 
			continue
		col_val = 0.2 + 0.7*seed/(num_odors - 1)
		color=plt.cm.inferno(col_val)
		lw = 2
		plt.plot(Tt, act/max(act), color=color, lw=lw, alpha=0.7)
		num_plots += 1
	
	plt.ylim(-0.05, 1.1)
	fig_name = '%s_int=%s_iM=%s_norm' % (iter_var_idxs, pulse_height, iM)
	save_fig(fig_name, subdir=data_flag, clear_plot=True)

# Plot, for each receptor unnormalized responses
for iM in range(obj.Mm):
	fig = fig_tuning_curves_pulse()
	for seed in range(num_odors):
		Yy_list = Yy_arr[seed]
		act = []
		for Yy in Yy_list:
			act.append(Yy[iM])
		col_val = 0 + 0.8*seed/(num_odors - 1)
		lw = 2
		color=plt.cm.inferno(col_val)
		plt.plot(Tt, act, color=color, lw=lw, alpha=0.7)
	
	plt.ylim(-5, 305)
	fig_name = '%s_int=%s_iM=%s' % (iter_var_idxs, pulse_height, iM)
	save_fig(fig_name, subdir=data_flag, clear_plot=True)
