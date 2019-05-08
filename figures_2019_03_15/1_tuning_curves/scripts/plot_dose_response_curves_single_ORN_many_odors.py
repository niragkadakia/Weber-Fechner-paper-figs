"""
Plot the dose responses of model neurons before and after adaptation


Created by Nirag Kadakia at 11:00 04-02-2019
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
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from encode_CS import single_encode_CS

data_flag = get_flag()
list_dict = read_specs_file(data_flag)
iter_vars = list_dict['iter_vars']

colors = plt.cm.inferno_r(sp.linspace(0.3, 0.5, 2))

# Background signal and sweep the signal for dose response
mu_Ss0s = iter_vars['mu_Ss0']
dSss = 10**sp.linspace(-1, 6, 100)

# How many individaul receptors to plot
iM_to_plot = 3
num_back_odors = 10
num_single_odors = 1
lw = 1.5
figsize = (3, 2)

# Doesn't matter what background is for unadapted
iter_var_idxs = [0]
vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)

# Unadapted tuning curve; just plot once
fig = plt.figure(figsize=figsize)
ax = plt.subplot(111)
obj = four_state_receptor_CS(**vars_to_pass)
list_dict['run_specs']['run_type'] = ['power_Kk']
obj = single_encode_CS(obj, list_dict['run_specs'])
vals = sp.zeros(len(dSss))
for iOdor in range(num_single_odors):
	for iS, dSs in enumerate(dSss):
		sp.random.seed(iOdor)
		obj.idxs = sp.random.randint(0, obj.Nn, obj.Kk)
		obj.Ss = sp.zeros(obj.Nn)
		obj.Ss[obj.idxs] = dSs
		obj.set_measured_activity()
		vals[iS] = obj.Yy[iM_to_plot]
	plt.plot(dSss, vals, color='k', lw=2, alpha=0.8)

# On same plot, Plot diff concentrations; same background
for iSs0, mu_Ss0 in enumerate(mu_Ss0s):
	for iOdor in range(num_single_odors):
		iter_var_idxs = [iSs0]
		vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
		obj = four_state_receptor_CS(**vars_to_pass)
		list_dict['run_specs']['run_type'] = ['power_Kk_adapted']
		obj.seed_Ss0 = iOdor
		obj.seed_dSs = iOdor
		obj = single_encode_CS(obj, list_dict['run_specs'])

		vals = sp.zeros(len(dSss))
		for iS, dSs in enumerate(dSss):
			obj.Ss[obj.idxs] = dSs
			obj.set_measured_activity()
			vals[iS] = obj.Yy[iM_to_plot]
		plt.plot(dSss, vals, color=colors[iSs0], lw=2, alpha=0.8)
		
plt.xscale('log')
plt.tick_params(axis='both', labelsize=10)
plt.yticks([0, 0.5, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
save_fig('dose_response_single_odor', subdir=data_flag)


# Plot diff concentrations; same background
fig = plt.figure(figsize = (2.8, 2.))
ax = plt.subplot(111)
for iSs0, mu_Ss0 in enumerate(mu_Ss0s):
	
	same_odor_vals = sp.zeros(len(dSss))
	for iOdor in range(num_back_odors):
		seed = iOdor

		iter_var_idxs = [iSs0]
		vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
		obj = four_state_receptor_CS(**vars_to_pass)
		obj.seed_Ss0 = iOdor
		obj.seed_dSs = iOdor
		list_dict['run_specs']['run_type'] = ['power_Kk_adapted']
		obj = single_encode_CS(obj, list_dict['run_specs'])

		vals = sp.zeros(len(dSss))
		
		idxs_bck = obj.idxs
		
		for iS, dSs in enumerate(dSss):
			
			if iOdor == 0:
				obj.Ss[idxs_bck] = dSs
				obj.set_measured_activity()
				same_odor_vals[iS] = obj.Yy[iM_to_plot]
		

			seed = 0
			sp.random.seed(seed)
			obj.idxs = sp.random.randint(0, obj.Nn, obj.Kk)
			obj.Ss = sp.zeros(obj.Nn)
			obj.Ss[obj.idxs] = dSs
			obj.set_measured_activity()
			vals[iS] = obj.Yy[iM_to_plot]
		plt.plot(dSss, vals, color=colors[iSs0], lw=lw, alpha=0.8)
	plt.plot(dSss, same_odor_vals, color='k', lw=2, alpha=0.8, ls='--')

plt.xscale('log')
plt.tick_params(axis='both', labelsize=10)
plt.yticks([0, 0.5, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
save_fig('dose_response_diff_bck', subdir=data_flag)