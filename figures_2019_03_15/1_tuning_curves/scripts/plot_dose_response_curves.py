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

colors = plt.cm.viridis(sp.linspace(0.2, 0.5, 2))

# Background signal and sweep the signal for dose response
mu_Ss0s = iter_vars['mu_Ss0']
dSss = 10**sp.linspace(-1, 6, 100)

# How many individaul receptors to plot
num_receps = 20

# Doesn't matter what background is for unadapted
iter_var_idxs = [0]
vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)

# Unadapted tuning curve; just plot once
obj = four_state_receptor_CS(**vars_to_pass)
list_dict['run_specs']['run_type'] = ['power_Kk']
obj = single_encode_CS(obj, list_dict['run_specs'])
vals = sp.zeros((len(dSss), obj.Mm))
for iS, dSs in enumerate(dSss):
	obj.Ss = sp.zeros(obj.Nn)
	obj.Ss[obj.idxs] = dSs
	obj.set_measured_activity()
	vals[iS, :] = obj.Yy

fig = plt.figure(figsize = (2.8, 2.))
ax = plt.subplot(111)
for iM in range(obj.Mm)[:num_receps]:
	plt.plot(dSss, vals[:, iM], color='k', lw=0.75, alpha=0.8)
plt.xscale('log')
plt.tick_params(axis='both', labelsize=12)
plt.yticks([0, 0.5, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
save_fig('dose_response_unadapted', subdir=data_flag)
	
# Plot diff concentrations; same background
fig = plt.figure(figsize = (2.8, 2.))
ax = plt.subplot(111)
for iSs0, mu_Ss0 in enumerate(mu_Ss0s):

	iter_var_idxs = [iSs0]
	vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
	obj = four_state_receptor_CS(**vars_to_pass)
	list_dict['run_specs']['run_type'] = ['power_Kk_adapted']
	obj = single_encode_CS(obj, list_dict['run_specs'])
	
	vals = sp.zeros((len(dSss), obj.Mm))
	for iS, dSs in enumerate(dSss):
		obj.Ss[obj.idxs] = dSs
		obj.set_measured_activity()
		vals[iS, :] = obj.Yy
	for iM in range(obj.Mm):
		plt.plot(dSss, vals[:, iM], color=colors[iSs0], lw=0.75, alpha=0.8)
	plt.plot(dSss, sp.median(vals, axis=1), color=colors[iSs0], lw=2)
plt.xscale('log')
plt.tick_params(axis='both', labelsize=12)
plt.yticks([0, 0.5, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
save_fig('dose_response_adapted', subdir=data_flag)


# Plot diff concentrations; same background
fig = plt.figure(figsize = (2.8, 2.))
ax = plt.subplot(111)
for iSs0, mu_Ss0 in enumerate(mu_Ss0s):

	print (mu_Ss0)
	
	iter_var_idxs = [iSs0]
	vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
	obj = four_state_receptor_CS(**vars_to_pass)
	list_dict['run_specs']['run_type'] = ['power_Kk_adapted']
	obj = single_encode_CS(obj, list_dict['run_specs'])
	
	idxs_bck = obj.idxs
	
	perfect = sp.zeros(len(dSss))
	vals = sp.zeros((len(dSss), obj.Mm))
	for iS, dSs in enumerate(dSss):
	
		# Assuming same odor
		obj.Ss = sp.zeros(obj.Nn)
		obj.Ss[idxs_bck] = dSs
		obj.set_measured_activity()
		perfect[iS] = sp.mean(obj.Yy)
		
		seed = 0
		sp.random.seed(seed)
		obj.idxs = sp.random.randint(0, obj.Nn, obj.Kk)
		obj.Ss = sp.zeros(obj.Nn)
		obj.Ss[obj.idxs] = dSs
		obj.set_measured_activity()
		vals[iS, :] = obj.Yy
	for iM in range(obj.Mm)[:num_receps]:
		plt.plot(dSss, vals[:, iM], color=colors[iSs0], lw=0.75, alpha=0.8)
	plt.plot(dSss, sp.median(vals, axis=1), color='k', lw=2.5, ls='--')
	
plt.xscale('log')
plt.tick_params(axis='both', labelsize=12)
plt.yticks([0, 0.5, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
save_fig('dose_response_adapted_diff_bck', subdir=data_flag)
