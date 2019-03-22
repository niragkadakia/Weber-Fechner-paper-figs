"""
Plot the estimations returned by the IHT algorithm. Only written
for a specific iterated variable list in the specs file. 

Created by Nirag Kadakia at 12:00 11-06-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import load_IHT_est, save_fig
from plot_formats import fig_errors_vs_Kk

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from four_state_receptor_CS import four_state_receptor_CS
from load_specs import read_specs_file
from analysis import binary_errors, binary_success


# Set the bounds for binary estimation errors
nonzero_bounds=[0.7, 1.3]
zero_bound=1./10.

# Grab saved IHT estimation file
data_flag = get_flag()
list_dict = read_specs_file(data_flag)
iter_vars = list_dict['iter_vars']
assert len(iter_vars) == 3, "Need 3 iter_vars"
iter_var_names = ['mu_Ss0', 'seed_dSs', 'Kk']
for iName, name in enumerate(iter_var_names):
	assert list(iter_vars.keys())[iName] == name, "%sth variable "\
		"must have name %s" % (iName, name)

x_est_arr, x_true_arr = load_IHT_est(data_flag)
succ_arr = sp.zeros(x_est_arr.shape[:-1])

it = sp.nditer(succ_arr, flags=['multi_index'])
while not it.finished:
	
	# Need to call dummy CS object and population variables to get error
	obj = four_state_receptor_CS()
	obj.dSs = x_true_arr[it.multi_index]
	obj.Nn = len(obj.dSs)
	obj.dSs_est = x_est_arr[it.multi_index]
	obj.mu_dSs = iter_vars['mu_Ss0'][it.multi_index[0]]
	obj.idxs = sp.where(obj.dSs != 0)
	
	errors = binary_errors(obj, nonzero_bounds, zero_bound)
	zero_errs = errors['errors_zero']
	nonzero_errs = errors['errors_nonzero']
	succ_arr[it.multi_index] = binary_success(zero_errs,nonzero_errs)
	it.iternext()

avg_succ = 100.*sp.mean(succ_arr, axis=1)

x = iter_vars['mu_Ss0']
y = iter_vars['Kk']

# Pcolormesh quirk -- add extra row or it cuts it off
x = sp.log(x)/sp.log(10)
delx = x[1] - x[0]
dely = y[1] - y[0]
x = sp.hstack((x, [x[-1] + delx]))
y = sp.hstack((y, [y[-1] + dely]))
X, Y = sp.meshgrid(x, y)

vminmax = [0, 101]
fig = plt.figure()
fig.set_size_inches(5, 4)
plt.pcolormesh(X, Y, avg_succ.T, cmap=plt.cm.hot, 
				rasterized=True, vmin=vminmax[0], vmax=vminmax[1])
plt.xticks([])
plt.yticks([])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)

subdir='%s' % data_flag
save_fig('IHT_est', subdir=subdir)