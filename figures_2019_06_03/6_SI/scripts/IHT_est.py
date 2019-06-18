"""
Use the iterative hard thresholding algorithm to get at the 

Created by Nirag Kadakia at 15:00 11-05-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
sys.path.append('../../shared_src')
from save_load_figure_data import save_IHT_est

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from kinetics import receptor_activity, linear_gain
from encode_CS import single_encode_CS


# Maximum iterations for termination
max_its = 10000

# Threshold for no change in subsequent iterations
del_thresh = 1e-7

# Threshold for step size reduction if delta is stuck; step factor reduction
thresh_del_ch = 1e-8
step_reduce_factor = 0.9

# Threshold for how many sparse components vs Kk
sparsity_thresh_factor = 2

# Step factor / Step factor of mu_Ss0
step_factor = 1./20

data_flag = get_flag()
list_dict = read_specs_file(data_flag)
iter_vars_dims = []
for iter_var in list_dict['iter_vars']:
	iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
it = sp.nditer(sp.zeros(iter_vars_dims), flags=['multi_index'])	

x_est_arr = None
x_true_arr = None

while not it.finished:
	
	# Initialize object
	vars_to_pass = compile_all_run_vars(list_dict, it.multi_index)
	obj = four_state_receptor_CS(**vars_to_pass)
	obj = single_encode_CS(obj, list_dict['run_specs'])
	
	# Re-calculate activity (remove L-NL filter portion for now)
	Yy = receptor_activity(obj.Ss, obj.Kk1, obj.Kk2, obj.eps) 
	
	if x_est_arr is None:
		arr_shape = tuple(iter_vars_dims) + (obj.Nn, )
		x_est_arr = sp.zeros(arr_shape)
		x_true_arr = sp.zeros(arr_shape)

	# Initial guess is the zero vector
	x_est = sp.zeros((1, obj.Nn))
	thresh_k = obj.Kk*sparsity_thresh_factor
	step = obj.mu_Ss0*step_factor
	
	delta = 1
	iT = 0
	while delta > del_thresh:
	
		# Linear matrix evaluated at previous step
		A = linear_gain(x_est[iT], obj.Kk1, obj.Kk2, obj.eps)
		
		# Normalize columns to Euclidean norm 1
		for iN in range(obj.Nn):
			norm  = (sp.sum(A[:, iN]**2.0))**0.5
			A[:, iN] = A[:, iN]/norm
		
		est_Yy = receptor_activity(x_est[iT], obj.Kk1, obj.Kk2, obj.eps)
		x_unfilt = x_est[iT] + step*sp.dot(A.T, (Yy - est_Yy))
		
		# Sparsify
		sorted_idxs = sp.argsort(x_unfilt)[::-1][:thresh_k]
		x_filt = sp.zeros(len(x_unfilt))
		x_filt[sorted_idxs] = x_unfilt[sorted_idxs]
		x_est = sp.vstack((x_est, x_filt))
		
		# Check if step size needs to adapt
		delta_new = sp.sum((x_est[iT - 1] - x_est[iT])**2.0)/(obj.mu_Ss0)
		if abs(delta_new - delta) < thresh_del_ch:
			step *= step_reduce_factor
		delta = delta_new
		
		print (it.multi_index, delta)
		iT += 1
		if (iT > max_its):
			break
	x_est_arr[it.multi_index] = x_est[-1]
	x_true_arr[it.multi_index] = obj.Ss
	it.iternext()

# Save estimates and true vectors to file
save_IHT_est(x_est_arr, x_true_arr, data_flag)
