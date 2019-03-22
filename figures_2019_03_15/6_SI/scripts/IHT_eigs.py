"""
Run a CS decoding run for one given index of a set of iterated
variables, iterating over many different signals. Basically a 
testing script.

Created by Nirag Kadakia at 13:00 09-01-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import scipy.linalg as LA
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from kinetics import linear_gain, receptor_activity
from encode_CS import single_encode_CS
	
	
# Number of signals to process
nSig = 1000
	
data_flags = get_flags()
list_dict = read_specs_file(data_flags[0])
assert len(list_dict['iter_vars']) == 2, "Only 2 iter var, mu_Ss0 and Kk"
mu_Ss0_vals = list_dict['iter_vars']['mu_Ss0']	
Kk_vals = list_dict['iter_vars']['Kk']	

colors = ['blue', 'orange']
lws = [2, 3]

for iK, Kk in enumerate(Kk_vals):
	for iMu, mu_Ss0 in enumerate(mu_Ss0_vals):
		
		fig = plt.figure()
		fig.set_size_inches(4, 3)
		plt.xticks(fontsize=16)
		plt.yticks(sp.arange(0, 50), fontsize=16)
		plt.xlim(-0.2, 5)
		plt.axvline(2, color='k')
		
		
		for iFlag, data_flag in enumerate(data_flags):
			list_dict = read_specs_file(data_flag)
			vars_to_pass = compile_all_run_vars(list_dict, [iMu, iK])
			
			eigs = []
			for iS in range(nSig):
				obj = four_state_receptor_CS(**vars_to_pass)
				obj.seed_dSs = iS
				obj = single_encode_CS(obj, list_dict['run_specs'])
				
				# Only get activity (not firing rate), and linear response
				# If activity is thresholded, remove linear response
				A = obj.Yy#receptor_activity(obj.Ss, obj.Kk1, obj.Kk2, obj.eps)
				R = obj.Rr#linear_gain(obj.Ss, obj.Kk1, obj.Kk2, obj.eps)
				for iM in range(obj.Mm):
					R[A < obj.NL_threshold, :] = 0.
				
				# Normalize columns
				for iN in range(obj.Nn):
					norm = (sp.sum(R[:, iN]**2.0))**0.5 
					if norm > 0:
						R[:, iN] = R[:, iN]/norm
				
				# Calculate absolute values of eigenvalues of m x k submatrix
				sp.random.seed()
				rand_idxs = sp.random.choice(obj.Nn, obj.Kk, replace=False)
				w, _ = LA.eig(sp.dot(R[:, rand_idxs].T, R[:, rand_idxs]))
				eigs.extend(abs(w))
			print (max(eigs))
			bins = sp.linspace(0, 10.0, 100)
			hist, _ = sp.histogram(eigs, bins=bins, normed=True)
			plt.plot((bins[1:] + bins[:-1])/2.0, hist, color=colors[iFlag], 
					 lw=lws[iFlag])
		
		for data_flag in data_flags:
			subdir='%s' % data_flag
			save_fig('Kk=%s_mu_Ss0=%s' % (Kk, mu_Ss0), subdir=subdir)

