"""
Find primacy sets for different background conditions

Created by Nirag Kadakia at 17:16 08-29-2018
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
from plot_formats import primacy_set_overlap_bckgrnd_ID, \
							primacy_set_overlap_bckgrnd_ID_consist
from save_load_figure_data import load_binary_errors, load_activities, \
									save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file


def primacy(data_flag, decoded_pct=75, primacy_sizes=sp.arange(1, 21, 1),
			rate_idx=0):
	"""
	Will calculate overlaps for all different primacy sets atop 
	different background odors
	"""
	
	prim_consist_pct = sp.zeros(len(primacy_sizes))
	
	for iP, num_prim in enumerate(primacy_sizes):

		list_dict = read_specs_file(data_flag)
		iter_vars_dims = []
		for iter_var in list_dict['iter_vars']:
			iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
		iter_vars = list_dict['iter_vars']
		iter_var_names = ['temporal_adaptation_rate', 'seed_dSs_2']
		for iName, name in enumerate(iter_var_names):
			assert list(iter_vars.keys())[iName] == name, "%sth variable "\
				"must have name %s" % (iName, name)
		num_rates = iter_vars_dims[0]
		num_odors = iter_vars_dims[1]
		Mm = list_dict['params']['Mm']
	
		# Combine nonzero and zero to get full decoded signal percentage
		binary_errors = load_binary_errors(data_flag)
		errors_nonzero = binary_errors['errors_nonzero']
		errors_zero = binary_errors['errors_zero']
		errors = 1.*binary_errors['errors_nonzero'][:, rate_idx, :]*\
						binary_errors['errors_zero'][:, rate_idx, :]/100
		Yy = load_activities(data_flag)
		
		# Hold all primacy sets for all odors for this iP
		concat_prim_sets = []
		
		# Total number of odors for which 75% accuracy actually reached.
		num_succ = 0
		
		# Get primacy ORNs for each odor signal
		for idSs in range(num_odors):
					
			# Get times at which the decoding percentage is above thresh
			acc_times = sp.where(errors[:, idSs] >= decoded_pct)[0]
			
			# First check if there is a point of 75% accuracy
			# Then find ordering of receptors that turned on before this time
			# Save orders in prim_ORNs, with inactives or inaccurates = nan.
			if len(acc_times) > 0:
				recep_on_times = []
				for iM in range(Mm):
					on_T = sp.where(Yy[:, rate_idx, idSs, iM] >= 15)[0]
					if len(on_T) > 0:
						if (on_T[0] <= acc_times[0]):
							recep_on_times.append(on_T[0])
						else:
							recep_on_times.append(1e5)
					else:
						recep_on_times.append(1e5)
				recep_order = sp.argsort(recep_on_times)
				concat_prim_sets.extend(recep_order[:num_prim])
				num_succ += 1
			
		concat_prim_sets = sp.array(concat_prim_sets)
		
		# Get percentage of time that each receptor is in the primacy set
		pct_in_primacy = []
		for iM in range(Mm):
			pct_in_primacy.append(sp.sum(1.*(concat_prim_sets == iM))/num_succ)
		
		# Plot each receptor, pct of odors for which in primacy set
		fig = primacy_set_overlap_bckgrnd_ID()
		
		# To calculate % consistency for primacy set in this iP
		pct = 0
		
		# Sort by most commonly in primacy set; color these purple
		sort_idxs = sp.argsort(pct_in_primacy)[::-1]
		sp.random.seed(iP)
		for iM in range(Mm):
			if iM in sort_idxs[:num_prim]:
				col = 'orange'
				shape = 'o'
				s = 80
				ec='purple'
				pct += pct_in_primacy[iM]
			else:
				col = 'orange'
				shape = 'o'
				s = 40
				ec = 'orange'
			plt.scatter(pct_in_primacy[iM] + sp.random.normal
						(0, (1 - pct_in_primacy[iM])*0.01), 
						sp.random.normal(0, 0.4), color=col, 
						#alpha=0.3 + 0.7*(Mm - iM)/Mm, 
						#s=10+iM*2, lw=1, marker=shape)
						alpha=0.7,
						s=s, lw=2, marker=shape,
						edgecolor=ec)
		save_fig('pct_of_time_in_primacy_set_p=%s' % num_prim, subdir=data_flag)
		
		prim_consist_pct[iP] = pct/num_prim
	
	fig = primacy_set_overlap_bckgrnd_ID_consist()
	plt.plot(primacy_sizes, prim_consist_pct, color='purple', lw=3)
	plt.xlim(0, num_prim + 1)
	save_fig('primacy_set_consistency', subdir=data_flag)
		
		
			
if __name__ == '__main__':
	data_flag = get_flag()
	primacy(data_flag)