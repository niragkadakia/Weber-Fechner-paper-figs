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
from plot_formats import fig_primacy_errors_time
from save_load_figure_data import load_binary_errors, load_activities, \
									save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file


def primacy(data_flags, 
			decoded_pct=75, 
			primacy_sizes=sp.arange(1, 51, 1),
			rate_idx=0):
	
	assert len(data_flags) == 2, "Need exactly 2 data flags"
	
	pct_overlap = sp.zeros(len(primacy_sizes))
	
	for nP, num_prim in enumerate(primacy_sizes):

		prim_ORNs = None
		for iFlag, data_flag in enumerate(data_flags):
		
			list_dict = read_specs_file(data_flag)
			iter_vars_dims = []
			for iter_var in list_dict['iter_vars']:
				iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
			iter_vars = list_dict['iter_vars']
			iter_var_names = ['temporal_adaptation_rate', 'seed_dSs']
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
			
			# New data structure to hold the primacy ORNs at first step
			if prim_ORNs is None:
				prim_ORNs = sp.empty((num_odors, num_prim, len(data_flags)))
			
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
					prim_ORNs[idSs, :, iFlag] = recep_order[:num_prim]
				else:
					prim_ORNs[idSs, :, iFlag] = sp.nan
		
		# For each primacy size, pct of receptor overlaps for each signal
		good_signals = 0
		for idSs in range(num_odors):
			ORNs_lo = prim_ORNs[idSs, :, 0]
			ORNs_hi = prim_ORNs[idSs, :, 1]
			nans = sp.sum(sp.isnan(ORNs_lo*ORNs_hi))
			if nans > 0:
				continue
			pct_overlap[nP] += len(sp.intersect1d(ORNs_lo, ORNs_hi))/num_prim
			good_signals += 1
		pct_overlap[nP] = 1.*pct_overlap[nP]/good_signals
		
	plt.plot(primacy_sizes, pct_overlap, color='k', marker='o')
	plt.ylim(0.5, 1.0)
	plt.show()
		
			
if __name__ == '__main__':
	data_flags = get_flags()
	primacy(data_flags)