"""
Plot mutual information of olfactory systems as a function of the 
stimulus magnitude.


Created by Nirag Kadakia at 12:40 04-28-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig
from plot_formats import fig_MI_trace

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file
from load_data import load_aggregated_entropy_objects


def plot_MI(data_flags, Kk_idxs = [[3, 0], [3, 3]], dual_color=0):
	"""
	Kk_idxs is a list of which of the iterated variables to plot
	(inner loop) for each data_flag (outer loop).
	"""

	cmaps = [plt.cm.Blues, plt.cm.Greens]
	colors = ['b', 'g']
	#cmaps = [plt.cm.Reds, plt.cm.Blues]
	#colors = ['r', 'b']
	
	assert len(data_flags) <= 2, "Only plot 2 data flags at a time"
	
	fig = fig_MI_trace()
	for iFlag, data_flag in enumerate(data_flags):
	
		list_dict = read_specs_file(data_flag)
		iter_vars = list_dict['iter_vars']
		iter_vars_dims = []
		for iter_var in list_dict['iter_vars']:
			iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	
		
		assert len(iter_vars) == 3, "Need 3 iter_vars"
		iter_var_names = ['mu_dSs', 'Kk_1', 'Kk_2']
		for iName, name in enumerate(iter_var_names):
			assert iter_vars.keys()[iName] == name, "%sth variable "\
				"must have name %s" % (iName, name)
		mu_Ss0_vals = iter_vars['mu_dSs']
		Kk_1_vals = iter_vars['Kk_1']
		Kk_2_vals = iter_vars['Kk_2']
		x = mu_Ss0_vals
		
		print ('Loading object list...'),
		CS_object_array = load_aggregated_entropy_objects(data_flag)
		print ('...loaded.')
		data = CS_object_array['entropy'][:, Kk_idxs[iFlag][0], 
											Kk_idxs[iFlag][1], :]
		
		# Plot 2 traces with corresponding colors; dash/solid if dual
		if iFlag == 1:
			if dual_color == True:
				print 'h'
				plt.plot(x, sp.average(data, axis=-1), lw=2,
					color=colors[1], linestyle='--', 
					zorder=1000 + iFlag + 2)
			else: 
				plt.plot(x, sp.average(data, axis=-1), lw=2, color=colors[iFlag],
					zorder=1000 + iFlag)
				plt.plot(x, data, lw=0.5, color=cmaps[iFlag](0.5), alpha=0.5)
		else:
			plt.plot(x, data, lw=0.5, color=cmaps[iFlag](0.5), alpha=0.5)
			plt.plot(x, sp.average(data, axis=-1), lw=2, color=colors[iFlag],
					zorder=1000 + iFlag)
		
		plt.ylim(0,  8.01)
		plt.xlim(1, 10**4)
		
	save_fig('MI_trace_%s_%s' % (data_flags, Kk_idxs), subdir=data_flag, 
									tight_layout=False)
		
				
if __name__ == '__main__':
	data_flags = get_flags()
	plot_MI(data_flags)