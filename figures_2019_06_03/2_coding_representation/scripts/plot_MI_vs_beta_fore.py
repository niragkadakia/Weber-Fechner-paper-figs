"""
Plot mutual information of olfactory systems as a function of
the degree of Weber Law adaptation. The MI is only calculated
for the foreground odor


Created by Nirag Kadakia at 12:40 06-03-2019
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
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects


def plot_MI_fore(data_flag, Kk_idxs=[4, 4], Tt_idxs=[0, 8, 15, 
				 40, 60, 80, 100]):
	"""
	Kk_idxs is a list of which of the iterated variables to plot
	(inner loop) for each data_flag (outer loop).
	"""
	
	# Should be the index right after the foreground odoor step
	Tt_idx = 104
	
	colors = plt.cm.viridis(sp.linspace(0.1, 0.7, len(Tt_idxs)))
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	

	assert len(iter_vars) == 3, "Need 3 iter_vars"
	iter_var_names = ['adaptive_beta_scaling_min', 'Kk_1', 'Kk_2']
	for iName, name in enumerate(iter_var_names):
		assert list(iter_vars.keys())[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	betas = iter_vars['adaptive_beta_scaling_min']
	Kk_1_vals = iter_vars['Kk_1']
	Kk_2_vals = iter_vars['Kk_2']

	print ('Loading object list...'),
	CS_object_array = load_aggregated_temporal_objects(data_flag)
	print ('...loaded.')
	
	# MI data vector has length [betas], (averaged over all receptors)
	data = sp.average(CS_object_array['entropy'][Tt_idx, :, Kk_idxs[0], 
					  Kk_idxs[1], :], axis=-1)
	
	fig = plt.figure(figsize=(6, 4))
	ax = plt.subplot()
	data = sp.average(CS_object_array['entropy']
					  [Tt_idx, :, Kk_idxs[0], Kk_idxs[1], :], axis=-1)
	
	plt.xticks(fontsize=18)
	plt.yticks(sp.arange(0, 6, 0.2), fontsize=18)
	plt.xlim(0, betas[-1]*1.05)
	plt.ylim(0.6, 2.51)
	plt.plot(betas, data, color='k', lw=2)
	
	# Only show every other tick label
	for label in ax.yaxis.get_ticklabels()[1::2]:
		label.set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	save_fig('MI_vs_beta_%s_%s' % (Kk_idxs, Tt_idx), subdir=data_flag)
	
				
if __name__ == '__main__':
	data_flag = get_flag()
	plot_MI_fore(data_flag)