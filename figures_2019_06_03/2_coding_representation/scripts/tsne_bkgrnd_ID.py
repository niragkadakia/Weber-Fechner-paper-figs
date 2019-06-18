"""
Run a dimensionality reduction on the ORN respose to static values of 
odors to attempt to cluster. Use diff background identities but
single concentration for both foreground and background.

Created by Nirag Kadakia at 17:26 19-22-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import matplotlib.pyplot as plt
import sys
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig
from plot_formats import fig_tnse

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
from entropy import response_entropy

def tsne(data_flag, cmap=plt.cm.inferno):

	# RUN ON 
	# jul18_tsne_rand_bkgrnd_WL
	# jul18_tsne_rand_bkgrnd_no_WL

	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	assert (len(iter_vars.keys()) == 1) and (list(iter_vars.keys())[0] 
		== 'mu_dSs'), "Only 1 iter var; must be `mu_dSs`"
	mu_dSs_vals = iter_vars['mu_dSs']
	num_intensities = len(mu_dSs_vals)
	num_signals = list_dict['params']['num_signals']
	
	# Struc for saved activities; order is odor intensity, odor ID, receptor ID
	Yys = sp.zeros((num_intensities, num_signals, list_dict['params']['Mm']))
		
	for idSs, mu_dSs in enumerate(mu_dSs_vals):
			
		iter_var_idxs = [idSs]	
		vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
		obj = response_entropy(**vars_to_pass)
		
		# Basic setting of Kk1, etc.; Ss, Yy, and eps will be overwritten below.
		obj.encode_power_Kk()
		
		# Set the signals and free energy, depending if adaptive or not.
		if 'run_type' in list_dict['run_specs'].keys():
			val = list_dict['run_specs']['run_type']
			if val[0] == 'entropy_calc':
				obj.encode_entropy_calc()
			elif val[0] == 'entropy_calc_adapted':
				obj.encode_entropy_calc_adapted()
			elif val[0] == 'entropy_calc_rand_bkgrnd':
				obj.encode_entropy_calc_rand_bkgrnd()
			elif val[0] == 'entropy_calc_adapted_rand_bkgrnd':
				obj.encode_entropy_calc_adapted_rand_bkgrnd()
			else:
				print ('`%s` run type not accepted for tsne calculation' % val[0])
				quit()
		else:
			print ('No entropy calculation run type specified, proceeding with' \
					'unadapted entropy calculation')
			obj.encode_entropy_calc()
		
		# Calculate the receptor response and save to array
		obj.set_mean_response_array()
		Yys[idSs, :, :] = obj.Yy.T
		
	# Do the dimensionality reduction with TSNE
	TSNE_func = TSNE(random_state=0)
	Yys_aggregated = Yys.reshape((-1, obj.Mm))
	reduced_idxs = (num_intensities, num_signals, 2)
	reduced_data = TSNE_func.fit_transform(Yys_aggregated).reshape(reduced_idxs)
	
	# Plot the distribution of the activities
	fig = plt.figure(figsize=(2, 2))
	ax = plt.subplot()
	ax.violinplot(Yys_aggregated.flatten(), showextrema=False)
	plt.ylim(0, 300)
	plt.yticks([0, 150, 300])
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	save_fig('tsne_Yy_Kk_1=%s_Kk_2=%s_nOdors=%s' % (obj.Kk_1, obj.Kk_2, 
			 num_signals), subdir=data_flag)
	
	# Project all points to 2D space; color by identity; size by intensity
	fig = fig_tnse()
	color_range = sp.linspace(0.1, 0.9, num_signals)
	marker_size_range = sp.linspace(15, 100, num_intensities)
	
	for iOdor in range(Yys.shape[1]):
		plt.scatter(reduced_data[:, iOdor, 0], reduced_data[:, iOdor, 1], 
					color=cmap(color_range[iOdor]), s=marker_size_range, 
					alpha=0.8)
	
	# Calculate silhouette score
	vals = sp.reshape(reduced_data, (num_intensities*num_signals, 2), order='f')
	labels = []
	for iD in range(num_signals):
		labels.extend([iD]*num_intensities)
	score = silhouette_score(vals, labels=labels)
	print ('Silhouette Score:', score)
	
	save_fig('tsne_Kk_1=%s_Kk_2=%s_nOdors=%s' % (obj.Kk_1, obj.Kk_2, 
				num_signals), subdir=data_flag)

	
if __name__ == '__main__':
	data_flag = get_flag()
	tsne(data_flag)
