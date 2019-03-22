"""
Run a dimensionality reduction on the ORN respose to static values of 
odors to attempt to cluster. Use diff background identities but
single concentration for both foreground and background.
This is only for illustration purposes

Created by Nirag Kadakia at 16:00 03-21-2019
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
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
	# tsne_ex
	
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
	TSNE_func = TSNE(learning_rate=200, random_state=0)
	Yys_aggregated = Yys.reshape((-1, obj.Mm))
	reduced_idxs = (num_intensities, num_signals, 2)
	reduced_data = TSNE_func.fit_transform(Yys_aggregated).reshape(reduced_idxs)
	
	# Just plot for first 2 of the background identities
	for iOdor in range(2):
		for iBck in range(Yys.shape[0]):
			dat = sp.asarray([Yys[iBck, iOdor, :]]*2).T
			plt.imshow(dat, aspect=0.4, vmin=0, vmax=250, 
						cmap=plt.cm.afmhot, interpolation='nearest')
			plt.xticks([])
			plt.yticks([])
			save_fig('tsne_ex_iOdor=%s_iBck=%s' % (iOdor, iBck), subdir=data_flag)
	
	# Colorbar as separate figure
	fig = plt.figure()
	fig.set_size_inches(1, 5)
	ax = fig.add_axes([0.1, 0.1, 0.3, 0.7])
	ticks = [0, 50, 100, 150, 200, 250]
	norm = mpl.colors.Normalize(vmin=0, vmax=250)
	cbar = mpl.colorbar.ColorbarBase(ax, cmap=plt.cm.afmhot,
							norm=norm, ticks=ticks,
							orientation='vertical')
	cbar.ax.tick_params(labelsize=35)
	save_fig('tsne_ex_colorbar', subdir=data_flag, tight_layout=False)

	# Just project points from first 2 background identities to 2D space
	fig = fig_tnse()
	for iOdor in range(2):
		plt.scatter(reduced_data[:, iOdor, 0], reduced_data[:, iOdor, 1], 
					color='k', s=100, alpha=0.8)
	save_fig('tsne_ex', subdir=data_flag)
	

if __name__ == '__main__':
	data_flag = get_flag()
	tsne(data_flag)
