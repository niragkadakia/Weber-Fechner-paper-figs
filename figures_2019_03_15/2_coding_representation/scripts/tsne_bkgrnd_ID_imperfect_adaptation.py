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
from utils import get_flags
from load_specs import read_specs_file, compile_all_run_vars
from entropy import response_entropy

# Run with 2 data flags
#   jul18_tsne_rand_bkgrnd_many_intensities_imperfect_WL_universal
#   jul18_tsne_rand_bkgrnd_many_intensities_imperfect_WL

iterations_per_beta = 1
cmap=plt.cm.inferno

data_flags = get_flags()
assert len(data_flags) == 2, "Only 2 data flags for universal and diverse"
scores_all = [[], []]

for iF, data_flag in enumerate(data_flags):

	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	assert (len(iter_vars.keys()) == 2) and (list(iter_vars.keys())[0] 
		== 'mu_dSs') and (list(iter_vars.keys())[1] == 'adaptive_beta_scaling_max'), \
		"Only 2 iter vars; must be `mu_dSs` and `adaptive_beta_scaling_max`"
	mu_dSs_vals = iter_vars['mu_dSs']
	betas = iter_vars['adaptive_beta_scaling_max']
	num_dSs = len(mu_dSs_vals)
	num_signals = list_dict['params']['num_signals']

	for iB, beta in enumerate(betas):

		scores = []

		for iT in range(iterations_per_beta):

			# Struc for saved activities; order is odor intensity, odor ID, receptor ID
			Yys = sp.zeros((num_dSs, num_signals, list_dict['params']['Mm']))
			Yys_old = sp.zeros((num_dSs, num_signals, list_dict['params']['Mm']))

			for idSs, mu_dSs in enumerate(mu_dSs_vals):

				iter_var_idxs = [idSs, iB]	
				vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
				obj = response_entropy(**vars_to_pass)

				# For each iT run, new signal seeds for both fore and back
				if idSs == 0:
					obj.seed_dSs = sp.random.seed(iT)
					obj.seed_dSs_2 = sp.random.seed(iT + 1000)

				# Basic setting of Kk1, etc.; Ss, Yy, and eps will be overwritten below.
				obj.encode_power_Kk()

				# First, set the foreground to zero; get random background
				# Store the fixed vals and seeds
				tmp_mu_dSs = obj.mu_dSs
				tmp_sigma_dSs = obj.sigma_dSs
				tmp_seed_Ss0 = obj.seed_Ss0
				obj.mu_dSs = 1e-5
				obj.sigma_dSs = 0
				sp.random.seed()
				obj.seed_Ss0 = sp.random.randint(1e7)
				obj.set_signal_array()
				obj.set_adapted_free_energy()
				obj.set_mean_response_array()
				Yys_old[idSs, :, :] = obj.Yy.T
				obj.mu_dSs = tmp_mu_dSs
				obj.sigma_dSs = tmp_sigma_dSs
				obj.seed_Ss0 = tmp_seed_Ss0

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
						print ('`%s` run type not accepted for tsne' % val[0])
						quit()
				else:
					print ('No entropy calculation run type specified, '\
							'proceeding with unadapted entropy calculation')
					obj.encode_entropy_calc()

				# Calculate the receptor response and save to array
				obj.set_mean_response_array()
				Yys[idSs, :, :] = obj.Yy.T

			print ('scores for this beta:', scores)
			print ('aggregated scores:', scores_all)

			# Do the dimensionality reduction with TSNE
			TSNE_func = TSNE(random_state=0)
			Yys_aggregated = Yys.reshape((-1, obj.Mm))
			Yys_old_aggregated = Yys_old.reshape((-1, obj.Mm))
			reduced_idxs = (num_dSs, num_signals, 2)
			reduced_data = TSNE_func.fit_transform(Yys_aggregated).reshape(reduced_idxs)

			# Plot only for 1 of the iterations; don't need to plot all
			if iT == 0:

				# Color by identity; size by intensity
				color_range = sp.linspace(0.1, 0.9, num_signals)
				marker_size_range = sp.linspace(15, 100, num_dSs)

				fig = fig_tnse()
				ax = plt.subplot()
				plt.xticks([])
				ax.set_xticks([])
				plt.yticks([])
				ax.set_yticks([])
				ax.spines['right'].set_visible(False)
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.spines['left'].set_visible(False)
				for iOdor in range(Yys.shape[1]):
					plt.scatter(reduced_data[:, iOdor, 0], reduced_data[:, iOdor, 1], 
								color=cmap(color_range[iOdor]), s=marker_size_range, 
								alpha=0.8)
				save_fig('tsne_Kk_1=%s_Kk_2=%s_nOdors=%s_beta=%.2f' % 
							(obj.Kk_1, obj.Kk_2, num_signals, beta), subdir=data_flag)


				fig = plt.figure(figsize=(2, 2))
				ax = plt.subplot()
				v1 = ax.violinplot(Yys_old_aggregated.flatten(), showextrema=False)
				for b in v1['bodies']:
					m = sp.mean(b.get_paths()[0].vertices[:, 0])
					b.get_paths()[0].vertices[:, 0] = \
						sp.clip(b.get_paths()[0].vertices[:, 0], -sp.inf, m)
				v2 = ax.violinplot(Yys_aggregated.flatten(), showextrema=False)
				for b in v2['bodies']:
					m = sp.mean(b.get_paths()[0].vertices[:, 0])
					b.get_paths()[0].vertices[:, 0] = \
						sp.clip(b.get_paths()[0].vertices[:, 0], m, sp.inf)
				plt.ylim(0, 300)
				plt.xticks([])
				plt.yticks([0, 150, 300])
				ax.spines['right'].set_visible(False)
				ax.spines['top'].set_visible(False)
				save_fig('tsne_Yy_both_Kk_1=%s_Kk_2=%s_nOdors=%s_beta=%.2f' % 
							(obj.Kk_1, obj.Kk_2, num_signals, beta), subdir=data_flag)


			# Calculate silhouette score
			vals = sp.reshape(reduced_data, (num_dSs*num_signals, 2), order='f')
			labels = []
			for iD in range(num_signals):
				labels.extend([iD]*num_dSs)
			scores.append(silhouette_score(vals, labels=labels))

		scores_all[iF].append(sp.mean(scores))

# Plot both universal and diverse
fig = fig_tnse()
fig.set_size_inches(6, 4)
ax = plt.subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.yticks(sp.arange(-0.2, 1, 0.1), ['-0.2', '', '0', '', '0.2', '', '0.4', 
		   '', '0.6', '', '0.8', '', '1.0' ])
plt.xlim(0, betas[-1]*1.05)
plt.ylim(-0.15, 0.7)
plt.plot(betas, scores_all[0], color='k')
plt.plot(betas/2.0, scores_all[1], color='k', ls='--')
save_fig('tsne_Kk_1=%s_Kk_2=%s_nOdors=%s' % (obj.Kk_1, obj.Kk_2, 
         num_signals), subdir=data_flag)
