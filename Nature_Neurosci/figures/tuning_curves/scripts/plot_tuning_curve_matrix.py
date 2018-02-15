"""
Plot tuning curve matrix for all signal concentrations

spec files used to create data (data_flags):
figure_tuning_curves_adaptation_Kk2
figure_tuning_curves_no_adaptation_Kk2

Current figure (Figure_tuning_curves.svg) uses odor_seed 12 and 22

Created by Nirag Kadakia at 17:30 02-14-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
from scipy.ndimage.filters import gaussian_filter
import sys
sys.path.append('../src')
import matplotlib.pyplot as plt
params = {'text.usetex': False, 'mathtext.fontset': 'dejavusans'}
plt.rcParams.update(params)
from matplotlib import cm
from utils import get_flag
from load_specs import read_specs_file, parse_iterated_vars, \
						parse_relative_vars
from utils import get_flags, merge_two_dicts
from save_load_data import load_tuning_curve, save_tuning_curve_matrix_fig
from figure_plot_formats import tuning_curve_matrix_subfigure


def plot_tuning_curve_matrix(data_flag):
	
	# Which diversity of the system to plot?
	sigma_Kk2_idx = 2
	
	# Which signal indices to plot?
	mu_Ss0_idxs = [3, 20]
	
	# Receptors to plot
	receptors_to_plot = sp.arange(40)
	
	tuning_curve_data = load_tuning_curve(data_flag)
	tuning_curve = tuning_curve_data['tuning_curve']
	
	"""
	cmaps_avgs = [plt.cm.Blues_r, plt.cm.Oranges_r, plt.cm.Greens_r]
	for idx, receptor_idx in enumerate(receptors_to_plot):
		#means = sp.average(tuning_curve[:, sigma_Kk2_idx, :, receptor_idx], 
		#					axis=1)
		#stds = sp.std(tuning_curve[5:15, sigma_Kk2_idx, :, receptor_idx], axis=1)
		#plt.plot(range(20), means, color='0')
		#color = 1.*idx/50.
		#plt.fill_between(range(20), means-stds, means+stds, 
		#					facecolor='%s' % color)
		plt.plot(range(10), stds)
	plt.show()
	quit()
	"""
	
	list_dict = read_specs_file(data_flag)
	for key in list_dict:
		exec("%s = list_dict[key]" % key)
	
	for idx, receptor_idx in enumerate(receptors_to_plot):
		fig = tuning_curve_matrix_subfigure()
		data = tuning_curve[mu_Ss0_idxs[0]:mu_Ss0_idxs[1], 
				sigma_Kk2_idx, :, receptor_idx]
		plt.pcolormesh(data.T, cmap=plt.cm.hot, rasterized=True)
		plt.colorbar()
		save_tuning_curve_matrix_fig(fig, sigma_Kk2_idx, 
										receptor_idx, data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_tuning_curve_matrix(data_flag)
	