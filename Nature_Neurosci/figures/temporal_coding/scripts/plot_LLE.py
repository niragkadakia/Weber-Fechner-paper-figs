"""
Plot a locally linear embedding of the activity vector 
in time.

Created by Nirag Kadakia at 17:00 02-08-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
from sklearn.manifold import LocallyLinearEmbedding
import sys
sys.path.append('../src')
import matplotlib.pyplot as plt
from load_specs import read_specs_file
from save_load_data import load_temporal_errors, save_est_signal_zeros_fig, \
							save_est_signal_nonzeros_fig
from figure_plot_formats import est_signal_zeros_subfigures, \
								est_signal_zeros_subfigures


def plot_LLE(data_flag, iter_vars_idx_to_plot=0, 
						avg_var_idx_to_plot=0, dts_to_plot=[1.0, 1.5], 
						zero_ylims=[-1, 1], nonzero_ylims=[-1, 1]):
	"""
	Plot estimated signal zero and nonzero components at two distinct 
	times in a temporal coding trace.
	
	Args:
		data_flag: string; data identifier.
		iter_vars_idx_to_plot: int; index of the iterated variable to be
			plotted. The iterated variable must be the first variable 
			of the iter_vars in specs.
		avg_var_idx_to_plot: int; index of the average variable  to be 
			plotted in successive plots. This variable must be the 
			second one of the iter_vars in specs.
		dts_to_plot: 2-entry list; times in the signal trace for which to plot.
			The nearest index of the time vector will be matched to this time.
		zero_ylims: 2-entry list; lower and upper limit of the signal 
			estimation plot for zero components. First value should be 
			negative, and second positive.
		nonzero_ylims: 2-entry list; lower and upper limit of the signal 
			estimation plot for nonzero components. First value should be 
			negative, and second positive.
		
	"""
	
	assert len(dts_to_plot) == 2, "dts_to_plot must be length-2 list"
	
	# Load data and get iterated variables and their dimensions
	data = load_temporal_errors(data_flag)
	list_dict = read_specs_file(data_flag)
	
	# Get min and max time indices
	Tt = data['Tt'] - data['Tt'][0]
	iT_lo = (sp.absolute(Tt - dts_to_plot[0])).argmin()
	iT_hi = (sp.absolute(Tt - dts_to_plot[1])).argmin()
	
	act = data['Yy']
	
	dYy_mean = sp.average(act[iT_lo:iT_hi, 9, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	dYy_std = sp.std(act[iT_lo:iT_hi, 9, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	plt.fill_between(Tt[iT_lo:iT_hi], dYy_mean - dYy_std, dYy_mean + dYy_std, color='r', alpha=0.5)
	plt.plot(Tt[iT_lo:iT_hi], dYy_mean, color='r')
	
	dYy_mean = sp.average(act[iT_lo:iT_hi, 10, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	dYy_std = sp.std(act[iT_lo:iT_hi, 10, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	plt.fill_between(Tt[iT_lo:iT_hi], dYy_mean - dYy_std, dYy_mean + dYy_std, color='b', alpha=0.5)
	plt.plot(Tt[iT_lo:iT_hi], dYy_mean, color='b')
	
	dYy_mean = sp.average(act[iT_lo:iT_hi, 11, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	dYy_std = sp.std(act[iT_lo:iT_hi, 11, :, :].reshape(iT_hi-iT_lo, -1), axis=1)
	plt.fill_between(Tt[iT_lo:iT_hi], dYy_mean - dYy_std, dYy_mean + dYy_std, color='g', alpha=0.5)
	plt.plot(Tt[iT_lo:iT_hi], dYy_mean, color='g')
	
	plt.show()
	
	""""
	for iT in sp.arange(iT_lo, iT_hi):
		LLE = LocallyLinearEmbedding(n_neighbors=2, n_components=2, reg = 0.01)
		Yy_1 = data['dYy'][iT:iT+4, iter_vars_idx_to_plot, 0, :]
		X_lle = LLE.fit_transform(Yy_1)
		plt.scatter(X_lle[:, 0], X_lle[:, 1], color='b')
	
	for iT in sp.arange(iT_lo, iT_hi):
		LLE = LocallyLinearEmbedding(n_neighbors=2, n_components=2, reg = 0.01)
		Yy_1 = data['dYy'][iT:iT+4, iter_vars_idx_to_plot, 8, :]
		X_lle = LLE.fit_transform(Yy_1)
		plt.scatter(X_lle[:, 0], X_lle[:, 1], color='r')
	plt.show()
	"""
	"""
	for sig in sp.arange(0, 10, 3):
		Yy = data['dYy'][iT_lo:iT_hi, iter_vars_idx_to_plot, sig, 0]
		plt.plot(Yy, color=plt.cm.Reds(0.2+1.*sig/15.))
		Yy = data['dYy'][iT_lo:iT_hi, iter_vars_idx_to_plot, sig, 8]
		plt.plot(Yy, color=plt.cm.Blues(0.2+1.*sig/15.))
	plt.show()
	"""
		
if __name__ == '__main__':
	data_flag = sys.argv[1]
	plot_LLE(data_flag, iter_vars_idx_to_plot=7, 
							avg_var_idx_to_plot=7, dts_to_plot=[18, 25],
							zero_ylims=[-0.01, 0.01], nonzero_ylims=[0, 0.12])