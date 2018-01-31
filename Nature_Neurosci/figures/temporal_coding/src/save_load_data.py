"""
Functions for loading data for paper figures.

Created by Nirag Kadakia at 20:30 10-04-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import cPickle
import shelve
import gzip
import os
from local_methods import def_data_dir
from local_methods import def_analysis_dir
import matplotlib.pyplot as plt

DATA_DIR = def_data_dir()
ANALYSIS_DIR = def_analysis_dir()		


def load_objects(obj_idx, data_flag):
	"""
	Load objects saved by save_objects in the save_data module. 

	Args:
		obj_idx: List whose indices corresponding to name of saved *.npz 
					file holding object array.
		data_flag: Data identifier for loading and saving.
	"""

	filename = '%s/objects/%s/%s.pklz' % (DATA_DIR, data_flag, obj_idx)
	with gzip.open(filename, 'rb') as f:
		obj = cPickle.load(f)
	
	return obj

def save_temporal_errors(data, data_flag):
	"""
	Save full aggregated error and success data for temporal decoding.
	
	Args:
		data: dictionary, holding errors, success, epsilons for each 
			iterated variable index
		data_flag: Data identifier for loading and saving.
	"""
	
	out_dir = '%s/analysis/%s' % (DATA_DIR, data_flag)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	filename = '%s/temporal_coding_figures_errors.pklz' % out_dir
	with gzip.open(filename, 'wb') as f:
		cPickle.dump(data, f, protocol = 2)
	print ('\nTemporal error data file saved to %s' % filename)
	
def load_temporal_errors(data_flag):
	"""
	Save full aggregated error and success data for temporal decoding.
	
	Args:
		data: dictionary, holding errors, success, epsilons for each 
			iterated variable index
		data_flag: Data identifier for loading and saving.
	"""
	
	filename = '%s/analysis/%s/temporal_coding_figures_errors.pklz' \
				% (DATA_DIR, data_flag)
	with gzip.open(filename, 'rb') as f:
		data_dict = cPickle.load(f)
	
	return data_dict	
	
def save_signal_trace_fig(fig, data_flag, xlims, dual=False):

	"""
	Save signal trace subfigures
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	
	if dual == False:
		file_str = 'signal_trace_%s' % xlims
	else:
		file_str = 'signal_trace_2_%s' % xlims
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_epsilon_trace_fig(fig, data_flag, xlims):

	"""
	Save signal trace subfigures
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'epsilon_trace_%s' % xlims
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_nonzero_errors_trace_fig(fig, data_flag, xlims, dual=False):

	"""
	Save signal trace subfigures
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	
	if dual == False:
		file_str = 'nonzero_errors_trace_%s' % xlims
	else:
		file_str = 'nonzero_errors_2_trace_%s' % xlims
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_zero_errors_trace_fig(fig, data_flag, xlims, dual=False):

	"""
	Save signal trace subfigures
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	
	if dual == False:
		file_str = 'zero_errors_trace_%s' % xlims
	else:
		file_str = 'zero_errors_2_trace_%s' % xlims
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_est_signal_zeros_fig(fig, data_flag, dts_to_plot, 
								avg_var_idx_to_plot, 
								iter_vars_idx_to_plot):

	"""
	Save signal estimation subfigures for zero components
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'est_signal_zeros_%s_%s_%s' % (dts_to_plot, 
												avg_var_idx_to_plot, 
												iter_vars_idx_to_plot)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_est_signal_nonzeros_fig(fig, data_flag, dts_to_plot, 
								avg_var_idx_to_plot, 
								iter_vars_idx_to_plot):

	"""
	Save signal estimation subfigures for nonzero components
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'est_signal_nonzeros_%s_%s_%s' % (dts_to_plot, 
												avg_var_idx_to_plot, 
												iter_vars_idx_to_plot)
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_num_wrong_zeros_vs_adapt_rate_fig(fig, data_flag, whf_thresh):
	"""
	Save performance vs adaptation rate figures for zero components
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'num_wrong_zeros_vs_adapt_rate_thresh=%s' % whf_thresh 
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_nonzero_pct_err_vs_adapt_rate_fig(fig, data_flag, whf_thresh):
	"""
	Save performance vs adaptation rate figures for zero components
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'nonzero_pct_err_vs_adapt_rate_thresh=%s' % whf_thresh 
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_plot_adapt_rate_vs_act_fig(fig, data_flag):
	"""
	Save adaptation rate versus activity level for ordered adaptation rates
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'adapt_rate_vs_act'
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_adapt_rate_discr_num_wrong_zero_3d_fig(fig, data_flag, 
										whf_threshs_1, whf_threshs_2):
	"""
	Save 3d plot of wrong nonzero elements as a function of 
	adaptataion rate and Kk_split
	"""
	
	if whf_threshs_1[1] > 50:
		whf_threshs_1[1] = sp.inf
	if whf_threshs_2[1] > 50:
		whf_threshs_2[1] = sp.inf
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'num_wrong_discr_threshs=[%.3f,%.3f],'\
		'[%.3f,%.3f]' % (whf_threshs_1[0], whf_threshs_1[1], 
		whf_threshs_2[0], whf_threshs_2[1])
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
def save_adapt_rate_discr_pct_error_nonzero_3d_fig(fig, data_flag, 
										whf_threshs_1, whf_threshs_2):
	"""
	Save 3d plot of pct error in nonzero elements as a function of 
	adaptataion rate and Kk_split
	"""
	
	out_dir = '%s/figures/temporal_coding/%s' % (ANALYSIS_DIR, data_flag)
	file_str = 'nonzero_pct_discr_threshs=[%.3f,%.3f],'\
		'[%.3f,%.3f]' % (whf_threshs_1[0], whf_threshs_1[1], 
		whf_threshs_2[0], whf_threshs_2[1])
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
	out_dir = '../subfigures/%s' % data_flag
	if not os.path.exists(out_dir): 
		os.makedirs(out_dir)
	
	filename = '%s/%s.png' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	filename = '%s/%s.svg' % (out_dir, file_str)
	plt.savefig(filename, bbox_inches = 'tight')
	
