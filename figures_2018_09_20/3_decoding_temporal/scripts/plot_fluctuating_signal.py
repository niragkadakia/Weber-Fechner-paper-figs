"""
Plot estimation error in temporal coding tasks in a false color plot.

Created by Nirag Kadakia at 22:00 04-16-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
from scipy.ndimage.filters import gaussian_filter
import matplotlib
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import load_binary_errors, load_success_ratios, \
									save_fig
from plot_formats import fig_signal_trace

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file
from load_data import load_aggregated_temporal_objects, \
						load_signal_trace_from_file


def plot_temporal_errors(data_flag, whiff_threshold=6, ylim=20, 
						 min_whf_dur=0.15):

	list_dict = read_specs_file(data_flag)
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))		
	iter_vars = list_dict['iter_vars']
	
	iter_var_names = ['temporal_adaptation_rate', 'seed_dSs']
	for iName, name in enumerate(iter_var_names):
		assert list(iter_vars.keys())[iName] == name, "%sth variable "\
			"must have name %s" % (iName, name)
	
	# Signal data can be loaded from specs file -- no need to open agg objs.
	sig_file = list_dict['fixed_vars']['signal_trace_file']
	sig_data = load_signal_trace_from_file(sig_file)
	Tt = sig_data[:, 0] - sig_data[0, 0]
	multiplier = list_dict['fixed_vars']['signal_trace_multiplier']
	offset = list_dict['fixed_vars']['signal_trace_offset']
	signal = (offset + sig_data[:, 1])*multiplier
		
	# Clip array to desired section
	xlims = (0, 1.0)
	xlim_idxs = [int(len(Tt)*xlims[0]), int(len(Tt)*xlims[1])]
	plot_range = range(xlim_idxs[0], xlim_idxs[1])
	Tt = Tt[plot_range]
	Tt = Tt - Tt[0]
	signal = signal[plot_range]
	
	# Signal with whiffs highlighted
	whf_bin = 1.*(signal > whiff_threshold)
	whf_beg = Tt[sp.where(sp.diff(whf_bin) == 1)[0]]
	whf_end = Tt[sp.where(sp.diff(whf_bin) == -1)[0]]
		
	if whf_bin[0] == 1:
		whf_beg = sp.hstack((sp.zeros(1), whf_beg))
	if len(whf_beg) > len(whf_end):
		whf_end = sp.hstack((whf_end, Tt[-1]))
	if len(whf_beg) < len(whf_end):
		whf_beg = sp.hstack((Tt[0], whf_beg))
	
	fig = fig_signal_trace()
	fig.set_size_inches(12, 2.3)
	plt.plot(Tt, signal, lw=2, color=plt.cm.Greys(0.9))
	
	# Bar above graph for each whiff
	for nWhf in range(len(whf_beg)):
		width = whf_end[nWhf] - whf_beg[nWhf]
		if width < min_whf_dur:
			width = min_whf_dur
		xy = (whf_beg[nWhf], ylim + 1)
		rect = matplotlib.patches.Rectangle(xy, width, ylim/10., 
				alpha=0.5, color='purple',  clip_on=False, lw=0)
		ax = plt.gca()
		ax.add_patch(rect)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.xticks(sp.arange(0, 100, 10), fontsize=18)
	plt.yticks(sp.arange(0, 300, 10), fontsize=18)
	plt.xlim(Tt[0], Tt[-1])
	plt.ylim(0, ylim)
	
	save_fig('signal_with_whiffs', subdir=data_flag)
	
	
if __name__ == '__main__':
	data_flag = get_flag()
	plot_temporal_errors(data_flag)