"""
Plot illustrateive firing rate responses to step stimuli. 

Created by Nirag Kadakia at 13:00 09-20-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import os
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from plot_formats import fig_tuning_curves_traces
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir, scripts_dir
sys.path.append(src_dir())
sys.path.append(scripts_dir())
from temporal_CS_run import temporal_CS_run
from utils import get_flag
from load_specs import read_specs_file

Mm_to_plot = [1, 5, 6, 43, 46, 34, 35]
ord_idx = 1000

data_flag = get_flag()
iter_var_idxs = list(map(int, sys.argv[2:]))
list_dict = read_specs_file(data_flag)
for key in list_dict:
	exec("%s = list_dict[key]" % key)
obj_list = temporal_CS_run(data_flag, iter_var_idxs, mu_dSs_multiplier=1/3., 
						   signal_window=[0, 2000], sigma_dSs_multiplier=1/9., 
						   save_data=False, decode=False)

# Hold activity vectors for each iM in Mm_to_plot; order by value at ord_idx
fig = fig_tuning_curves_traces()
act_vecs = []
mags = []
for iM in Mm_to_plot:
	act = []
	for obj in obj_list:
		act.append(obj.Yy[iM])
	act_vecs.append(act)
	mags.append(act[ord_idx])
	
# Plot in order of magnitude of last value
sort_idxs = sp.argsort(mags)[::-1]
print (sort_idxs)
for idx, iM in enumerate(Mm_to_plot):
	act_to_plot = act_vecs[sort_idxs[idx]]
	col_val = 0.2 + 0.7*idx/(len(Mm_to_plot) - 1)
	color=plt.cm.inferno_r(col_val)
	lw = 2.5 #- 2*idx/(len(Mm_to_plot) - 1)
	plt.plot(obj.signal_trace_Tt, act_to_plot, lw=lw, color=color)
plt.show()
#fig_name = '%s_Mm=%s' % (iter_var_idxs, Mm_to_plot)
#save_fig(fig_name, subdir=data_flag, clear_plot=True, tight_layout=True)

