"""
Plot the spread of adapted background values.

Created by Nirag Kadakia at 14:33 11-14-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from four_state_receptor_CS import four_state_receptor_CS
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
#from encode_CS import single_encode_CS
from lin_alg_structs import random_matrix

data_flag = get_flag()
iter_var_idxs = list(map(int, sys.argv[2:]))

list_dict = read_specs_file(data_flag)
vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
obj = four_state_receptor_CS(**vars_to_pass)
	
# Set the temporal signal array from file; truncate to signal window
#obj.mu_dSs_2 = 0
#obj.sigma_dSs_2 = 0
#obj = single_encode_CS(obj, list_dict['run_specs'])

activity_stats = [obj.adapted_activity_mu, obj.adapted_activity_sigma]
adapted_activity = random_matrix([obj.Mm], params=activity_stats, 
									seed=obj.seed_adapted_activity)

fig = plt.figure()
fig.set_size_inches(6, 1)
plt.yticks([])
plt.xticks(fontsize=14)
sizes = sp.linspace(30, 60, obj.Mm)
plt.scatter(adapted_activity*obj.firing_max, [1]*obj.Mm, 
			color='purple', alpha=0.6, s=sizes)
subdir = data_flag
save_fig('adapted_A0', subdir=subdir)


		