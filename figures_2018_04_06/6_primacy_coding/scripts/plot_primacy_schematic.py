"""
Subfigures to illustrate the primacy coding idea.

Created by Nirag Kadakia at 10:16 05-16-2018
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
from plot_formats import fig_primacy_schematic
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flag
from load_specs import read_specs_file, compile_all_run_vars
from four_state_receptor_CS import four_state_receptor_CS
from encode_CS import single_encode_CS


def plot_primacy_schematic(data_flag, seed_Kk2=16, vmin=10, vmax=30):
	
	list_dict = read_specs_file(data_flag)
	iter_vars = list_dict['iter_vars']
	Nn = list_dict['params']['Nn']
	Mm = list_dict['params']['Mm']
	iter_vars_dims = []
	for iter_var in list_dict['iter_vars']:
		iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	
	
	assert len(iter_vars_dims) == 2, "Need 2 	iterated vars; mu_Ss0 and Kk"
	assert iter_vars.keys()[0] == 'mu_Ss0', "Need mu_Ss0 as first iter var"
	assert iter_vars.keys()[1] == 'Kk', "Need mu_Ss0 as first iter var"
	
	for Kk in iter_vars['Kk']:
		for mu_Ss0 in iter_vars['mu_Ss0']:
			
			# Dummy iter_vars: both will be overwritten
			iter_var_idxs = sp.zeros(len(iter_vars)).astype(int)
			vars_to_pass = compile_all_run_vars(list_dict, iter_var_idxs)
			
			# Overwrite signal strength and complexity manually
			if Kk == 1:
				vars_to_pass['manual_dSs_idxs'] = sp.array([0])
			else:
				sp.random.seed(0)
				vars_to_pass['manual_dSs_idxs'] = \
					sp.random.choice(range(Nn), size=int(Kk), replace=False)
			vars_to_pass['mu_Ss0'] = mu_Ss0
			
			# Set threshold at 10 Hz
			vars_to_pass['NL_threshold'] = 0.05
			vars_to_pass['seed_Kk2'] = seed_Kk2
			obj = four_state_receptor_CS(**vars_to_pass)
			obj = single_encode_CS(obj, list_dict['run_specs'])
			
			fig = fig_primacy_schematic()
			plt.imshow(obj.Yy.reshape(7, 5), cmap=plt.cm.Purples, 
						vmin=vmin, vmax=vmax)
			save_fig('primacy_schematic_Kk=%s_mu_Ss0=%1.3f' % (Kk, mu_Ss0), 
						subdir=data_flag)
	

	# Separate figure for colorbar
	fig = plt.figure()
	fig.set_size_inches(1, 5)
	ax1 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
	norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	ticks = [vmin, vmax]
	cbar = mpl.colorbar.ColorbarBase(ax1, cmap=plt.cm.Purples,
							norm=norm, ticks=ticks, orientation='vertical')
	cbar.ax.tick_params(labelsize=35)
	cbar.ax.set_yticklabels([r'<%s' % vmin, r'>%s' % vmax])
	cbar.ax.yaxis.set_ticks_position('right')
	save_fig('primacy_schematic_Kk=%s_mu_Ss0=%1.2f' % (Kk, mu_Ss0), 
						subdir=data_flag, tight_layout=False)
	
		
if __name__ == '__main__':
	data_flag = get_flag()
	plot_primacy_schematic(data_flag)