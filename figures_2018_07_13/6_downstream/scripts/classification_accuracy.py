"""
Plot error of classification tasks calculated from tensor flow 
optimization, as a function of the number of signals

Created by Nirag Kadakia at 12:40 07-15-2018
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""


import scipy as sp
import sys
import gc
import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig
from plot_formats import fig_classification_accuracy

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir
sys.path.append(src_dir())
from utils import get_flags
from load_specs import read_specs_file
from load_data import load_aggregated_object_list


def plot_accuracy_vs_num_signals(data_flags):

	fig = fig_classification_accuracy()
	
	cmap = plt.cm.viridis
	colors = [0, .25, 0.55, 0.8]
	lws = [3, 3, 3, 3]
	lss = ['-', '-', '-', '-']
	
	for iFlag, data_flag in enumerate(data_flags):
		list_dict = read_specs_file(data_flag)
		iter_vars = list_dict['iter_vars']
		iter_vars_dims = []
		for iter_var in list_dict['iter_vars']:
			iter_vars_dims.append(len(list_dict['iter_vars'][iter_var]))	
		
		assert len(iter_vars) == 2, "Need 2 iter_vars"
		iter_var_names = ['Jj_mask_seed', 'num_signals']
		for iName, name in enumerate(iter_var_names):
			assert list(iter_vars.keys())[iName] == name, "%sth variable "\
				"must have name %s" % (iName, name)
		Jj_mask_seed = iter_vars['Jj_mask_seed']
		num_signals = iter_vars['num_signals']
		
		obj_list = load_aggregated_object_list(iter_vars_dims, data_flag)
		obj_list = sp.reshape(obj_list, iter_vars_dims)
		avg_accuracies = sp.zeros(len(num_signals))
		for iNum in range(len(num_signals)):
			sum_acc = 0
			for seed in range(len(Jj_mask_seed)):
				obj = obj_list[seed, iNum]
				
				# Get number of input patterns for which the index of the 
				# maximum probability (among labels) is the same as the label.
				# `1` is perfect, `0.5` is chance.
				acc = sp.sum(sp.argmax(obj.test_data_calc, axis=1) == 
								sp.argmax(obj.test_data_labels, axis=1))/\
								obj.test_data_calc.shape[0]
				sum_acc += acc
			avg_accuracies[iNum] = sum_acc/len(Jj_mask_seed)
			
		plt.plot(num_signals, avg_accuracies, color=cmap(colors[iFlag]),
					lw=lws[iFlag], linestyle=lss[iFlag])
		
		del (obj_list)
		del (obj)
		del (acc)
		del (list_dict)
		gc.collect()
	
	save_fig('classification_accuracy', subdir=data_flag)
					
if __name__ == '__main__':
	data_flags = get_flags()
	plot_accuracy_vs_num_signals(data_flags)