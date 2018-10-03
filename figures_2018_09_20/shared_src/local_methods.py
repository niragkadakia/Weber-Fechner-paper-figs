"""
Miscellaneous functions that depend on the machine or local
parameters or setups. 

Requires that the paths are changed, and that the name of this file is changed
to "local_methods"

Created by Evan Cudone at 6:38 11-30-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

def src_dir():	
	""" 
	Define where the source code is for CS-variability-adaptation
	"""
	
	data_dir = 'C:/Users/nk479/Documents/code/CS-variability-adaptation/src'
	
	return data_dir
	
def def_data_dir():
	"""
	Define a data directory here for all input and output.
	"""
	
	data_dir = "C:/Users/nk479/Dropbox (emonetlab)/users/" \
				"nirag_kadakia/data/CS-variability-adaptation"

	return data_dir
	
def def_figure_analysis_data_dir():
	"""
	Define a data directory here for all input and output.
	"""
	
	data_dir = "C:/Users/nk479/Dropbox (emonetlab)/users/" \
				"nirag_kadakia/data/CS-variability-adaptation/" \
				"analysis_for_paper"

	return data_dir
	