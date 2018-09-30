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
import matplotlib
import matplotlib.pyplot as plt
sys.path.append('../../shared_src')
from save_load_figure_data import save_fig

# The location of the source code for CS-variability-adaptation is listed
# in the ../../shared_src/local_methods file within src_dir()
from local_methods import src_dir, def_data_dir
sys.path.append(src_dir())


def plot_signal_stats(sig_file, whf_thresh=1, sig_offset=0.0, sig_mult=100):
	"""
	Signal multiplier should come from appropriate specs file. 
	For manuscript, using jul18_power_law_EA_temporal_bkgrnd_step_mult=0.3.txt, 
	for which the multiplier is 150
	"""
		
	data_dir = def_data_dir()
	data = sp.loadtxt('%s/signal_traces/%s.dat' % (data_dir, sig_file))
	
	# Clip array to desired section
	Tt = data[:, 0]
	signal = data[:, 1]
	signal += sig_offset
	signal *= sig_mult
	
	# Signal with whiffs highlighted
	whf_bin = 1.*(signal > whf_thresh)
	whf_beg = Tt[sp.where(sp.diff(whf_bin) == 1)[0]]
	whf_end = Tt[sp.where(sp.diff(whf_bin) == -1)[0]]
	
	if whf_bin[0] == 1:
		whf_beg = sp.hstack((sp.zeros(1), whf_beg))
	if len(whf_beg) > len(whf_end):
		whf_end = sp.hstack((whf_end, Tt[-1]))
	if len(whf_beg) < len(whf_end):
		whf_beg = sp.hstack((Tt[0], whf_beg))
	
	whf_durs = []
	for nWhf in range(len(whf_beg)):
		width = whf_end[nWhf] - whf_beg[nWhf]
		whf_durs.append(width)
	print (len(whf_durs))
	hist, bins = sp.histogram(sp.log(whf_durs)/sp.log(10), bins=50, density=1)
	log_prob = sp.log(hist)/sp.log(10)
	plt.scatter((bins[1:] + bins[:-1])/2., log_prob, marker='+')
	
	# -3/2 fit line
	xr = sp.arange(-3.0, -1.25, 0.01)
	yint = -3.4
	plt.plot(xr, xr*-1.5 + yint, color='r')
	
	plt.xlim(-3.5, -1.0)
	#plt.ylim(-3.0, -1.5)
	plt.show()
	
	
	
if __name__ == '__main__':
	plot_signal_stats(sys.argv[1])