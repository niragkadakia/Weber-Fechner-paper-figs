"""
Plot step stimulus cartoon for figure.


Created by Nirag Kadakia at 14:37 04-09-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
sys.path.append('../../shared_src')
from plot_formats import fig_temporal_kernel
from save_load_figure_data import save_fig

stim_val = 0.7

fig = fig_temporal_kernel()
fig.set_size_inches(10, 2)
stim_Tt = sp.arange(-2, 7, 1e-3)
stim = stim_val*10*(stim_Tt > 0)
stim = gaussian_filter(stim, sigma=30)
ax = plt.axes(frameon=False)
plt.plot(stim, lw=5, color='k')
plt.ylim(-0.2, 12)
ax.axis('off')
save_fig('step_stim_cartoon_%s' % stim_val, subdir='./')
