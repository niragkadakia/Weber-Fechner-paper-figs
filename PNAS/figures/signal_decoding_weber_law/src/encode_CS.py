"""
Run CS encoding and decoding via four_state_receptor_CS; single iteration.


Created by Nirag Kadakia at 10:30 09-05-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

from four_state_receptor_CS import four_state_receptor_CS


def single_encode_CS(vars_to_pass=dict(), run_specs=dict()):
	"""
	Run CS encoding and decoding via four_state_receptor_CS; single iteration.
	
	Optional args:
		vars_to_pass: dictionary; any overriden arguments to 
						four_state_receptor_CS.
		run_specs: dictionary; parameters of the run.
	"""		
	
	a = four_state_receptor_CS(**vars_to_pass)
	
	if 'run_type' in run_specs.keys():
		val = run_specs['run_type']
		if val[0] == 'normal_Kk':
			a.encode_normal_Kk()
		elif val[0] == 'uniform_Kk':
			a.encode_uniform_Kk()
		elif val[0] == 'mixture_Kk':
			a.encode_mixture_Kk()
		elif val[0] == 'normal_activity':
			a.encode_normal_activity()
		elif val[0] == 'uniform_activity':
			a.encode_uniform_activity()
		elif val[0] == 'normal_activity_fixed_Kk2':
			override_parameters = dict()
			override_parameters['mu_Ss0'] = float(val[1])
			override_parameters['mu_eps'] = float(val[2])
			a.encode_normal_activity(**override_parameters)
		elif val[0] == 'adapted_normal_activity':
			a.encode_adapted_normal_activity()
		elif val[0] == 'normal_activity_mixture':
			a.encode_normal_activity_mixture()
		elif val[0] == 'normal_signal_adapted_energy':
			a.encode_normal_signal_adapted_energy()
		elif val[0] == 'manual_signal_adapted_energy':
			a.encode_manual_signal_adapted_energy()
		else:
			print ('Run specification %s not recognized' % val[0])
			quit()
	else:
		print ('No run type specified, proceeding with normal_activity')
		a.encode_normal_activity()
	
	return a