"""
Object to encode and decode sparse signals using compressed sensing
but with passing a sparse odor signal through a sensory system 
described by a 4-state receptor system. Off and on states 
are distinguished here, and binding kinetics are assumed fast 
enough to leave near the steady state limit. 

The response matrix is assumed to be a linearization of the full
nonlinear response. This linearization is essentially the matrix of 
inverse disassociation constants. This script tests the decoding 
fidelity for various choices of the  mean value of the inverse K. 

Created by Nirag Kadakia at 23:30 07-31-2017
This work is licensed under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
International License. 
To view a copy of this license, 
visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
"""

import scipy as sp
import sys
sys.path.append('../src')
from lin_alg_structs import random_matrix, sparse_vector, \
							sparse_vector_bkgrnd, manual_sparse_vector
from kinetics import bkgrnd_activity, linear_gain, receptor_activity, \
						free_energy, Kk2_samples, Kk2_eval_normal_activity, \
						Kk2_eval_exponential_activity, \
						Kk2_eval_uniform_activity
from optimize import decode_CS, decode_nonlinear_CS
from utils import clip_array


INT_PARAMS = ['Nn', 'Kk', 'Mm', 'seed_Ss0', 'seed_dSs', 'seed_Kk1', 
				'seed_Kk2', 'seed_receptor_activity', 'estimate_full_signal']


class four_state_receptor_CS:	
	"""	
	Object for encoding and decoding a four-state receptor model 
	using compressed sensing
	"""

	def __init__(self, **kwargs):
	
		# Set system parameters
		self.Nn = 50
		self.Kk = 5
		self.Mm = 20

		# Set random seeds
		self.seed_Ss0 = 1
		self.seed_dSs = 1
		self.seed_Kk1 = 1
		self.seed_Kk2 = 1
		self.seed_eps = 1
		self.seed_receptor_activity = 1
		self.seed_adapted_activity = 1
		
		# Randomly chosen sparse signal background and fluctuations
		self.mu_Ss0 = 1.
		self.sigma_Ss0 = 0.001
		self.mu_dSs = 0.3
		self.sigma_dSs = 0.1
		
		# Nonsparse signal deviation, distributed normally over N components
		self.width_dSs = 5
		self.center_dSs = 0
		
		# Manual signals
		self.manual_dSs_idxs = sp.array([0])
		
		# All K1 and K2 from a single Gaussian distribution
		self.mu_Kk1 = 1e4
		self.sigma_Kk1 = 1e3
		self.mu_Kk2 = 1e-3
		self.sigma_Kk2 = 1e-4
		
		# All K1 and K2 from a mixture of 2 Gaussians
		self.mu_Kk1_1 = 1e4
		self.sigma_Kk1_1 = 1e3
		self.mu_Kk1_2 = 1e4
		self.sigma_Kk1_2 = 1e3
		self.mu_Kk2_1 = 1e-3
		self.sigma_Kk2_1 = 1e-4
		self.mu_Kk2_2 = 1e-3
		self.sigma_Kk2_2 = 1e-4
		self.Kk1_p = 0.5
		self.Kk2_p = 0.5
		
		# All K1 and K2 from a single uniform distribution
		self.uniform_Kk1_lo = 1e2
		self.uniform_Kk1_hi = 1e3
		self.uniform_Kk2_lo = 1e-3
		self.uniform_Kk2_hi = 1e-2
		
		# K1 and K2: each receptor from a distinct Gaussian, 
		# with uniform prior on means and sigmas
		self.mu_Kk1_hyper_lo = 1e4
		self.mu_Kk1_hyper_hi = 1e4
		self.sigma_Kk1_hyper_lo = 1e3
		self.sigma_Kk1_hyper_hi = 1e3
		self.mu_Kk2_hyper_lo = 1e-3
		self.mu_Kk2_hyper_hi = 1e-3
		self.sigma_Kk2_hyper_lo = 1e-4
		self.sigma_Kk2_hyper_hi = 1e-4
				
		# Fixed activity distributions for adapted individual odorant response, 
		# where each activity is chosen from mixture of 2 Gaussians
		self.activity_p = 0.5
		self.receptor_tuning_mixture_mu_1 = 0.5
		self.receptor_tuning_mixture_mu_2 = 0.5
		self.receptor_tuning_mixture_sigma_1 = 0.1
		self.receptor_tuning_mixture_sigma_2 = 0.1
		
		# Fixed activity distributions for adapted individual odorant response, 
		# where each activity is chosen from a single uniform distribution
		self.uniform_activity_lo = 1e-5
		self.uniform_activity_hi = 1e-4
 		
		# Fixed activity distributions for adapted individual odorant response, 
		# where each activity is chosen from a normal distribution depending on
		# the receptor. The means for each receptor are chosen normally and the 
		# sigmas are chosen uniformly. 
		self.receptor_tuning_mu_hyper_mu = 0.2
		self.receptor_tuning_mu_hyper_sigma = 0
		self.receptor_tuning_sigma_hyper_lo = 0
		self.receptor_tuning_sigma_hyper_hi = 0.3
		
		# Random free energy statistics
		self.mu_eps = 5.0
		self.sigma_eps = 0.0
		
		# Tuned free energy maximum
		self.normal_eps_tuning_prefactor = 0.0
		self.normal_eps_tuning_width = 1.0
		self.WL_scaling = 0.0
		
		# Fix tuning curve statistics for adapted full activity
		self.adapted_activity_mu = 0.5
		self.adapted_activity_sigma = 0.01
		
		# Estimate full signal or just mu_dSs above background?
		self.estimate_full_signal = True 
		
		# Overwrite variables with passed arguments	
		for key in kwargs:
			if key in INT_PARAMS:
				exec ('self.%s = int(kwargs[key])' % key)
			else:
				exec ('self.%s = kwargs[key]' % key)
			
	def set_sparse_signals(self):
		"""Set random sparse signals"""
	
		self.params_dSs = [self.mu_dSs, self.sigma_dSs]
		self.params_Ss0 = [self.mu_Ss0, self.sigma_Ss0]
		self.dSs, self.idxs = sparse_vector([self.Nn, self.Kk], 
											self.params_dSs, 
											seed=self.seed_dSs)
		
		# Ss0 is the ideal (learned) background stimulus without noise
		self.Ss0, self.Ss0_noisy = sparse_vector_bkgrnd([self.Nn, self.Kk], 
														self.idxs, 
														self.params_Ss0, 
														seed=self.seed_Ss0)
		
		# The true signal, including background noise
		self.Ss = self.dSs + self.Ss0_noisy
	
	def set_manual_signals(self):
		"""
		Set manually-selected sparse signals. 
		"""
		
		self.params_dSs = [self.mu_dSs, self.sigma_dSs]
		self.params_Ss0 = [self.mu_Ss0, self.sigma_Ss0]
		self.idxs = self.manual_dSs_idxs
		self.dSs = manual_sparse_vector(self.Nn, self.manual_dSs_idxs, 
										self.params_dSs, seed=self.seed_dSs)
		
		# Ss0 is the ideal (learned) background stimulus without noise
		self.Ss0, self.Ss0_noisy = sparse_vector_bkgrnd([self.Nn, self.Kk], 
														self.manual_dSs_idxs, 
														self.params_Ss0, 
														seed=self.seed_Ss0)
		
		# The true signal, including background noise
		self.Ss = self.dSs + self.Ss0_noisy
	
	def set_adapted_free_energy(self):
		"""
		Set free energy based on adapted activity activity.
		"""
		
		activity_stats = [self.adapted_activity_mu, self.adapted_activity_sigma]
		adapted_activity = random_matrix([self.Mm], params=activity_stats, 
									seed=self.seed_adapted_activity)
		self.eps = free_energy(self.Ss0, self.Kk1, self.Kk2, adapted_activity)
	
	def set_normal_free_energy(self):
		"""
		Set free energy as a function of odorant; normal tuning curve.
		"""
		
		self.eps_base = self.mu_eps + self.normal_eps_tuning_prefactor* \
						sp.exp(-(1.*sp.arange(self.Mm))**2.0/(2.0* \
						self.normal_eps_tuning_width)**2.0)
		
		self.eps = self.WL_scaling*sp.log(self.mu_dSs) + self.eps_base
		
	def set_random_free_energy(self):
		"""
		Set free energy as random vector.
		"""
		
		self.eps = random_matrix([self.Mm], [self.mu_eps, self.sigma_eps], 
									seed = self.seed_eps)

	def set_mixture_Kk(self, clip=True):
		"""
		Set K1 and K2 matrices where each receptor response is chosen from 
		a Gaussian mixture with stats mu_Kk1_1, sigma_Kk2_1, 
		mu_Kk1_2, sigma_Kk2_2.
		"""
		
		assert 0 <= self.Kk1_p <= 1., "Kk1 Mixture ratio must be between 0 and 1"
		assert 0 <= self.Kk2_p <= 1., "Kk2 Mixture ratio must be between 0 and 1"
		
		self.Kk1 = sp.zeros((self.Mm, self.Nn))
		self.Kk2 = sp.zeros((self.Mm, self.Nn))
		
		num_comp1 = int(self.Kk1_p*self.Mm)
		num_comp2 = self.Mm - num_comp1
		params_Kk1_1 = [self.mu_Kk1_1, self.sigma_Kk1_1]
		params_Kk1_2 = [self.mu_Kk1_2, self.sigma_Kk1_2]
		self.Kk1[:num_comp1, :] = random_matrix([num_comp1, self.Nn], 
										params_Kk1_1, seed = self.seed_Kk1)
		self.Kk1[num_comp1:, :] = random_matrix([num_comp2, self.Nn], 
										params_Kk1_2, seed = self.seed_Kk1)
		
		num_comp1 = int(self.Kk2_p*self.Mm)
		num_comp2 = self.Mm - num_comp1
		params_Kk2_1 = [self.mu_Kk2_1, self.sigma_Kk2_1]
		params_Kk2_2 = [self.mu_Kk2_2, self.sigma_Kk2_2]
		
		self.Kk2[:num_comp1, :] = random_matrix([num_comp1, self.Nn], 
										params_Kk2_1, seed = self.seed_Kk2)		
		self.Kk2[num_comp1:, :] = random_matrix([num_comp2, self.Nn], 
										params_Kk2_2, seed = self.seed_Kk2)
		
		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
			
									
	def set_normal_Kk(self, clip=True):	
		"""
		Set K1 and K2 where each receptor from a distinct Gaussian with 
		uniform prior on means and sigmas (mu_Kk1_hyper_lo, mu_Kk1_hyper_hi, 
		sigma_Kk1_hyper_lo, sigma_Kk1_hyper_hi, etc.)
		"""
		
		Kk1_mus = random_matrix([self.Mm], params=[self.mu_Kk1_hyper_lo, 
							self.mu_Kk1_hyper_hi], sample_type='uniform',
							seed=self.seed_Kk1)
		Kk1_sigmas = random_matrix([self.Mm], params=[self.sigma_Kk1_hyper_lo, 
							self.sigma_Kk1_hyper_hi], sample_type='uniform',
							seed=self.seed_Kk1)
		Kk2_mus = random_matrix([self.Mm], params=[self.mu_Kk2_hyper_lo, 
							self.mu_Kk2_hyper_hi], sample_type='uniform',
							seed=self.seed_Kk2)
		Kk2_sigmas = random_matrix([self.Mm], params=[self.sigma_Kk2_hyper_lo,
							self.sigma_Kk2_hyper_hi], sample_type='uniform',
							seed=self.seed_Kk2)
		
		self.Kk1 = random_matrix([self.Mm, self.Nn], [Kk1_mus, Kk1_sigmas], 
									sample_type='rank2_row_gaussian', 
									seed = self.seed_Kk1)
		self.Kk2 = random_matrix([self.Mm, self.Nn], [Kk2_mus, Kk2_sigmas],
									sample_type='rank2_row_gaussian', 
									seed = self.seed_Kk2)
		

		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
	
	def set_uniform_ordered_Kk(self, clip=True):
		"""
		Set K1 and K2 where each receptor from same Gaussian, and the 
		tuning curves are ordered such that row one is centered at N1, etc.
		"""
		self.Kk1 = random_matrix([self.Mm, self.Nn], [self.uniform_Kk1_lo, 
								self.uniform_Kk1_hi], sample_type='uniform',
								seed = self.seed_Kk1)
		
		self.Kk2 = random_matrix([self.Mm, self.Nn], [self.uniform_Kk2_lo, 
								self.uniform_Kk2_hi], sample_type='uniform', 
								seed = self.seed_Kk2)
		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
							
		for iM in range(self.Mm):
			self.Kk2[iM, :] = sp.sort(self.Kk2[iM, :])
			tmp = sp.hstack((self.Kk2[iM, :], self.Kk2[iM, ::-1]))
			self.Kk2[iM, :] = tmp[::2]
			self.Kk2[iM, :] = sp.roll(self.Kk2[iM, :], iM*int(self.Nn/self.Mm))
			
	def set_uniform_Kk(self, clip=True):	
		"""
		K1 and K2 are chosen from a uniform distribution.
		"""
		
		self.Kk1 = random_matrix([self.Mm, self.Nn], [self.uniform_Kk1_lo, 
								self.uniform_Kk1_hi], sample_type='uniform',
								seed = self.seed_Kk1)
		self.Kk2 = random_matrix([self.Mm, self.Nn], [self.uniform_Kk2_lo, 
								self.uniform_Kk2_hi], sample_type='uniform', 
								seed = self.seed_Kk2)
		
	def set_Kk2_normal_activity(self, clip=True, **kwargs):
		"""
		Fixed activity distributions for adapted individual odorant response, 
		where each activity is chosen from a normal distribution depending on
		the receptor. The means for each receptor are chosen normally and the 
		sigmas are chosen uniformly. Kk2 is then derived from these choices,
		while Kk1 is chosen from a single Gaussian distribution. The function
		allows for overriding of mu_eps and mu_Ss0 if system has adapted to 
		another background state.
		"""
		
		matrix_shape = [self.Mm, self.Nn]
		
		params_Kk1 = [self.mu_Kk1, self.sigma_Kk1]
		self.Kk1 = random_matrix(matrix_shape, params_Kk1, seed=self.seed_Kk1)
	
		mu_Ss0 = self.mu_Ss0
		mu_eps = self.mu_eps
		for key in kwargs:
			exec ('%s = kwargs[key]' % key)
		
		mu_stats = [self.receptor_tuning_mu_hyper_mu, 
					self.receptor_tuning_mu_hyper_sigma]
		sigma_stats = [self.receptor_tuning_sigma_hyper_lo, 
						self.receptor_tuning_sigma_hyper_hi]
		activity_mus = random_matrix([self.Mm], params=mu_stats, 
										sample_type='normal',
										seed=self.seed_receptor_activity)
		activity_sigmas = random_matrix([self.Mm], params=sigma_stats, 
										sample_type='uniform',
										seed=self.seed_receptor_activity)
		
		self.Kk2 = Kk2_eval_normal_activity(matrix_shape, activity_mus, 
										activity_sigmas, mu_Ss0, mu_eps, 
										self.seed_Kk2)
		
		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
	
	def set_Kk2_uniform_activity(self, clip=True, **kwargs):
		"""
		Fixed activity distributions for adapted individual odorant response, 
		where each activity is chosen from a single uniform distribution.
		"""
		
		matrix_shape = [self.Mm, self.Nn]
		
		params_Kk1 = [self.mu_Kk1, self.sigma_Kk1]
		self.Kk1 = random_matrix(matrix_shape, params_Kk1, seed=self.seed_Kk1)
	
		params_Kk2 = [self.uniform_activity_lo, self.uniform_activity_hi]
		self.Kk2 = Kk2_eval_uniform_activity(matrix_shape, params_Kk2, 
											self.mu_Ss0, self.mu_eps, 
											self.seed_Kk2)
		
		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
	
	def set_Kk2_normal_activity_mixture(self, clip=True, **kwargs):
		"""
		Fixed activity distributions for adapted individual odorant response, 
		where each activity is chosen from a Gaussian mixture.
		"""
		
		matrix_shape = [self.Mm, self.Nn]
		
		params_Kk1 = [self.mu_Kk1, self.sigma_Kk1]
		self.Kk1 = random_matrix(matrix_shape, params_Kk1, seed=self.seed_Kk1)
	
		mu_Ss0 = self.mu_Ss0
		mu_eps = self.mu_eps
		for key in kwargs:
			exec ('%s = kwargs[key]' % key)
		
		assert 0 <= self.activity_p <= 1., "Mixture ratio must be between 0 and 1"
		
		activity_mus = sp.zeros(self.Mm)
		activity_sigmas = sp.zeros(self.Mm)
		
		num_comp1 = int(self.activity_p*self.Mm)
		num_comp2 = self.Mm - num_comp1
		
		activity_mus[:num_comp1] = self.receptor_tuning_mixture_mu_1
		activity_mus[num_comp1:] = self.receptor_tuning_mixture_mu_2
		activity_sigmas[:num_comp1] = self.receptor_tuning_mixture_sigma_1
		activity_sigmas[num_comp1:] = self.receptor_tuning_mixture_sigma_2
		
		self.Kk2 = Kk2_eval_normal_activity(matrix_shape, activity_mus, 
											activity_sigmas, mu_Ss0, mu_eps, 
											self.seed_Kk2)
		
		if clip == True:
			array_dict = clip_array(dict(Kk1 = self.Kk1, Kk2 = self.Kk2))
			self.Kk1 = array_dict['Kk1']
			self.Kk2 = array_dict['Kk2']
	
	
	def set_measured_activity(self):
		# True receptor activity
		self.Yy = receptor_activity(self.Ss, self.Kk1, self.Kk2, self.eps)
		
		# Learned background activity only utilizes average background signal 
		self.Yy0 = bkgrnd_activity(self.Ss0, self.Kk1, self.Kk2, self.eps)
		
		# Measured response above background
		self.dYy = self.Yy - self.Yy0

	def set_linearized_response(self):
		# Linearized response can only use the learned background, or ignore
		# that knowledge
		self.Rr = linear_gain(self.Ss0, self.Kk1, self.Kk2, self.eps)
			
	def encode_normal_activity(self, **kwargs):
		# Run all functions to encode when activity is normally distributed.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_Kk2_normal_activity(**kwargs)
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_uniform_activity(self, **kwargs):
		# Run all functions to encode when activity is uniformly distributed.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_Kk2_uniform_activity()
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_normal_activity_mixture(self, **kwargs):
		# Run all functions to encode when activity arises from a mixture.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_Kk2_normal_activity_mixture(**kwargs)
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_normal_Kk(self):
		# Run all functions to encode when K matrices are Gaussian.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_normal_Kk()
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_uniform_Kk(self):
		# Run all functions to encode when K matrices are uniform.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_uniform_Kk()
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_mixture_Kk(self):
		# Run all functions to encode when full activity of each receptor 
		# is from a Gaussian mixture.
		self.set_sparse_signals()
		self.set_random_free_energy()
		self.set_mixture_Kk()
		self.set_measured_activity()
		self.set_linearized_response()
		
	def encode_adapted_normal_activity(self):
		# Run all functions to encode when full activity of each receptor 
		# is normal.
		self.set_sparse_signals()
		self.set_normal_Kk()
		self.set_adapted_free_energy()
		self.set_measured_activity()
		self.set_linearized_response()

	def encode_normal_signal_adapted_energy(self):
		self.set_normal_signals()
		self.set_uniform_Kk()
		self.set_normal_free_energy()
		self.set_measured_activity()
		self.set_linearized_response()
	
	def encode_manual_signal_adapted_energy(self):
		self.set_manual_signals()
		self.set_uniform_Kk()
		self.set_normal_free_energy()
		self.set_measured_activity()
		self.set_linearized_response()
	
	def decode(self):
		self.dSs_est = decode_CS(self.Rr, self.dYy)	
		
	def decode_nonlinear(self):
		self.dSs_est = decode_nonlinear_CS(self)