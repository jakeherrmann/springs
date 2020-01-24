from .agent import Agent

import springs_interface as spr
import math

class Agent_Fibroblast(Agent):
	def __init__(self, spring_index=None, spring=None):
		self.actions = [ self.maintain ]
		self.spring_index = spring_index
		self.spring = spring
		self.strain = 0.0
		self.strain_energy_rate = 0.0
		self.degrade_rate_basal = 0.100 ;
		self.repair_rate_basal  = 0.100 ;
		self.degrade_coefs = [ 1.0000 , 10.4350 ,  3.5173 , 13.0130 , 0.3709 ]
		self.repair_coefs  = [ 0.3968 ,  4.5944 , 45.2080 ]

	def repair(self, strain_energy_rate):
		repair_rate = ( self.repair_coefs[0]
			          + (self.repair_coefs[1]-self.repair_coefs[0])*(-math.expm1(-self.repair_coefs[2]*strain_energy_rate)) )
		repair_rate *= self.repair_rate_basal
		return repair_rate

	def degrade(self, strain):
		degrade_rate = ( self.degrade_coefs[0]*math.exp(-self.degrade_coefs[1]*strain)
			           + self.degrade_coefs[2]/(1.0+math.exp(-self.degrade_coefs[3]*(strain-self.degrade_coefs[4]))) )
		degrade_rate *= self.degrade_rate_basal
		return degrade_rate

	def maintain(self, time_step):
		if not self.spring.broken:
			K = self.spring.stiffness_tension[0]
			if K>0.0:
				# strain = self.spring.strain
				# if self.strain_prev is None:
				# 	self.strain_prev = strain
				# strain_energy_rate = ((strain-self.strain_prev)/time_step)**2 #incorrect
				D = self.degrade(self.strain)
				R = self.repair(self.strain_energy_rate)
				delta_K = R - D*math.sqrt(K)
				delta_K = delta_K if delta_K > -K else -K
				# self.strain_prev = strain
				ratio_K = 1.0 + delta_K/K
				self.spring.stiffness_tension = [ k*ratio_K for k in self.spring.stiffness_tension ]
