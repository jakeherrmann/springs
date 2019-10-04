from .agent import Agent

import springs_interface as spr
import math

class Agent_Fibroblast(Agent):
	def __init__(self, location=None):
		self.actions = [ self.maintain ]
		self.location = location
		self.strain_prev = None
		self.degrade_rate_basal = 0.004 ;
		self.repair_rate_basal  = 0.004 ;
		self.degrade_coefs = [ 1.0000 , 10.4350 , 3.5173 , 13.0130 , 0.3709 ]
		self.repair_coefs  = [ 1.0000 ,  2.5173 , 1.2500 ]

	def repair(self, strain_energy_rate):
		repair_rate = ( self.repair_coefs[0]
			          + self.repair_coefs[1]*(-math.expm1(-self.repair_coefs[2]*strain_energy_rate)) )
		repair_rate *= self.repair_rate_basal
		return repair_rate

	def degrade(self, strain):
		degrade_rate = ( self.degrade_coefs[0]*math.exp(-self.degrade_coefs[1]*strain)
			           + self.degrade_coefs[2]/(1.0+math.exp(-self.degrade_coefs[3]*(strain-self.degrade_coefs[4]))) )
		degrade_rate *= self.degrade_rate_basal
		return degrade_rate

	def maintain(self, time_step):
		if not self.location.broken:
			K = self.location.stiffness_tension[0]
			if K>0.0:
				strain = self.location.strain
				if self.strain_prev is None:
					self.strain_prev = strain
				strain_energy_rate = ((strain-self.strain_prev)/time_step)**2

				D = self.degrade(strain)
				R = self.repair(strain_energy_rate)
				delta_K = R - D*math.sqrt(K)
				delta_K = delta_K if delta_K > -K else -K

				self.strain_prev = strain

				ratio_K = 1.0 + delta_K/K
				self.location.stiffness_tension = [ k*ratio_K for k in self.location.stiffness_tension ]
