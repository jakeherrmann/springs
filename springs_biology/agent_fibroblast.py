from .agent import Agent

import springs_interface as spr
import math

class Agent_Fibroblast(Agent):

	enzyme_degrade_coefs = [ 1.0000 , 10.4350 ,  3.5173 , 13.0130 , 0.3709 ]

	def __init__(self, parent=None, wall=None):
		super().__init__(
			parent=parent,
			wall=wall)
		self.actions = [ self.maintain ]
		#
		self.degrade_variable = None
		self.repair_variable = None
		self.strain = 0.0
		self.strain_rate = 0.0
		self.strain_energy_rate = 0.0
		self.stiffness = 0.0
		#
		self.degrade_rate_multiplier = 0.050
		self.repair_rate_multiplier  = 0.050
		self.degrade_coefs = []
		self.repair_coefs  = []

	def set_behavior(self, degrade_behavior, repair_behavior):
		SIGMIN =  1.0
		SIGMAX = 10.0
		EXPMIN = (2.0*SIGMIN)-SIGMAX
		EXPMAX = 2.0*(SIGMAX-SIGMIN)
		if degrade_behavior=='Constant':
			degrade_variable = None;
			degrade_coefs = [-1.0,+1.0,SIGMIN,SIGMIN]
		else:
			degrade_behavior = degrade_behavior.split('_')
			if   degrade_behavior[0]=='StrainAverage'   : degrade_variable = 'strain'            ; SATVAL = 1.0 ;
			elif degrade_behavior[0]=='StrainRate'      : degrade_variable = 'strain_rate'       ; SATVAL = 0.2 ;
			elif degrade_behavior[0]=='StrainEnergyRate': degrade_variable = 'strain_energy_rate'; SATVAL = 0.1 ;
			elif degrade_behavior[0]=='Stiffness'       : degrade_variable = 'stiffness'         ; SATVAL = 3.0 ;
			if   degrade_behavior[1]=='ExponentialIncreasing': degrade_coefs = [-SATVAL,+SATVAL,EXPMIN,SIGMAX]
			elif degrade_behavior[1]=='ExponentialDecreasing': degrade_coefs = [-SATVAL,+SATVAL,EXPMAX,SIGMIN]
			elif degrade_behavior[1]=='SigmoidalIncreasing'  : degrade_coefs = [0.2*SATVAL,SATVAL,SIGMIN,SIGMAX]
			elif degrade_behavior[1]=='SigmoidalDecreasing'  : degrade_coefs = [0.2*SATVAL,SATVAL,SIGMAX,SIGMIN]
		SIGMIN =  1.0
		SIGMAX = 10.0
		EXPMIN = (2.0*SIGMIN)-SIGMAX
		EXPMAX = 2.0*(SIGMAX-SIGMIN)
		if repair_behavior=='Constant':
			repair_variable = None;
			repair_coefs = [-1.0,+1.0,SIGMIN,SIGMIN]
		else:
			repair_behavior = repair_behavior.split('_')
			if   repair_behavior[0]=='StrainAverage'   : repair_variable = 'strain'            ; SATVAL = 1.0 ;
			elif repair_behavior[0]=='StrainRate'      : repair_variable = 'strain_rate'       ; SATVAL = 0.2 ;
			elif repair_behavior[0]=='StrainEnergyRate': repair_variable = 'strain_energy_rate'; SATVAL = 0.1 ;
			elif repair_behavior[0]=='Stiffness'       : repair_variable = 'stiffness'         ; SATVAL = 3.0 ;
			if   repair_behavior[1]=='ExponentialIncreasing': repair_coefs = [-SATVAL,+SATVAL,EXPMIN,SIGMAX]
			elif repair_behavior[1]=='ExponentialDecreasing': repair_coefs = [-SATVAL,+SATVAL,EXPMAX,SIGMIN]
			elif repair_behavior[1]=='SigmoidalIncreasing'  : repair_coefs = [0.2*SATVAL,SATVAL,SIGMIN,SIGMAX]
			elif repair_behavior[1]=='SigmoidalDecreasing'  : repair_coefs = [0.2*SATVAL,SATVAL,SIGMAX,SIGMIN]

	@staticmethod
	def sigmoid(coefs, arg):
		# coefs = [X0,X1,Y0,Y1] for X and Y saturation limits of sigmoid curve
		arg_normalized = ( arg - coefs[0] -0.5*(coefs[1]-coefs[0]) )/(0.5*(coefs[1]-coefs[0]))
		sig_normalized = 1.0/(1.0+math.exp(-6.0*arg_normalized))
		sig = coefs[2] + (coefs[3]-coefs[2])*sig_normalized
		return sig

	@classmethod
	def enzyme_activity(cls, strain):
		enzyme_degrade_rate_low_strain  = cls.enzyme_degrade_coefs[0]*math.exp(-cls.enzyme_degrade_coefs[1]*strain)
		enzyme_degrade_rate_high_strain = cls.enzyme_degrade_coefs[2]/(1.0+math.exp(-cls.enzyme_degrade_coefs[3]*(strain-cls.enzyme_degrade_coefs[4])))
		enzyme_degrade_rate = enzyme_degrade_rate_low_strain + enzyme_degrade_rate_high_strain
		return enzyme_degrade_rate

	def repair(self):
		if self.repair_variable is not None and hasattr(self,self.repair_variable):
			repair_rate = Agent_Fibroblast.sigmoid( self.repair_coefs , getattr(self,self.repair_variable) )
		else:
			repair_rate = 1.0
		repair_rate *= self.repair_rate_multiplier
		return repair_rate

	def degrade(self):
		if self.degrade_variable is not None and hasattr(self,self.degrade_variable):
			degrade_rate = Agent_Fibroblast.sigmoid( self.degrade_coefs , getattr(self,self.degrade_variable) )
		else:
			degrade_rate = 1.0
		degrade_rate *= Agent_Fibroblast.enzyme_activity(self.strain)
		degrade_rate *= self.degrade_rate_multiplier
		return degrade_rate

	def maintain(self, time_step):
		#TODO better to make strain and strain energy rate properties of wall that can be sensed
		for spring in self.wall.structure.springs:
			if not spring.broken:
				K = spring.stiffness_tension[0]
				if K>0.0:
					# strain = spring.strain
					# if self.strain_prev is None:
					# 	self.strain_prev = strain
					# strain_energy_rate = ((strain-self.strain_prev)/time_step)**2 #incorrect
					D = self.degrade()
					R = self.repair()
					delta_K = R - D*math.sqrt(K)
					delta_K = delta_K if delta_K > -K else -K
					# self.strain_prev = strain
					ratio_K = 1.0 + delta_K/K
					spring.stiffness_tension = [ k*ratio_K for k in spring.stiffness_tension ]
