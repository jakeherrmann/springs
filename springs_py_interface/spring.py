
from .node import Node
from .springFileIO import FileVariable, FileFormat
import math

class Spring:
	def __init__(self, node_start_index=None, node_end_index=None):
		self.node_start_index = node_start_index
		self.node_end_index   = node_end_index
		self.force_length_type_tension = 1
		self.force_length_type_compression = 0
		self.rest_length = 1.0
		self.force_length_parameters_tension = [1.0]
		self.force_length_parameters_compression = []
		#
		self.node_start = None
		self.node_end = None
		self.strain = 0.0
		self.force = 0.0
		self.energy = 0.0
		self.effective_spring_constant = 0.0
		self.broken = False
		self.repairable = True
		#
		self.adjacent_nodes = []
		self.adjacent_springs = []
		self.structures = []

	@staticmethod
	def get_file_format(precision, use_solver_format=False):
		if precision=='float':
			FPP = 'f'
		else:
			FPP = 'd'
		file_format = FileFormat()
		file_format.fixed_format = False
		file_format.variables.append( FileVariable('node_start_index'                   ,'I', 1, int   ,True) )
		file_format.variables.append( FileVariable('node_end_index'                     ,'I', 1, int   ,True) )
		file_format.variables.append( FileVariable('force_length_type_tension'          ,'B', 1, int   ,True) )
		file_format.variables.append( FileVariable('force_length_type_compression'      ,'B', 1, int   ,True) )
		file_format.variables.append( FileVariable('rest_length'                        ,FPP, 1, float ,True) )
		file_format.variables.append( FileVariable('force_length_parameters_tension'    ,FPP,-1,[float],True) )
		file_format.variables.append( FileVariable('force_length_parameters_compression',FPP,-1,[float],True) )
		if not use_solver_format:
			file_format.variables.append( FileVariable('strain'                   ,FPP,1,float,True) )
			file_format.variables.append( FileVariable('force'                    ,FPP,1,float,True) )
			file_format.variables.append( FileVariable('energy'                   ,FPP,1,float,True) )
			file_format.variables.append( FileVariable('effective_spring_constant',FPP,1,float,True) )
			file_format.variables.append( FileVariable('broken'                   ,'?',1,bool ,True) )
			file_format.variables.append( FileVariable('repairable'               ,'?',1,bool ,True) )
		return file_format

	def calc_force(self):
		self.strain = 0.0
		self.force = 0.0
		self.energy = 0.0
		self.effective_spring_constant = 0.0
		if self.broken:
			return
		else:
			position_start = self.node_start.position
			position_end   = self.node_end.position
			delta_position = [ x_end-x_start for x_start, x_end in zip(position_start, position_end) ]
			length = math.sqrt(sum([ dx*dx for dx in delta_position ]))
			delta_length = length - self.rest_length
			self.strain = delta_length / self.rest_length
			#
			if delta_length == 0.0:
				self.effective_spring_constant = self.spring_constant_at_rest()
				return
			elif delta_length > 0.0:
				DL = delta_length
				force_length_type = self.force_length_type_tension
				force_length_parameters = self.force_length_parameters_tension
				force_sign = +1.0
			elif delta_length < 0.0:
				DL = abs(delta_length)
				force_length_type = self.force_length_type_compression
				force_length_parameters = self.force_length_parameters_compression
				force_sign = -1.0
			#
			if force_length_type == 0: #none
				return
			elif force_length_type == 1: #polynomial
				self.force                     = sum([ k*math.pow(DL,p)         for p, k in enumerate(force_length_parameters,start=1) ])
				self.energy                    = sum([ k*math.pow(DL,p+1)/(p+1) for p, k in enumerate(force_length_parameters,start=1) ])
				self.effective_spring_constant = sum([ k*math.pow(DL,p-1)/p     for p, k in enumerate(force_length_parameters,start=1) ])
			elif force_length_type == 2: #exponential
				self.force                     = sum([ k*math.expm1(p*DL)          for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
				self.energy                    = sum([ k*((math.expm1(p*DL)/p)-DL) for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
				self.effective_spring_constant = sum([ k*p*math.exp(p*DL)          for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
			elif force_length_type == 3: #powerlaw
				self.force                     = sum([ k*math.pow(DL,p)             for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
				self.energy                    = sum([ k*math.pow(DL,p+1.0)/(p+1.0) for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
				self.effective_spring_constant = sum([ k*math.pow(DL,p-1.0)/p       for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
			self.force *= force_sign
		return

	def spring_constant_at_rest(self):
		if self.broken:
			return 0.0
		else:
			if self.force_length_type_tension != 0 and len(self.force_length_parameters_tension) > 0:
				force_length_parameters = self.force_length_parameters_tension
				force_length_type       = self.force_length_type_tension
			elif self.force_length_type_compression != 0 and len(self.force_length_parameters_compression) > 0:
				force_length_parameters = self.force_length_parameters_compression
				force_length_type       = self.force_length_type_compression
			else:
				return 0.0
			#
			if force_length_type == 0: #none
				return 0.0
			elif force_length_type == 1: #polynomial
				return force_length_parameters[0]
			elif force_length_type == 2: #exponential
				return sum([ k*p for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
			elif force_length_type == 3: #powerlaw
				return sum([ 0.0 if p>1.0 else k if p==1.0 else 1.0e20 for k, p in zip(force_length_parameters[0::2],force_length_parameters[1::2]) ])
			else:
				return 0.0

	def modify_spring_constant(self, multiplier, modify_tension=True, modify_compression=True):
		if modify_tension:
			if self.force_length_type_tension == 0: #none
				pass
			elif self.force_length_type_tension == 1: #polynomial
				self.force_length_parameters_tension = [ k*multiplier for k in self.force_length_parameters_tension ]
			elif self.force_length_type_tension == 2: #exponential
				self.force_length_parameters_tension = [ k*multiplier if i%2==0 else k for i, k in enumerate(self.force_length_parameters_tension) ]
			elif self.force_length_type_tension == 3: #powerlaw
				self.force_length_parameters_tension = [ k*multiplier if i%2==0 else k for i, k in enumerate(self.force_length_parameters_tension) ]
		if modify_compression:
			if self.force_length_type_compression == 0: #none
				pass
			elif self.force_length_type_compression == 1: #polynomial
				self.force_length_parameters_compression = [ k*multiplier for k in self.force_length_parameters_compression ]
			elif self.force_length_type_compression == 2: #exponential
				self.force_length_parameters_compression = [ k*multiplier if i%2==0 else k for i, k in enumerate(self.force_length_parameters_compression) ]
			elif self.force_length_type_compression == 3: #powerlaw
				self.force_length_parameters_compression = [ k*multiplier if i%2==0 else k for i, k in enumerate(self.force_length_parameters_compression) ]
		return
