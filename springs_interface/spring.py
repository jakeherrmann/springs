
from .node import Node
from .springFileIO import FileVariable, FileFormat
import math

class Spring:
	def __init__(self, node_start_index=None, node_end_index=None):
		self.node_start_index = node_start_index
		self.node_end_index   = node_end_index
		self.rest_length = 1.0
		self.stiffness_tension = [1.0]
		self.stiffness_compression = []
		self.compressible = False
		#
		self.node_start = None
		self.node_end = None
		self.force = 0.0
		self.strain = 0.0
		self.broken = False
		self.repairable = True
		#
		self.adjacent_nodes = []
		self.adjacent_springs = []
		self.structures = []

	@staticmethod
	def get_file_format(precision, num_stiffness_tension, num_stiffness_compression, use_solver_format=False):
		if precision=='float':
			FPP = 'f'
		else:
			FPP = 'd'
		NST = num_stiffness_tension
		NSC = num_stiffness_compression
		file_format = FileFormat()
		file_format.fixed_format = True
		file_format.variables.append( FileVariable('node_start_index'     ,'I',1  , int   ,True) )
		file_format.variables.append( FileVariable('node_end_index'       ,'I',1  , int   ,True) )
		file_format.variables.append( FileVariable('rest_length'          ,FPP,1  , float ,True) )
		file_format.variables.append( FileVariable('stiffness_tension'    ,FPP,NST,[float],True) )
		file_format.variables.append( FileVariable('stiffness_compression',FPP,NSC,[float],True) )
		file_format.variables.append( FileVariable('compressible'         ,'?',1  , bool  ,True) )
		if not use_solver_format:
			file_format.variables.append( FileVariable('force'     ,FPP,1,float,True) )
			file_format.variables.append( FileVariable('strain'    ,FPP,1,float,True) )
			file_format.variables.append( FileVariable('broken'    ,'?',1,bool ,True) )
			file_format.variables.append( FileVariable('repairable','?',1,bool ,True) )
		return file_format

	def effective_stiffness(self):
		if self.broken:
			return 0.0
		else:
			position_start = self.node_start.position
			position_end   = self.node_end.position
			delta_position = [ x_end-x_start for x_start, x_end in zip(position_start, position_end) ]
			length = math.sqrt(sum([ dx*dx for dx in delta_position ]))
			delta_length = length - self.rest_length
			total_stiffness = 0.0
			if delta_length >= 0.0:
				delta_length_pow = 1.0
				for k in self.stiffness_tension:
					total_stiffness += k * delta_length_pow
					delta_length_pow *= delta_length
			elif self.compressible:
				delta_length = abs(delta_length)
				delta_length_pow = 1.0
				for k in self.stiffness_compression:
					total_stiffness += k * delta_length_pow
					delta_length_pow *= delta_length
			return total_stiffness

	def calc_force(self):
		if self.broken:
			self.strain = 0.0
			self.force = 0.0
		else:
			position_start = self.node_start.position
			position_end   = self.node_end.position
			delta_position = [ x_end-x_start for x_start, x_end in zip(position_start, position_end) ]
			length = math.sqrt(sum([ dx*dx for dx in delta_position ]))
			delta_length = length - self.rest_length
			self.strain = delta_length / self.rest_length
			self.force = 0.0
			if delta_length >= 0.0:
				delta_length_pow = delta_length
				for k in self.stiffness_tension:
					self.force += k * delta_length_pow
					delta_length_pow *= delta_length
			elif self.compressible:
				delta_length = abs(delta_length)
				delta_length_pow = delta_length
				for k in self.stiffness_compression:
					self.force += k * delta_length_pow
					delta_length_pow *= delta_length

