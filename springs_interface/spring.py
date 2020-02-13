
from .node import Node
from .springFileIO import FileVariable, FileFormat
import math

class Spring:
	def __init__(self, node_start=None, node_end=None):
		self.node_start = node_start
		self.node_end   = node_end
		self.rest_length = 1.0
		self.stiffness_tension = [1.0]
		self.stiffness_compression = []
		self.compressible = False
		#
		self.node_start_pointer = None
		self.node_end_pointer = None
		self.force = 0.0
		self.strain = 0.0
		self.broken = False
		self.repairable = True
		#
		self.adjacent_nodes_pointers = []
		self.adjacent_springs_pointers = []

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
		file_format.variables.append( FileVariable('node_start'           ,'I',1  , int   ,True) )
		file_format.variables.append( FileVariable('node_end'             ,'I',1  , int   ,True) )
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

	def calc_force(self):
		if self.broken:
			self.strain = 0.0
			self.force = 0.0

		else:
			position_start = self.node_start_pointer.position
			position_end   = self.node_end_pointer.position
			delta_position = [ x_end-x_start for x_start, x_end in zip(position_start, position_end) ]
			length = math.sqrt(sum([ dx*dx for dx in delta_position ]))
			delta_length = length - self.rest_length
			self.strain = delta_length / self.rest_length
			self.force = 0.0
			if delta_length > 0.0:
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

