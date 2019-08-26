
from .node import Node
from .springFileIO import FileVariable, FileFormat

class Spring:

	def __init__(self, node_start=None, node_end=None):
		self.node_start = node_start
		self.node_end   = node_end
		self.rest_length = 1.0
		self.stiffness_tension = [1.0]
		self.stiffness_compression = []
		self.compressible = False
		#
		self.force = 0.0
		self.strain = 0.0
		self.broken = False
		self.repairable = True

	@staticmethod
	def get_file_format(precision, num_stiffness_tension, num_stiffness_compression):
		if precision=='float':
			FPP = 'f'
		else:
			FPP = 'd'
		NST = num_stiffness_tension
		NSC = num_stiffness_compression
		file_format = FileFormat()
		file_format.variables.append( FileVariable('node_start'           ,'I',1  ,int    ,True) )
		file_format.variables.append( FileVariable('node_end'             ,'I',1  ,int    ,True) )
		file_format.variables.append( FileVariable('rest_length'          ,FPP,1  ,float  ,True) )
		file_format.variables.append( FileVariable('stiffness_tension'    ,FPP,NST,[float],True) )
		file_format.variables.append( FileVariable('stiffness_compression',FPP,NSC,[float],True) )
		file_format.variables.append( FileVariable('compressible'         ,'?',1  ,bool   ,True) )
		return file_format