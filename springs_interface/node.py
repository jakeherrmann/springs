
from .springFileIO import FileVariable, FileFormat

class Node:
	def __init__(self, num_dimensions):
		self.position   = [0.0] * num_dimensions
		self.force      = [0.0] * num_dimensions
		self.fixed      = False
		self.referenced = True
		#
		self.adjacent_nodes = []
		self.adjacent_springs = []
		self.structures = []

	@staticmethod
	def get_file_format(num_dimensions, precision, use_solver_format=False):
		if precision=='float':
			FPP = 'f'
		else:
			FPP = 'd'
		ND = num_dimensions
		file_format = FileFormat()
		file_format.fixed_format = True
		file_format.variables.append( FileVariable('position',FPP,ND,[float],True) )
		file_format.variables.append( FileVariable('force'   ,FPP,ND,[float],True) )
		file_format.variables.append( FileVariable('fixed'   ,'?',1 , bool  ,True) )
		if not use_solver_format:
			file_format.variables.append( FileVariable('referenced','?',1,bool,True) )
		return file_format