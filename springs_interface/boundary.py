
from .node import Node
from .springFileIO import FileVariable, FileFormat

class Boundary:

	def __init__(self, num_dimensions=2):
		self.nodes_indexes = []
		self.fixed = False
		self.force_directions = []
		self.force_magnitudes = []
		self.outward_direction  = [0.0] * num_dimensions
		self.displacements = []
		#
		self.nodes = None

	@staticmethod
	def get_file_format(num_dimensions):
		ND = num_dimensions
		file_format = FileFormat()
		file_format.fixed_format = False
		file_format.variables.append( FileVariable('nodes_indexes'    ,'I',-1, [int]   ,True) )
		file_format.variables.append( FileVariable('fixed'            ,'?', 1,  bool   ,True) )
		file_format.variables.append( FileVariable('force_directions' ,'d',-1,[[float]],True) )
		file_format.variables.append( FileVariable('force_magnitudes' ,'d',-1, [float] ,True) )
		file_format.variables.append( FileVariable('outward_direction','d',ND, [float] ,True) )
		file_format.variables.append( FileVariable('displacements'    ,'d',-1,[[float]],True) )
		return file_format

	def reset_conditions(self, which_conditions):
		if which_conditions=='all':
			which_conditions = ('force_magnitudes', 'displacements')
		elif isinstance(which_conditions, str):
			which_conditions = (which_conditions,)
		for c in which_conditions:
			if hasattr(self, c):
				setattr(self, c, [])

	def apply_conditions(self):
		for n, node in enumerate(self.nodes):
			node.fixed = node.fixed or self.fixed
		if self.force_magnitudes:
			for n, node in enumerate(self.nodes):
				force = [ d*self.force_magnitudes[n] for d in self.force_directions[n] ]
				node.force = [ fn+fb for fn, fb in zip(node.force, force) ]
		if self.displacements:
			for node, displacement in zip(self.nodes, self.displacements):
				node.position = [ pn+db for pn, db in zip(node.position, displacement) ]
				