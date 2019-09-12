
from .node         import Node
from .spring       import Spring
from .springFileIO import FileVariable, FileFormat

class Structure:
	def __init__(self, nodes=None, springs=None):
		self.nodes   = nodes   if nodes   is not None else []
		self.springs = springs if springs is not None else []
		#
		self.nodes_pointer = None
		self.springs_pointer = None

	@staticmethod
	def get_file_format():
		file_format = FileFormat()
		file_format.fixed_format = False
		file_format.variables.append( FileVariable('nodes'  ,'I',-1,[int],True) )
		file_format.variables.append( FileVariable('springs','I',-1,[int],True) )
		return file_format

	def connect_structures_nodes_springs(self, spring_network):
		self.nodes_pointer   = [ spring_network.nodes[n]   for n in self.nodes   ]
		self.springs_pointer = [ spring_network.springs[s] for s in self.springs ]
