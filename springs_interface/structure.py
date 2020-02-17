
from .node         import Node
from .spring       import Spring
from .springFileIO import FileVariable, FileFormat

class Structure:
	def __init__(self, nodes_indexes=None, springs_indexes=None):
		self.nodes_indexes   = nodes_indexes   if nodes_indexes   is not None else []
		self.springs_indexes = springs_indexes if springs_indexes is not None else []
		#
		self.nodes = None
		self.springs = None
		self.adjacent_structures = []

	@staticmethod
	def get_file_format():
		file_format = FileFormat()
		file_format.fixed_format = False
		file_format.variables.append( FileVariable('nodes_indexes'  ,'I',-1,[int],True) )
		file_format.variables.append( FileVariable('springs_indexes','I',-1,[int],True) )
		return file_format