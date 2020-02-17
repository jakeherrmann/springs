
from .structure    import Structure
from .springFileIO import FileVariable, FileFormat

class StructureGroup:
	def __init__(self, structures_indexes=None):
		self.structures_indexes = structures_indexes if structures_indexes is not None else []
		#
		self.structures = None

	@staticmethod
	def get_file_format():
		file_format = FileFormat()
		file_format.fixed_format = False
		file_format.variables.append( FileVariable('structures_indexes','I',-1,[int],True) )
		return file_format