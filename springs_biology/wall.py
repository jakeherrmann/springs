import springs_interface as spr

class Wall:
	def __init__(self, structure=None, alveolus=None):
		self.index = None
		self.structure = structure
		self.thickness = None
		self.alveolus = alveolus
		self.agents = []
		self.adjacent_walls = [] # this is difference from adjacent_structures!!  this can distinguish one side vs another.