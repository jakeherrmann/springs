
from .node import Node

class Boundary:

	def __init__(self, num_dimensions=2):
		self.nodes = []
		self.fixed = False
		self.force_directions = []
		self.force_magnitudes = []
		self.displacement_direction  = [0.0] * num_dimensions
		self.displacement_magnitudes = []
		#
		self.nodes_pointer = None

	def apply_conditions(self):
		for n, node in enumerate(self.nodes_pointer):
			node.fixed = node.fixed or self.fixed
		if self.force_magnitudes:
			for n, node in enumerate(self.nodes_pointer):
				force = [ d*self.force_magnitudes[n] for d in self.force_directions[n] ]
				node.force = [ fn+fb for fn, fb in zip(node.force, force) ]
		if self.displacement_magnitudes:
			for n, node in enumerate(self.nodes_pointer):
				displacement = [ d*self.displacement_magnitudes[n] for d in self.displacement_direction ]
				node.position = [ pn+db for pn, db in zip(node.position, displacements) ]
				