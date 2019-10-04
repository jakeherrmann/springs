from .agent_fibroblast import Agent_Fibroblast

import springs_interface as spr

class Lung:
	def __init__(self,
			spring_network=None,
			walls=None,
			agents=None,
			spring_break_variable=None,
			spring_break_threshold=0.0):
		self.net = spring_network
		self.walls = walls
		self.agents = agents
		self.spring_break_variable  = spring_break_variable
		self.spring_break_threshold = spring_break_threshold
		self.alveoli = None
		self.fibroblasts = None

	def save(self, save_dir):
		self.net.write_spring_network(save_dir)
		file_name = save_dir/'lung_walls.dat'
		file_format = spr.Structure.get_file_format()
		file_format.write_binary_file(file_name, self.walls)

	def load(self, load_dir):
		self.net.read_spring_network(load_dir)
		file_name = load_dir/'lung_walls.dat'
		file_format = spr.Structure.get_file_format()
		file_format.read_binary_file(file_name, self.walls)

	def stretch(self, stretch_increment, dimensions='all', boundary_indexes='all'):
		self.net.apply_stretch(stretch_increment, dimensions=dimensions, boundary_indexes=boundary_indexes)
		self.net.solve()
		self.net.calc_spring_force()
		if self.spring_break_variable is not None and self.spring_break_threshold is not None:
			self.net.break_spring(self.spring_break_variable, self.spring_break_threshold)

	def agent_actions(self, time_step):
		if self.agents is not None:
			for agent in self.agents:
				agent.do_actions(time_step)
			self.net.break_spring('stiffness_tension', 0.0, relop='<=')

	def add_fibroblast_every_spring(self):
		self.agents = []
		for spring in self.net.springs:
			self.agents.append( Agent_Fibroblast(location=spring) )

